import numpy as np
import pytest
from pathlib import Path
from surfgeopy import (
    ImplicitSurface,
    IntegrationConfig,
    IntegrationResult,
    LevelSetSurface,
    RECURSIVE_NODES_GAUSS_LEGENDRE,
    ProjectionResult,
    ReferenceQuadrature,
    SurfaceMesh,
    affine_triangle_points,
    integrate,
    integration,
    make_reference_quadrature,
    pullback,
    pushforward,
    quadrule_on_simplex,
    simplex_barycentric_coordinates,
    subdivide,
    SimpleImplicitSurfaceProjection,
)

MESH_PATH = Path(__file__).parent / "mesh_test" / "sphere_N=104.mat"

class TestSurfgeopyFunctions:

    def test_pullback(self):
        unisolvent_nodes_triangle = np.array([[0.0, 0.0], [0.0, 1.0], [0.5, 0.5], [1.0, 0.0]])
        result_square = pullback(unisolvent_nodes_triangle)
        expected_result_square = np.array([[-1, -1], [-1, 1], [1, 1], [1, -1]])
        assert np.array_equal(result_square, expected_result_square)

    def test_subdivide(self):
        vertices = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        faces = np.array([[0, 1, 2]])
        expected_vertices = np.array([
            [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
            [0.5, 0.5, 0.], [0.5, 0., 0.5], [0., 0.5, 0.5]
        ])
        expected_faces = np.array([
            [0, 3, 4], [3, 1, 5], [4, 5, 2], [3, 5, 4]
        ])
        result = subdivide(vertices, faces)
        np.testing.assert_array_equal(result[0], expected_vertices)
        np.testing.assert_array_equal(result[1], expected_faces)

    def test_pushforward(self):
        unisolvent_nodes_square = np.array([[-1, -1], [-1, 1], [1, 1], [1, -1]])
        result_triangle = pushforward(unisolvent_nodes_square)
        expected_result_triangle = np.array([[0.0, 0.0], [0.0, 1.0], [0.5, 0.5], [1.0, 0.0]])
        assert np.array_equal(result_triangle, expected_result_triangle)

    def test_simplex_barycentric_coordinates(self):
        points = np.array([[0.0, 0.0], [0.25, 0.5], [1.0, 0.0]])
        expected = np.array([[1.0, 0.0, 0.0], [0.25, 0.25, 0.5], [0.0, 1.0, 0.0]])
        np.testing.assert_allclose(simplex_barycentric_coordinates(points), expected)

    def test_affine_triangle_points(self):
        vertices = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        barycentric_points = np.array([[0.25, 0.25, 0.5]])
        expected = np.array([[0.5, 2.0, 0.0]])
        np.testing.assert_allclose(affine_triangle_points(vertices, barycentric_points), expected)

    def test_implicit_surface_project_does_not_mutate_input(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        surface = ImplicitSurface(zero_levelset_function, gradient_function)
        point = np.array([0.5, 0.5, 0.5])
        projected = surface.project(point)

        np.testing.assert_allclose(point, np.array([0.5, 0.5, 0.5]))
        assert np.abs(zero_levelset_function(projected)) < 1e-16

    def test_implicit_surface_project_with_info(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        surface = ImplicitSurface(zero_levelset_function, gradient_function)
        result = surface.project_with_info(np.array([0.5, 0.5, 0.5]))

        assert isinstance(result, ProjectionResult)
        assert result.converged
        assert result.iterations > 0
        assert result.residual < 1e-14

    def test_implicit_surface_projection_reports_failure(self):
        surface = ImplicitSurface(lambda _: 1.0, lambda _: np.zeros(3))
        result = surface.project_with_info(np.array([0.0, 0.0, 0.0]))

        assert not result.converged
        assert result.residual == 1.0

    def test_reference_quadrature_factory(self):
        quadrature = make_reference_quadrature(4, "Pull_back_Gauss")
        assert isinstance(quadrature, ReferenceQuadrature)
        assert quadrature.size == 6
        assert quadrature.evaluation_points.shape == (6, 2)
        assert quadrature.weight_scale(0) > 0

    def test_recursivenodes_reference_quadrature_factory(self):
        quadrature = make_reference_quadrature(3, RECURSIVE_NODES_GAUSS_LEGENDRE)

        assert isinstance(quadrature, ReferenceQuadrature)
        assert quadrature.size == 9
        assert quadrature.reference_points.shape == (9, 2)
        np.testing.assert_allclose(np.sum(quadrature.weights), 0.5)
        assert np.all(quadrature.reference_points >= 0.0)
        assert np.all(np.sum(quadrature.reference_points, axis=1) <= 1.0)

    def test_surface_mesh_from_mat(self):
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        assert mesh.n_vertices > 0
        assert mesh.n_faces > 0
        assert mesh.vertices.shape[1] == 3
        assert mesh.faces.shape[1] >= 3

    def test_high_level_integrate_api(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(
            interpolation_degree=6,
            refinement_level=1,
            integration_degree=14,
            quadrature_rule="Gauss_Legendre",
        )

        result = integrate(surface, lambda _: 1.0, config)

        assert isinstance(result, IntegrationResult)
        assert result.values.shape == (mesh.n_faces,)
        assert result.points.shape[1] == 3
        assert result.weights.shape == result.points.shape[:1]
        assert result.n_quadrature_points == result.points.shape[0]
        assert np.abs(4 * np.pi - result.total) < 1e-10
        assert float(result) == result.total

    def test_integration_config_validates_degrees(self):
        with pytest.raises(ValueError, match="interpolation_degree"):
            IntegrationConfig(interpolation_degree=0)
        with pytest.raises(ValueError, match="refinement_level"):
            IntegrationConfig(interpolation_degree=1, refinement_level=-1)
        with pytest.raises(ValueError, match="integration_degree"):
            IntegrationConfig(interpolation_degree=1, integration_degree=0)

    @pytest.mark.parametrize("deg, expected_num_weights", [
        (1, 1), (2, 3), (3, 4), (4, 6), (5, 7), (6, 11), (7, 12), (8, 16), (9, 19),
        (10, 24), (11, 27), (12, 32), (13, 36), (14, 42), (15, 46), (16, 52), (17, 57),
        (18, 66), (19, 70), (20, 78), (21, 85), (22, 93), (23, 100), (24, 109), (25, 117)
    ])
    def test_quadrule_on_simplex(self, deg, expected_num_weights):
        weights, _ = quadrule_on_simplex(deg)
        assert len(weights) == expected_num_weights

    def test_closest_point(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        x0 = np.array([0.5, 0.5, 0.5])
        pnts_p = SimpleImplicitSurfaceProjection(zero_levelset_function, gradient_function, x0)
        assert np.abs(zero_levelset_function(pnts_p) < 1e-16)

    def test_integration_pull_back_gauss(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh_path = str(MESH_PATH)
        interpolation_degree = 8
        lp_degree = float("inf")
        refinement_level = 1
        constant_function = lambda _: 1
        result = integration(zero_levelset_function, gradient_function, mesh_path, interpolation_degree, lp_degree,
                             refinement_level, constant_function)
        assert np.abs(4 * np.pi - np.sum(result)) < 1e-12 

    def test_integration_gauss_legendre(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh_path = str(MESH_PATH)
        interpolation_degree = 6
        lp_degree = float("inf")
        refinement_level = 1
        integration_degree = 14
        quadrature_rule = "Gauss_Legendre"
        constant_function = lambda _: 1
        result = integration(zero_levelset_function, gradient_function, mesh_path, interpolation_degree, lp_degree,
                             refinement_level, constant_function, deg_integration=integration_degree,
                             quadrature_rule=quadrature_rule)
        assert np.abs(4 * np.pi - np.sum(result)) < 1e-10

if __name__ == '__main__':
    pytest.main()
