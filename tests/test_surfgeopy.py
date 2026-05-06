import numpy as np
import pytest
from pathlib import Path
from surfgeopy import (
    AdaptiveIntegrationResult,
    DiagnosticIntegrationResult,
    ImplicitSurface,
    IntegrationConfig,
    IntegrationResult,
    LevelSetSurface,
    RECURSIVE_NODES_GAUSS_LEGENDRE,
    ProjectionResult,
    ReferenceQuadrature,
    SurfaceGeometryResult,
    SurfaceMesh,
    affine_triangle_points,
    adaptive_integrate,
    integrate,
    integrate_with_diagnostics,
    integration,
    make_reference_quadrature,
    pullback,
    pushforward,
    quadrule_on_simplex,
    simplex_barycentric_coordinates,
    surface_geometry,
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

    def test_integrate_with_diagnostics_api(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(
            interpolation_degree=4,
            refinement_level=0,
            integration_degree=8,
            quadrature_rule="Gauss_Legendre",
        )

        diagnostics = integrate_with_diagnostics(
            surface,
            lambda _: 1.0,
            config,
            interpolation_degree_step=2,
            integration_degree_step=2,
            relative_tolerance=1.0e-6,
        )

        assert isinstance(diagnostics, DiagnosticIntegrationResult)
        assert isinstance(diagnostics.base_result, IntegrationResult)
        assert isinstance(diagnostics.enriched_result, IntegrationResult)
        assert diagnostics.enriched_result.config.interpolation_degree == 6
        assert diagnostics.enriched_result.config.integration_degree == 10
        assert diagnostics.total == diagnostics.enriched_result.total
        assert diagnostics.n_quadrature_points == diagnostics.enriched_result.n_quadrature_points
        assert diagnostics.absolute_error_estimate >= 0.0
        assert diagnostics.relative_error_estimate >= 0.0
        assert diagnostics.local_absolute_errors.shape == (mesh.n_faces,)
        assert diagnostics.max_local_error >= 0.0
        assert diagnostics.target_reached in {True, False}
        assert "Estimated abs. error" in diagnostics.summary()
        assert np.abs(4 * np.pi - diagnostics.total) < 1e-6

    def test_integrate_with_diagnostics_validates_enrichment(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(interpolation_degree=4)

        with pytest.raises(ValueError, match="interpolation_degree_step"):
            integrate_with_diagnostics(surface, lambda _: 1.0, config, interpolation_degree_step=-1)
        with pytest.raises(ValueError, match="integration_degree_step"):
            integrate_with_diagnostics(surface, lambda _: 1.0, config, integration_degree_step=-1)
        with pytest.raises(ValueError, match="at least one"):
            integrate_with_diagnostics(
                surface,
                lambda _: 1.0,
                config,
                interpolation_degree_step=0,
                integration_degree_step=0,
            )

    def test_adaptive_integrate_refines_marked_faces(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(
            interpolation_degree=4,
            refinement_level=0,
            integration_degree=8,
            quadrature_rule="Gauss_Legendre",
        )

        adaptive = adaptive_integrate(
            surface,
            lambda _: 1.0,
            config,
            absolute_tolerance=0.0,
            relative_tolerance=0.0,
            max_iterations=2,
            marking_fraction=0.1,
        )

        assert isinstance(adaptive, AdaptiveIntegrationResult)
        assert adaptive.n_iterations == 2
        assert adaptive.history[0].n_faces == mesh.n_faces
        assert adaptive.history[0].n_marked_faces >= 1
        assert adaptive.n_faces > mesh.n_faces
        assert adaptive.absolute_error_estimate >= 0.0
        assert adaptive.relative_error_estimate >= 0.0
        assert "Final faces" in adaptive.summary()
        assert np.abs(4 * np.pi - adaptive.total) < 1e-5

    def test_adaptive_integrate_validates_parameters(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(interpolation_degree=4)

        with pytest.raises(ValueError, match="max_iterations"):
            adaptive_integrate(surface, lambda _: 1.0, config, max_iterations=0)
        with pytest.raises(ValueError, match="marking_fraction"):
            adaptive_integrate(surface, lambda _: 1.0, config, marking_fraction=0.0)
        with pytest.raises(ValueError, match="min_marked_faces"):
            adaptive_integrate(surface, lambda _: 1.0, config, min_marked_faces=0)

    def test_surface_geometry_api_on_sphere(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(
            interpolation_degree=6,
            refinement_level=0,
            integration_degree=8,
            quadrature_rule="Gauss_Legendre",
        )

        geometry = surface_geometry(surface, config)

        assert isinstance(geometry, SurfaceGeometryResult)
        assert geometry.n_faces == mesh.n_faces
        assert geometry.n_points == geometry.points.shape[0]
        assert geometry.weights.shape == (geometry.n_points,)
        assert geometry.tangent_u.shape == geometry.points.shape
        assert geometry.tangent_v.shape == geometry.points.shape
        assert geometry.normal.shape == geometry.points.shape
        assert geometry.metric_tensor.shape == (geometry.n_points, 2, 2)
        assert geometry.second_fundamental_form.shape == (geometry.n_points, 2, 2)
        assert geometry.area_density.shape == (geometry.n_points,)
        assert geometry.mean_curvature.shape == (geometry.n_points,)
        assert geometry.gaussian_curvature.shape == (geometry.n_points,)
        assert np.all(geometry.area_density > 0.0)
        np.testing.assert_allclose(np.linalg.norm(geometry.points, axis=1), 1.0, atol=1e-10)
        np.testing.assert_allclose(np.linalg.norm(geometry.normal, axis=1), 1.0, atol=1e-12)
        assert np.all(np.isfinite(geometry.mean_curvature))
        assert np.all(np.isfinite(geometry.gaussian_curvature))
        assert np.abs(np.median(geometry.gaussian_curvature) - 1.0) < 1e-3
        assert np.abs(np.median(np.abs(geometry.mean_curvature)) - 1.0) < 1e-3

    def test_surface_geometry_requires_second_order_interpolation(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = IntegrationConfig(interpolation_degree=1)

        with pytest.raises(ValueError, match="at least 2"):
            surface_geometry(surface, config)

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
