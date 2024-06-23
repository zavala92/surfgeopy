import numpy as np
import pytest
from surfgeopy import pullback, subdivide, pushforward, quadrule_on_simplex, SimpleImplicitSurfaceProjection, integration

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
        mesh_path = "mesh_test/sphere_N=104.mat"
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
        mesh_path = "mesh_test/sphere_N=104.mat"
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
