import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.polynomial.chebyshev import chebfit, chebval
import scipy.io
from scipy import special
from surfgeopy import (
    DiagnosticIntegrationResult,
    DEFAULT_QUADRATURE_RULE,
    CompiledIntegrationScheme,
    ImplicitSurface,
    IndicatorRefinementResult,
    IntegrationConfig,
    IntegrationResult,
    LaplaceSingleLayerOperator,
    LevelSetSurface,
    MODEPY_GRUNDMANN_MOELLER,
    MODEPY_SIMPLEX_RULES,
    MODEPY_VIOREANU_ROKHLIN,
    MODEPY_XIAO_GIMBUTAS,
    PatchGeometry,
    ProjectionResult,
    ReferenceQuadrature,
    ScreenedLaplaceBeltramiParametrixOperator,
    ScreenedParametrixConfig,
    ScreenedParametrixResult,
    SurfaceGeometryResult,
    SurfaceMesh,
    SingularDiagnosticResult,
    SingularIntegralResult,
    SingularIntegrationConfig,
    affine_triangle_points,
    chebyshev_coefficients_2d,
    chebyshev_derivative_coefficients_2d,
    chebyshev_lobatto_nodes,
    chebyshev_values_2d,
    compile_integration,
    dense_chebyshev_coefficients_2d,
    integrate,
    integrate_with_diagnostics,
    integration,
    laplace_single_layer_potential,
    make_reference_quadrature,
    project_triangle_nodes,
    pullback,
    pushforward,
    quadrule_on_simplex,
    read_mesh_data,
    simplex_barycentric_coordinates,
    screened_laplace_beltrami_parametrix,
    surface_geometry,
    refine_by_indicator,
    subdivide,
    subdivide_conforming,
    SimpleImplicitSurfaceProjection,
)

MESH_PATH = Path(__file__).parent / "mesh_test" / "sphere_N=104.mat"

class TestSurfgeopyFunctions:

    def test_top_level_import_is_lazy(self):
        code = (
            "import sys\n"
            "import surfgeopy\n"
            "print(surfgeopy.__version__)\n"
            "print('numba' in sys.modules)\n"
            "print('minterpy' in sys.modules)\n"
            "print('modepy' in sys.modules)\n"
        )
        completed = subprocess.run(
            [sys.executable, "-c", code],
            cwd=Path(__file__).parents[1],
            capture_output=True,
            text=True,
            check=True,
        )

        version, has_numba, has_minterpy, has_modepy = completed.stdout.strip().splitlines()
        assert version == "1.0.0a0"
        assert has_numba == "False"
        assert has_minterpy == "False"
        assert has_modepy == "False"

    def test_module_cli_version(self):
        completed = subprocess.run(
            [sys.executable, "-m", "surfgeopy", "--version"],
            cwd=Path(__file__).parents[1],
            capture_output=True,
            text=True,
            check=True,
        )

        assert completed.stdout.strip() == "surfgeopy 1.0.0a0"

    def test_module_cli_doctor(self):
        completed = subprocess.run(
            [sys.executable, "-m", "surfgeopy", "doctor"],
            cwd=Path(__file__).parents[1],
            capture_output=True,
            text=True,
            check=True,
        )

        assert "surfgeopy 1.0.0a0" in completed.stdout
        assert "python " in completed.stdout
        assert "numpy:" in completed.stdout
        assert "numba:" in completed.stdout

    def test_module_cli_demo_smoke(self):
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "surfgeopy",
                "demo",
                "--mesh-refinement-level",
                "0",
                "--interpolation-degree",
                "2",
                "--integration-degree",
                "3",
            ],
            cwd=Path(__file__).parents[1],
            capture_output=True,
            text=True,
            check=True,
        )

        assert "surfgeopy sphere demo" in completed.stdout
        assert "area" in completed.stdout
        assert "quad points" in completed.stdout

    def test_compiled_integration_scheme_reuses_quadrature(self):
        surface = LevelSetSurface.unit_sphere(mesh_refinement_level=0)
        config = IntegrationConfig(
            interpolation_degree=2,
            integration_degree=3,
            quadrature_rule="Gauss_Legendre",
        )

        scheme = compile_integration(surface, config)
        direct = integrate(surface, lambda point: point[2] ** 2, config)
        compiled = scheme.integrate(lambda point: point[2] ** 2)

        assert isinstance(scheme, CompiledIntegrationScheme)
        assert scheme.n_faces == surface.mesh.n_faces
        assert scheme.n_quadrature_points == direct.n_quadrature_points
        np.testing.assert_allclose(compiled.values, direct.values)
        assert compiled.total == pytest.approx(direct.total)

        calls = {"count": 0}

        def vectorized_integrand(points):
            calls["count"] += 1
            return points[:, 2] ** 2

        vectorized = scheme.integrate(vectorized_integrand, vectorized=True)

        assert calls["count"] == 1
        np.testing.assert_allclose(vectorized.values, direct.values)
        assert vectorized.total == pytest.approx(direct.total)

        values_result = scheme.integrate_values(scheme.points[:, 2] ** 2)
        np.testing.assert_allclose(values_result.values, direct.values)

    def test_compiled_integration_scheme_validates_point_values(self):
        surface = LevelSetSurface.unit_sphere(mesh_refinement_level=0)
        config = IntegrationConfig(
            interpolation_degree=2,
            integration_degree=3,
            quadrature_rule="Gauss_Legendre",
        )
        scheme = compile_integration(surface, config)

        with pytest.raises(ValueError, match="point_values"):
            scheme.integrate_values(np.ones(scheme.n_quadrature_points + 1))

    def test_fft_chebyshev_coefficients_match_dense_solve(self):
        degree = 8
        nodes = chebyshev_lobatto_nodes(degree)
        uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
        values = np.empty((degree + 1, degree + 1, 3), dtype=float)
        values[..., 0] = uu**3 + 2.0 * uu * vv - 0.25 * vv**2
        values[..., 1] = np.cos(2.0 * uu) + vv**4
        values[..., 2] = uu**2 * vv - vv**3

        fft_coefficients = chebyshev_coefficients_2d(values)
        dense_coefficients = dense_chebyshev_coefficients_2d(values)

        np.testing.assert_allclose(fft_coefficients, dense_coefficients, atol=2.0e-12)

        samples_u = np.array([-0.75, -0.1, 0.35, 0.9])
        samples_v = np.array([-0.6, 0.2, 0.5, 0.8])
        fft_values = chebyshev_values_2d(fft_coefficients, samples_u, samples_v)
        dense_values = chebyshev_values_2d(dense_coefficients, samples_u, samples_v)
        np.testing.assert_allclose(fft_values, dense_values, atol=2.0e-12)

    def test_fft_chebyshev_derivatives_match_polynomial_derivatives(self):
        degree = 6
        nodes = chebyshev_lobatto_nodes(degree)
        uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
        values = uu**3 + 2.0 * uu * vv - 0.25 * vv**2
        coefficients = chebyshev_coefficients_2d(values)
        du_coefficients = chebyshev_derivative_coefficients_2d(coefficients, axis=0)
        dv_coefficients = chebyshev_derivative_coefficients_2d(coefficients, axis=1)

        samples_u = np.array([-0.8, -0.2, 0.3, 0.75])
        samples_v = np.array([-0.4, 0.1, 0.6, 0.9])
        du = chebyshev_values_2d(du_coefficients, samples_u, samples_v)
        dv = chebyshev_values_2d(dv_coefficients, samples_u, samples_v)

        np.testing.assert_allclose(du, 3.0 * samples_u**2 + 2.0 * samples_v, atol=2.0e-12)
        np.testing.assert_allclose(dv, 2.0 * samples_u - 0.5 * samples_v, atol=2.0e-12)

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

    def test_subdivide_conforming_splits_neighbor_hanging_edges(self):
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        faces = np.array([[0, 1, 2], [0, 2, 3]])

        refined_vertices, refined_faces = subdivide_conforming(vertices, faces, face_index=np.array([0]))

        assert refined_vertices.shape == (7, 3)
        assert refined_faces.shape == (6, 3)
        assert any(set(face) == {0, 5, 3} for face in refined_faces)
        assert any(set(face) == {5, 2, 3} for face in refined_faces)

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

    def test_project_triangle_nodes_requires_convergence(self):
        surface = ImplicitSurface(lambda _: 1.0, lambda _: np.zeros(3))
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        barycentric_points = np.array([[1.0, 0.0, 0.0]])

        with pytest.raises(RuntimeError, match="Projection failed to converge"):
            project_triangle_nodes(surface, vertices, barycentric_points)

    def test_reference_quadrature_factory(self):
        quadrature = make_reference_quadrature(4, "Pull_back_Gauss")
        assert isinstance(quadrature, ReferenceQuadrature)
        assert quadrature.size == 6
        assert quadrature.evaluation_points.shape == (6, 2)
        assert quadrature.weight_scale(0) > 0

    def test_default_quadrature_rule_is_vioreanu_rokhlin(self):
        config = IntegrationConfig(interpolation_degree=4)
        quadrature = make_reference_quadrature(4)

        assert DEFAULT_QUADRATURE_RULE == MODEPY_VIOREANU_ROKHLIN
        assert config.quadrature_rule == MODEPY_VIOREANU_ROKHLIN
        assert quadrature.rule == MODEPY_VIOREANU_ROKHLIN

    @pytest.mark.parametrize("rule", MODEPY_SIMPLEX_RULES)
    def test_modepy_simplex_reference_quadrature_factory(self, rule):
        quadrature = make_reference_quadrature(4, rule)

        assert isinstance(quadrature, ReferenceQuadrature)
        assert quadrature.reference_points.shape[1] == 2
        assert quadrature.evaluation_points.shape == quadrature.reference_points.shape
        np.testing.assert_allclose(np.sum(quadrature.weights), 0.5)
        assert np.all(quadrature.reference_points >= -1.0e-14)
        assert np.all(np.sum(quadrature.reference_points, axis=1) <= 1.0 + 1.0e-14)
        assert quadrature.weight_scale(0) > 0

    def test_modepy_quadrature_constants_are_public(self):
        assert MODEPY_XIAO_GIMBUTAS in MODEPY_SIMPLEX_RULES
        assert MODEPY_GRUNDMANN_MOELLER in MODEPY_SIMPLEX_RULES
        assert MODEPY_VIOREANU_ROKHLIN in MODEPY_SIMPLEX_RULES

    def test_vioreanu_rokhlin_unavailable_degree_suggests_xiao_gimbutas(self):
        with pytest.raises(ValueError, match="ModePy_XiaoGimbutas"):
            make_reference_quadrature(21, MODEPY_VIOREANU_ROKHLIN)

    def test_surface_mesh_from_mat(self):
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        assert mesh.n_vertices > 0
        assert mesh.n_faces > 0
        assert mesh.vertices.shape[1] == 3
        assert mesh.faces.shape[1] >= 3

    def test_surface_mesh_icosphere_constructor(self):
        mesh = SurfaceMesh.icosphere(refinement_level=1, radius=2.0)

        assert mesh.n_vertices == 42
        assert mesh.n_faces == 80
        np.testing.assert_allclose(np.linalg.norm(mesh.vertices, axis=1), 2.0)
        assert np.min(mesh.faces) >= 0
        assert np.max(mesh.faces) < mesh.n_vertices

        with pytest.raises(ValueError, match="refinement_level"):
            SurfaceMesh.icosphere(refinement_level=-1)
        with pytest.raises(ValueError, match="radius"):
            SurfaceMesh.icosphere(radius=0.0)

    def test_level_set_surface_sphere_constructor(self):
        surface = LevelSetSurface.sphere(
            radius=2.0,
            center=np.array([1.0, -1.0, 0.5]),
            mesh_refinement_level=0,
        )

        assert surface.mesh.n_faces == 20
        point = np.array([3.0, -1.0, 0.5])
        assert surface.level_set(point) == pytest.approx(0.0)
        np.testing.assert_allclose(surface.gradient(point), np.array([4.0, 0.0, 0.0]))
        np.testing.assert_allclose(
            np.linalg.norm(surface.mesh.vertices - np.array([1.0, -1.0, 0.5]), axis=1),
            2.0,
        )

        unit_surface = LevelSetSurface.unit_sphere(mesh_refinement_level=0)
        assert unit_surface.mesh.n_faces == 20
        assert unit_surface.level_set(np.array([1.0, 0.0, 0.0])) == pytest.approx(0.0)

    def test_read_mesh_data_uses_named_mat_arrays(self, tmp_path):
        mesh_path = tmp_path / "named_mesh.mat"
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        faces = np.array([[1, 2, 3]])
        scipy.io.savemat(
            mesh_path,
            {
                "metadata": np.array([[123.0]]),
                "faces": faces,
                "vertices": vertices,
            },
        )

        loaded_vertices, loaded_faces = read_mesh_data(str(mesh_path))

        np.testing.assert_allclose(loaded_vertices, vertices)
        np.testing.assert_array_equal(loaded_faces, np.array([[0, 1, 2]]))

    def test_read_mesh_data_supports_existing_xs_surfs_names(self):
        vertices, faces = read_mesh_data(str(MESH_PATH))

        assert vertices.shape == (54, 3)
        assert faces.shape == (104, 3)
        assert np.min(faces) == 0
        assert np.max(faces) == 53

    def test_read_mesh_data_rejects_invalid_face_data(self, tmp_path):
        mesh_path = tmp_path / "bad_mesh.mat"
        scipy.io.savemat(
            mesh_path,
            {
                "vertices": np.array([
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]),
                "faces": np.array([[1.0, 2.5, 3.0]]),
            },
        )

        with pytest.raises(ValueError, match="integer face indices"):
            read_mesh_data(str(mesh_path))

    def test_surface_mesh_from_gmsh_high_order_triangle(self, tmp_path):
        mesh_path = tmp_path / "quadratic_triangle.msh"
        mesh_path.write_text(
            "\n".join(
                [
                    "$MeshFormat",
                    "2.2 0 8",
                    "$EndMeshFormat",
                    "$Nodes",
                    "6",
                    "1 0 0 0",
                    "2 1 0 0",
                    "3 0 1 0",
                    "4 0.5 0 0",
                    "5 0.5 0.5 0",
                    "6 0 0.5 0",
                    "$EndNodes",
                    "$Elements",
                    "1",
                    "1 9 2 11 22 1 2 3 4 5 6",
                    "$EndElements",
                ]
            ),
            encoding="utf-8",
        )

        linear_mesh = SurfaceMesh.from_gmsh(str(mesh_path))
        high_order_mesh = SurfaceMesh.from_gmsh(str(mesh_path), preserve_order=True)

        assert linear_mesh.faces.shape == (1, 3)
        assert high_order_mesh.faces.shape == (1, 6)
        np.testing.assert_array_equal(linear_mesh.faces, np.array([[0, 1, 2]]))
        np.testing.assert_array_equal(high_order_mesh.faces, np.array([[0, 1, 2, 3, 4, 5]]))
        np.testing.assert_allclose(
            high_order_mesh.vertices[high_order_mesh.faces[0, 3:]],
            np.array([
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.0],
            ]),
        )

    def test_surface_mesh_to_gmsh_round_trips_high_order_triangle(self, tmp_path):
        mesh = SurfaceMesh(
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.0, 0.5, 0.0],
                ]
            ),
            np.array([[0, 1, 2, 3, 4, 5]]),
        )
        high_order_path = tmp_path / "high_order.msh"
        linear_path = tmp_path / "linearized.msh"

        mesh.to_gmsh(str(high_order_path))
        mesh.to_gmsh(str(linear_path), linearize=True)

        round_trip = SurfaceMesh.from_gmsh(str(high_order_path), preserve_order=True)
        linearized = SurfaceMesh.from_gmsh(str(linear_path), preserve_order=True)

        np.testing.assert_allclose(round_trip.vertices, mesh.vertices)
        np.testing.assert_array_equal(round_trip.faces, mesh.faces)
        np.testing.assert_array_equal(linearized.faces, np.array([[0, 1, 2]]))

    def test_surface_mesh_from_gmsh4_skips_non_triangular_blocks(self, tmp_path):
        mesh_path = tmp_path / "gmsh4_triangle.msh"
        mesh_path.write_text(
            "\n".join(
                [
                    "$MeshFormat",
                    "4.1 0 8",
                    "$EndMeshFormat",
                    "$Nodes",
                    "1 3 1 3",
                    "2 1 0 3",
                    "1",
                    "2",
                    "3",
                    "0 0 0",
                    "1 0 0",
                    "0 1 0",
                    "$EndNodes",
                    "$Elements",
                    "2 2 1 2",
                    "1 1 1 1",
                    "1 1 2",
                    "2 1 2 1",
                    "2 1 2 3",
                    "$EndElements",
                ]
            ),
            encoding="utf-8",
        )

        mesh = SurfaceMesh.from_gmsh(str(mesh_path))

        np.testing.assert_allclose(
            mesh.vertices,
            np.array([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]),
        )
        np.testing.assert_array_equal(mesh.faces, np.array([[0, 1, 2]]))

    def test_surface_mesh_gmsh_validates_triangle_element_types(self, tmp_path):
        mesh_path = tmp_path / "linear_triangle.msh"
        SurfaceMesh(
            np.array([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]),
            np.array([[0, 1, 2]]),
        ).to_gmsh(str(mesh_path))

        with pytest.raises(ValueError, match="Unsupported triangular Gmsh element type"):
            SurfaceMesh.from_gmsh(str(mesh_path), element_types=[3])

        unsupported = SurfaceMesh(
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            np.array([[0, 1, 2, 3]]),
        )
        with pytest.raises(ValueError, match="Supported node counts"):
            unsupported.to_gmsh(str(tmp_path / "quad.msh"))

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

    def test_refine_by_indicator_marks_large_face_center_values(self):
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        mesh = SurfaceMesh(vertices, faces)
        surface = LevelSetSurface(mesh, lambda x: x[2], lambda _: np.array([0.0, 0.0, 1.0]))

        refined = refine_by_indicator(
            surface,
            lambda x: x[0],
            max_iterations=2,
            threshold_fraction=0.5,
            use_absolute=False,
        )

        assert isinstance(refined, IndicatorRefinementResult)
        assert refined.n_iterations == 2
        assert refined.initial_n_faces == 2
        assert refined.initial_n_vertices == 4
        assert refined.n_vertices == 7
        assert refined.n_refined_iterations == 1
        assert refined.total_marked_faces == 5
        assert refined.face_growth_factor == pytest.approx(3.0)
        assert refined.vertex_growth_factor == pytest.approx(7.0 / 4.0)
        assert refined.max_marked_fraction == pytest.approx(2.0 / 3.0)
        assert refined.stop_reason == "max_iterations"
        np.testing.assert_array_equal(refined.face_counts, np.array([2, 6, 6]))
        np.testing.assert_array_equal(refined.vertex_counts, np.array([4, 7, 7]))
        np.testing.assert_allclose(refined.max_indicators, np.array([2.0 / 3.0, 5.0 / 6.0]))
        np.testing.assert_allclose(refined.marked_fractions, np.array([0.5, 2.0 / 3.0]))
        assert refined.history[0].n_faces == 2
        assert refined.history[0].n_faces_after == 6
        assert refined.history[0].n_vertices == 4
        assert refined.history[0].n_vertices_after == 7
        assert refined.history[0].n_new_faces == 4
        assert refined.history[0].n_new_vertices == 3
        assert refined.history[0].n_marked_faces == 1
        assert refined.history[0].n_unmarked_faces == 1
        assert refined.history[0].did_refine
        np.testing.assert_array_equal(refined.history[0].marked_faces, np.array([0]))
        np.testing.assert_allclose(refined.history[0].indicator_values, np.array([2.0 / 3.0, 1.0 / 3.0]))
        assert refined.history[0].threshold == pytest.approx(1.0 / 3.0)
        assert refined.history[0].min_indicator == pytest.approx(1.0 / 3.0)
        assert refined.history[0].mean_indicator == pytest.approx(0.5)
        assert refined.history[0].median_indicator == pytest.approx(0.5)
        assert refined.history[0].std_indicator == pytest.approx(1.0 / 6.0)
        assert refined.history[-1].stop_reason == "max_iterations"
        assert refined.history[-1].n_marked_faces == 4
        assert not refined.history[-1].did_refine
        assert refined.n_faces == 6
        assert "Stop reason" in refined.summary()
        assert "Final faces" in refined.summary()
        assert "Last mean indicator" in refined.summary()

    def test_refine_by_indicator_reports_zero_indicator_stop(self):
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        surface = LevelSetSurface(
            SurfaceMesh(vertices, np.array([[0, 1, 2]])),
            lambda x: x[2],
            lambda _: np.array([0.0, 0.0, 1.0]),
        )

        refined = refine_by_indicator(surface, lambda _: 0.0, max_iterations=3)

        assert refined.n_iterations == 1
        assert refined.n_faces == 1
        assert refined.n_vertices == 3
        assert refined.stop_reason == "zero_indicator"
        assert refined.history[0].n_marked_faces == 0
        assert refined.history[0].marked_fraction == 0.0
        assert refined.history[0].n_new_faces == 0
        assert refined.history[0].n_new_vertices == 0
        assert not refined.history[0].did_refine
        np.testing.assert_allclose(refined.mean_indicators, np.array([0.0]))

    def test_refine_by_indicator_validates_parameters(self):
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        surface = LevelSetSurface(
            SurfaceMesh(vertices, np.array([[0, 1, 2]])),
            lambda x: x[2],
            lambda _: np.array([0.0, 0.0, 1.0]),
        )

        with pytest.raises(ValueError, match="max_iterations"):
            refine_by_indicator(surface, lambda _: 1.0, max_iterations=0)
        with pytest.raises(ValueError, match="threshold_fraction"):
            refine_by_indicator(surface, lambda _: 1.0, threshold_fraction=0.0)
        with pytest.raises(ValueError, match="non-finite"):
            refine_by_indicator(surface, lambda _: np.nan)

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
        assert geometry.n_patches == mesh.n_faces
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

        patch = geometry.patch(0)
        assert isinstance(patch, PatchGeometry)
        assert patch.index == 0
        assert patch.n_points == geometry.offsets[1] - geometry.offsets[0]
        np.testing.assert_allclose(patch.points, geometry.points[:patch.n_points])
        assert patch.area > 0.0
        assert patch.radius > 0.0
        assert patch.center.shape == (3,)
        assert len(geometry.patches) == geometry.n_patches
        assert geometry.patch_areas.shape == (geometry.n_patches,)
        assert geometry.patch_centers.shape == (geometry.n_patches, 3)
        assert geometry.patch_radii.shape == (geometry.n_patches,)
        np.testing.assert_allclose(np.sum(geometry.patch_areas), np.sum(geometry.weights))
        assert np.all(geometry.patch_radii > 0.0)
        with pytest.raises(IndexError, match="patch index"):
            geometry.patch(geometry.n_patches)

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

    def test_legacy_integration_uses_default_vioreanu_rokhlin(self):
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

    def test_laplace_single_layer_prototype_on_unit_sphere(self):
        zero_levelset_function = lambda x: x[0]**2 + x[1]**2 + x[2]**2 - 1
        gradient_function = lambda x: np.array([2*x[0], 2*x[1], 2*x[2]])
        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = SingularIntegrationConfig(
            interpolation_degree=4,
            regular_order=8,
            smooth_degree=4,
            moment_order=18,
            correction_order=10,
            near_threshold=1.5,
        )
        assert config.singular_model == "curvature"
        with pytest.raises(ValueError, match="singular_model"):
            SingularIntegrationConfig(singular_model="linear")
        with pytest.raises(ValueError, match="evaluation_strategy"):
            SingularIntegrationConfig(evaluation_strategy="qbx")

        target = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
        result = laplace_single_layer_potential(surface, target, config=config)

        assert isinstance(result, SingularIntegralResult)
        assert result.values.shape == (1,)
        assert result.qbx_target_flags is None
        assert result.qbx_radii is None
        assert result.qbx_orders is None
        assert result.qbx_convergence_ratios is None
        assert result.qbx_error_estimates is None
        assert result.near_panel_counts[0] > 0
        assert result.singular_panel_counts[0] > 0
        assert abs(result.values[0] - 1.0) < 0.15

        near_target = 1.001 * target
        near_result = LaplaceSingleLayerOperator(surface, config).evaluate(near_target)

        assert near_result.near_panel_counts[0] > 0
        assert near_result.singular_panel_counts[0] == 0
        assert abs(near_result.values[0] - 1.0 / 1.001) < 1.0e-3

        hybrid_config = SingularIntegrationConfig(
            interpolation_degree=4,
            regular_order=8,
            smooth_degree=4,
            moment_order=18,
            correction_order=10,
            near_threshold=1.5,
            evaluation_strategy="hybrid_qbx",
            qbx_order=8,
            qbx_quadrature_order=20,
            qbx_radius_factor=0.5,
        )
        hybrid_result = LaplaceSingleLayerOperator(surface, hybrid_config).evaluate(target)

        assert hybrid_result.qbx_target_flags.shape == (1,)
        assert hybrid_result.qbx_target_flags[0]
        assert hybrid_result.qbx_radii.shape == (1,)
        assert hybrid_result.qbx_orders.shape == (1,)
        assert hybrid_result.qbx_convergence_ratios.shape == (1,)
        assert hybrid_result.qbx_error_estimates.shape == (1,)
        assert hybrid_result.qbx_radii[0] > 0.0
        assert hybrid_result.qbx_orders[0] == hybrid_config.qbx_order
        assert hybrid_result.qbx_convergence_ratios[0] >= 0.0
        assert hybrid_result.qbx_error_estimates[0] >= 0.0
        assert hybrid_result.near_panel_counts[0] > 0
        assert hybrid_result.singular_panel_counts[0] > 0
        assert abs(hybrid_result.values[0] - 1.0) < 5.0e-5

        operator = LaplaceSingleLayerOperator(surface, config)
        operator_result = operator.evaluate(target)

        assert operator.n_patches == mesh.n_faces
        np.testing.assert_allclose(operator_result.values, result.values)

        diagnostic = operator.evaluate_with_diagnostics(target, target_tolerance=1.0e-8)

        assert isinstance(diagnostic, SingularDiagnosticResult)
        assert diagnostic.values.shape == (1,)
        assert diagnostic.base_values.shape == (1,)
        assert diagnostic.enriched_values.shape == (1,)
        assert diagnostic.tail_indicators.shape == (1,)
        assert diagnostic.max_tail_indicators.shape == (1,)
        assert diagnostic.absolute_error_estimates[0] >= 0.0
        assert diagnostic.tail_indicators[0] >= 0.0
        assert diagnostic.max_tail_indicators[0] >= 0.0
        assert diagnostic.tail_indicators[0] >= diagnostic.max_tail_indicators[0]
        assert diagnostic.near_panel_counts[0] > 0
        assert diagnostic.singular_panel_counts[0] > 0
        assert diagnostic.corrected_panel_counts[0] > 0
        assert diagnostic.recommended_config.interpolation_degree > config.interpolation_degree
        assert diagnostic.n_base_patches == mesh.n_faces
        assert diagnostic.n_enriched_patches_built <= diagnostic.n_base_patches
        assert 0.0 < diagnostic.enriched_patch_fraction <= 1.0

    def test_screened_laplace_beltrami_parametrix_on_unit_sphere(self):
        alpha = 0.1
        degree = (-1.0 + np.sqrt(1.0 - 4.0 * alpha)) / 2.0

        def exact_green(cosine):
            return -special.lpmv(0, degree, -cosine) / (4.0 * np.sin(np.pi * degree))

        def singular_green(cosine):
            chord = np.sqrt(np.maximum(2.0 * (1.0 - cosine), np.finfo(float).tiny))
            distance = chord * (1.0 + chord**2 / 24.0 + 3.0 * chord**4 / 640.0)
            return special.k0(np.sqrt(alpha) * distance) / (2.0 * np.pi)

        cosine_samples = np.cos(np.linspace(np.pi, 1.0e-4, 1000))
        remainder_coefficients = chebfit(
            cosine_samples,
            exact_green(cosine_samples) - singular_green(cosine_samples),
            12,
        )

        def smooth_remainder(source, target):
            cosine = np.clip(np.dot(source, target), -1.0, 1.0)
            return float(chebval(cosine, remainder_coefficients))

        def zero_levelset_function(x):
            return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1

        def gradient_function(x):
            return np.array([2 * x[0], 2 * x[1], 2 * x[2]])

        def density(point):
            z = point[2]
            y1 = z
            y2 = 0.5 * (3.0 * z**2 - 1.0)
            return y1 + 0.25 * y2

        mesh = SurfaceMesh.from_mat(str(MESH_PATH))
        surface = LevelSetSurface(mesh, zero_levelset_function, gradient_function)
        config = ScreenedParametrixConfig(
            interpolation_degree=3,
            regular_order=8,
            smooth_degree=4,
            moment_order=14,
            near_threshold=1.5,
        )
        with pytest.raises(ValueError, match="alpha"):
            ScreenedLaplaceBeltramiParametrixOperator(surface, 0.0, config)
        with pytest.raises(ValueError, match="singular_model"):
            ScreenedParametrixConfig(singular_model="linear")

        target = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
        result = screened_laplace_beltrami_parametrix(
            surface,
            target,
            alpha,
            density,
            smooth_remainder,
            config,
        )

        z = target[2]
        exact_value = z / (alpha + 2.0) + 0.25 * (0.5 * (3.0 * z**2 - 1.0)) / (
            alpha + 6.0
        )

        assert isinstance(result, ScreenedParametrixResult)
        assert result.values.shape == (1,)
        assert result.near_panel_counts[0] > 0
        assert result.singular_panel_counts[0] > 0
        assert result.tail_indicators[0] >= 0.0
        assert abs(result.values[0] - exact_value) < 2.0e-2

        operator = ScreenedLaplaceBeltramiParametrixOperator(
            surface,
            alpha,
            config,
            smooth_remainder,
        )
        assert operator.n_patches == mesh.n_faces
        np.testing.assert_allclose(operator.evaluate(target, density).values, result.values)

    def test_screened_laplace_beltrami_convergence_benchmark_smoke(self):
        module_path = (
            Path(__file__).parents[1]
            / "examples"
            / "pde"
            / "screened_laplace_beltrami_sphere_convergence.py"
        )
        spec = importlib.util.spec_from_file_location(
            "screened_laplace_beltrami_sphere_convergence",
            module_path,
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)

        points = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
        np.testing.assert_allclose(
            module.manufactured_rhs(points, 0.1),
            np.array([2.1, 0.0]),
        )

        geometry_rows = module.run_geometry_study(
            MESH_PATH,
            interpolation_degree=3,
            integration_degree=3,
            refinement_levels=(0, 1),
        )

        assert len(geometry_rows) == 2
        assert np.isfinite(geometry_rows[0].l2_error)
        assert np.isfinite(geometry_rows[1].l2_error)
        assert geometry_rows[1].l2_error < geometry_rows[0].l2_error

        remainder_rows = module.run_remainder_study(degrees=(2, 4))

        assert len(remainder_rows) == 2
        assert np.isfinite(remainder_rows[0].max_error)
        assert remainder_rows[1].max_error < remainder_rows[0].max_error

if __name__ == '__main__':
    pytest.main()
