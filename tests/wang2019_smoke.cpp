#include "rar/Wang2019Remesher.h"
#include <CGAL/boost/graph/helpers.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace rar;
    using namespace rar::wang2019;
    const Point a(0, 0, 0), b(1, 0, 0), c(0.5, std::sqrt(3.0)/2, 0);
    const Angles equilateral = triangle_angles(a, b, c);
    if (std::abs(equilateral.min - 60) > 1e-8 ||
        std::abs(equilateral.max - 60) > 1e-8)
        throw std::runtime_error("Angle calculations are incorrect");

    // A closed 2-manifold tetrahedron, no dependence on external model files.
    Mesh mesh;
    const auto v0 = mesh.add_vertex(Point(0, 0, 0));
    const auto v1 = mesh.add_vertex(Point(3, 0, 0));
    const auto v2 = mesh.add_vertex(Point(0, 1, 0));
    const auto v3 = mesh.add_vertex(Point(0, 0, 1));
    if (mesh.add_face(v0, v2, v1) == Mesh::null_face() ||
        mesh.add_face(v0, v1, v3) == Mesh::null_face() ||
        mesh.add_face(v1, v2, v3) == Mesh::null_face() ||
        mesh.add_face(v2, v0, v3) == Mesh::null_face())
        throw std::runtime_error("Could not construct tetrahedron");
    Options o;
    o.rounds = 3;
    o.smoothing_steps = 1;
    o.budget_fraction = 0.5;
    o.protect_features = false;
    o.verbose = false;
    const Statistics stats = run(mesh, o);
    if (!mesh.is_valid() || !CGAL::is_triangle_mesh(mesh) ||
        stats.vertices_after == 0 || stats.faces_after == 0)
        throw std::runtime_error("Output mesh invalid or empty");
    std::cout << "Wang2019 smoke PASS: " << stats.vertices_after
              << " vertices, " << stats.faces_after << " triangles\n";
}
