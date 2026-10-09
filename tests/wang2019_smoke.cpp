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
    // This reproduces the CGAL Euler delta semantics independently from
    // angle-policy acceptance. One split creates +1 vertex / +2 triangles.
    {
        Mesh probe(mesh);
        const auto h = halfedge(*edges(probe).first, probe);
        const auto f0 = face(h, probe);
        const auto f1 = face(opposite(h, probe), probe);
        const auto oldV = probe.number_of_vertices(), oldF = probe.number_of_faces();
        const Point midpoint = CGAL::midpoint(
            probe.point(source(h, probe)), probe.point(target(h, probe)));
        const auto new_h = CGAL::Euler::split_edge(h, probe);
        const auto vm = target(new_h, probe);
        probe.point(vm) = midpoint;
        if (!split_quad(probe, f0, vm) || !split_quad(probe, f1, vm))
            throw std::runtime_error("Quad retriangulation failed");
        verify_delta(probe, oldV, oldF, +1, +2, "smoke split");

        bool collapsed = false;
        for (const auto e : edges(probe)) {
            if (is_border(e, probe) ||
                !CGAL::Euler::does_satisfy_link_condition(e, probe))
                continue;
            const auto vc = probe.number_of_vertices(), fc = probe.number_of_faces();
            CGAL::Euler::collapse_edge(e, probe);
            verify_delta(probe, vc, fc, -1, -2, "smoke collapse");
            collapsed = true;
            break;
        }
        if (!collapsed) throw std::runtime_error("No legal test collapse");
    }
    Options o;
    o.rounds = 3;
    o.smoothing_steps = 1;
    o.budget_fraction = 0.5;
    o.protect_features = false;
    o.verbose = false;
    const Statistics stats = run(mesh, o);
    if (!mesh.is_valid() || !CGAL::is_triangle_mesh(mesh) ||
        stats.vertices_before != stats.vertices_after ||
        stats.faces_before != stats.faces_after ||
        stats.insertions != stats.collapses ||
        stats.vertices_after == 0 || stats.faces_after == 0)
        throw std::runtime_error("Output mesh invalid or empty");
    std::cout << "Wang2019 smoke PASS: " << stats.vertices_after
              << " vertices, " << stats.faces_after << " triangles\n";
}
