#include <CGAL/Polygon_mesh_processing/remesh.h>

#include "rar/FeatureConstraints.h"
#include "rar/OutputNaming.h"
#include "rar/PropertyMapSizingField.h"
#include "rar/RemeshBackend.h"

#include <CGAL/boost/graph/generators.h>
#include <CGAL/squared_distance_3.h>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using rar::Mesh;
using rar::Point;

void require(bool condition, const char* message)
{
    // Unlike assert(), these checks remain active in Release builds.
    if (!condition) throw std::runtime_error(message);
}

Mesh open_grid()
{
    Mesh mesh;
    const std::array<double, 5> coordinates{0., 0.08, 0.3, 1., 2.};
    std::array<std::array<Mesh::Vertex_index, 5>, 5> v;
    for (std::size_t y = 0; y < 5; ++y) {
        for (std::size_t x = 0; x < 5; ++x) {
            const double px = coordinates[x], py = coordinates[y];
            v[y][x] = mesh.add_vertex(Point(px, py, 0.12 * px * py));
        }
    }
    for (std::size_t y = 0; y < 4; ++y) {
        for (std::size_t x = 0; x < 4; ++x) {
            require(mesh.add_face(v[y][x], v[y][x+1], v[y+1][x+1]) != Mesh::null_face(), "grid face");
            require(mesh.add_face(v[y][x], v[y+1][x+1], v[y+1][x]) != Mesh::null_face(), "grid face");
        }
    }
    return mesh;
}

Mesh sharp_cube()
{
    Mesh mesh;
    std::vector<Mesh::Vertex_index> v;
    for (const Point& p : {Point(0,0,0), Point(1,0,0), Point(1,1,0), Point(0,1,0),
                          Point(0,0,1), Point(1,0,1), Point(1,1,1), Point(0,1,1)}) {
        v.push_back(mesh.add_vertex(p));
    }
    const int triangles[][3] = {
        {0,2,1}, {0,3,2}, {4,5,6}, {4,6,7}, {0,1,5}, {0,5,4},
        {1,2,6}, {1,6,5}, {2,3,7}, {2,7,6}, {3,0,4}, {3,4,7}
    };
    for (const auto& t : triangles) {
        require(mesh.add_face(v[t[0]], v[t[1]], v[t[2]]) != Mesh::null_face(), "cube face");
    }
    return mesh;
}

Mesh closed_surface()
{
    Mesh mesh;
    CGAL::make_icosahedron(mesh, Point(0,0,0), 1.);
    for (const auto v : mesh.vertices()) {
        const Point p = mesh.point(v);
        mesh.point(v) = Point(1.3 * p.x(), p.y(), 0.8 * p.z());
    }
    return mesh;
}

void run(Mesh& mesh, rar::RemeshBackend backend, bool features,
         bool project, bool uniform, unsigned int iterations, unsigned int relax)
{
    auto constraints = rar::make_feature_constraints(mesh, features, 50.);
    std::vector<std::pair<Mesh::Vertex_index, Point>> fixed;
    for (const auto v : mesh.vertices()) {
        if (constraints.vertex_map[v]) fixed.emplace_back(v, mesh.point(v));
    }
    if (features) {
        require(constraints.feature_edge_count == 12, "cube sharp edges not detected");
        require(fixed.size() == 8, "cube corners not fixed");
    }
    auto np = CGAL::parameters::number_of_iterations(iterations)
        .number_of_relaxation_steps(relax)
        .edge_is_constrained_map(constraints.edge_map)
        .vertex_is_constrained_map(constraints.vertex_map)
        .collapse_constraints(false)
        .relax_constraints(true)
        .do_project(project);

    if (uniform) {
        double length = 0.45; // Exercise the scalar convenience overload too.
        rar::dispatch_remeshing(faces(mesh), length, mesh, backend, np);
    } else {
        auto lengths = mesh.add_property_map<Mesh::Vertex_index, double>("v:test_length", 0.5).first;
        for (const auto v : mesh.vertices()) {
            lengths[v] = 0.45 + 0.12 * std::abs(CGAL::to_double(mesh.point(v).x()));
        }
        rar::PropertyMapSizingField sizing(lengths, 0.5);
        rar::dispatch_remeshing(faces(mesh), sizing, mesh, backend, np);
    }
    for (const auto& entry : fixed) {
        require(!mesh.is_removed(entry.first), "fixed corner removed");
        require(mesh.point(entry.first) == entry.second, "fixed corner moved");
    }
    mesh.collect_garbage();
    require(CGAL::is_valid_polygon_mesh(mesh), "invalid result mesh");
    require(CGAL::is_triangle_mesh(mesh), "non-triangular result mesh");
    for (const auto v : mesh.vertices()) {
        const Point p = mesh.point(v);
        require(std::isfinite(CGAL::to_double(p.x())) &&
                std::isfinite(CGAL::to_double(p.y())) &&
                std::isfinite(CGAL::to_double(p.z())), "non-finite position");
    }
}

void compare(const Mesh& a, const Mesh& b)
{
    require(a.number_of_vertices() == b.number_of_vertices(), "vertex count differs");
    require(a.number_of_edges() == b.number_of_edges(), "edge count differs");
    require(a.number_of_faces() == b.number_of_faces(), "face count differs");
    auto ac = a.property_map<Mesh::Edge_index, bool>("e:sharp_feature").value();
    auto bc = b.property_map<Mesh::Edge_index, bool>("e:sharp_feature").value();
    for (const auto v : a.vertices()) {
        // Fixtures have scale O(1); tolerance is 1e-12 in coordinate distance.
        require(CGAL::to_double(CGAL::squared_distance(a.point(v), b.point(v))) <= 1e-24,
                "vertex position differs");
    }
    for (const auto h : a.halfedges()) {
        require(a.target(h) == b.target(h) && a.next(h) == b.next(h) &&
                a.face(h) == b.face(h), "halfedge connectivity differs");
    }
    for (const auto e : a.edges()) require(ac[e] == bc[e], "edge constraint differs");
}

void test_naming()
{
    rar::RemeshConfig cfg;
    cfg.input_path = "test.off";
    require(cfg.remesh_backend == rar::RemeshBackend::CGAL, "default backend changed");
    const auto original = rar::make_auto_output_path(cfg);
    require(original.find("__backend-cgal") != std::string::npos, "missing CGAL suffix");
    cfg.remesh_backend = rar::RemeshBackend::Custom;
    const auto custom = rar::make_auto_output_path(cfg);
    require(custom != original && custom.find("__backend-custom") != std::string::npos,
            "backends share automatic output name");
}
} // namespace

int main()
{
    try {
        test_naming();
        const std::array<Mesh, 3> inputs{closed_surface(), open_grid(), sharp_cube()};
        unsigned int cases = 0;
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            for (bool project : {false, true}) {
                for (bool uniform : {false, true}) {
                    for (unsigned int iterations : {1u, 5u}) {
                        Mesh original = inputs[i], custom = inputs[i];
                        const unsigned int relax = iterations == 1 ? 0 : 3;
                        run(original, rar::RemeshBackend::CGAL, i == 2, project, uniform, iterations, relax);
                        run(custom, rar::RemeshBackend::Custom, i == 2, project, uniform, iterations, relax);
                        compare(original, custom);
                        require(CGAL::is_closed(custom) == CGAL::is_closed(inputs[i]), "boundary status changed");
                        require(custom.number_of_vertices() != inputs[i].number_of_vertices(), "fixture did not exercise remeshing");
                        ++cases;
                        std::cout << "fixture=" << i << " project=" << project
                                  << " uniform=" << uniform << " iterations=" << iterations
                                  << " vertices=" << custom.number_of_vertices() << '\n';
                    }
                }
            }
        }
        // Both copies can also accept an empty face range.
        Mesh original, custom;
        double length = 0.5;
        const auto np = CGAL::parameters::number_of_iterations(1);
        rar::dispatch_remeshing(faces(original), length, original, rar::RemeshBackend::CGAL, np);
        rar::dispatch_remeshing(faces(custom), length, custom, rar::RemeshBackend::Custom, np);
        require(original.is_empty() && custom.is_empty(), "empty input changed");
        std::cout << "Passed " << cases << " backend parity cases plus empty input and naming checks.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
