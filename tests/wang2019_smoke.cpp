#include "rar/Wang2019Remesher.h"
#include <CGAL/boost/graph/helpers.h>
#include <cmath>
#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>

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
    // Exercise the actual angle-driven insertion policy, not just raw Euler.
    {
        Mesh probe(mesh);
        Tree reference(faces(probe).first, faces(probe).second, probe);
        reference.accelerate_distance_queries();
        Options insertion;
        insertion.max_angle = 86.0;
        insertion.protect_features = false;
        insertion.project = false;
        RejectionStats reasons;
        std::size_t angle_flips = 0;
        const auto before_v = probe.number_of_vertices();
        std::vector<Face> candidates;
        for (const Face face : faces(probe)) candidates.push_back(face);
        bool success = false;
        for (const Face face : candidates) {
            if (insert_at_large_angle(probe, face, reference,
                                      insertion, {}, &reasons, &angle_flips)) {
                success = true;
                break;
            }
        }
        // Low-valence tetrahedra need not have a feasible split+flip plan.
        // A rejected proposal must be a true no-op; an accepted proposal
        // must preserve connectivity and add exactly one live vertex.
        if (probe.number_of_vertices() != before_v + (success ? 1u : 0u) ||
            !probe.is_valid() || !CGAL::is_triangle_mesh(probe))
            throw std::runtime_error("Angle-driven proposal violated a mesh invariant");
        std::cout << "Tetra split accepted=" << success << " reasons:";
        for (std::size_t i=0; i<reasons.rejected.size(); ++i)
            std::cout << " " << cause_name(i) << "=" << reasons.rejected[i];
        std::cout << '\n';
    }
    // An interior obtuse patch must exercise the actual angle-driven insertion
    // policy. Unlike the tetrahedron, its interior vertices have legal
    // alternative local connectivities.
    {
        Mesh grid;
        std::array<std::array<Vertex, 5>, 5> vv;
        for (int y=0;y<5;++y)
            for (int x=0;x<5;++x) {
                const double px=(x==2 && y==2) ? 2.65 : static_cast<double>(x);
                const double py=(x==2 && y==2) ? 2.10 : static_cast<double>(y);
                vv[y][x]=grid.add_vertex(Point(px,py,0));
            }
        for (int y=0;y<4;++y)
            for (int x=0;x<4;++x) {
                const Vertex a=vv[y][x], b=vv[y][x+1],
                             c=vv[y+1][x+1], d=vv[y+1][x];
                if (grid.add_face(a,b,c)==Mesh::null_face() ||
                    grid.add_face(a,c,d)==Mesh::null_face())
                    throw std::runtime_error("Failed to generate planar grid");
            }
        if (!grid.is_valid() || !CGAL::is_triangle_mesh(grid))
            throw std::runtime_error("Invalid planar insertion test grid");
        Tree reference(faces(grid).first,faces(grid).second,grid);
        reference.accelerate_distance_queries();
        Options trial;
        trial.min_angle=30.0;
        trial.max_angle=90.0;
        trial.protect_features=false;
        trial.project=false;
        RejectionStats rejected;
        bool inserted=false;
        std::size_t angle_flips=0;
        // Verify Sec. 4.2.4 threshold behavior independently of collapsing:
        // beta_min is preserved while enough triangles already violate it,
        // and is raised only when more deletion candidates are needed.
        std::vector<Face> sorted_faces;
        for (const Face f : faces(grid)) sorted_faces.push_back(f);
        std::sort(sorted_faces.begin(), sorted_faces.end(),
                  [&](Face x, Face y) {
                      return angles_of(grid,x).min<angles_of(grid,y).min;
                  });
        const double min_angle=angles_of(grid,sorted_faces[0]).min;
        const double unchanged=small_angle_threshold_for_count(grid,sorted_faces,60.0,1);
        const double relaxed=small_angle_threshold_for_count(grid,sorted_faces,1.0,1);
        const double none=small_angle_threshold_for_count(grid,sorted_faces,30.0,0);
        if (unchanged!=60.0 || !(relaxed>min_angle) || none!=30.0)
            throw std::runtime_error("Incorrect beta_min candidate relaxation");
        const auto n0=grid.number_of_vertices();
        const auto f0=grid.number_of_faces();
        std::vector<Face> bad;
        for (const Face f : faces(grid))
            if (angles_of(grid,f).max>trial.max_angle) bad.push_back(f);
        for (const Face f : bad) {
            if (insert_at_large_angle(grid,f,reference,trial,{},
                                      &rejected,&angle_flips)) {
                inserted=true;
                break;
            }
        }
        if (!inserted) {
            std::ostringstream os;
            os << "All planar-grid insertion candidates rejected:";
            for (std::size_t i=0;i<rejected.rejected.size();++i)
                os << ' ' << cause_name(i) << '=' << rejected.rejected[i];
            throw std::runtime_error(os.str());
        }
        if (!grid.is_valid() || !CGAL::is_triangle_mesh(grid) ||
            grid.number_of_vertices()!=n0+1 ||
            grid.number_of_faces()!=f0+2)
            throw std::runtime_error("Angle-driven planar insertion corrupted topology");
        // Exercise the complete staged pipeline on a small open surface.
        // In non-strict mode, verify that reported insert/collapse counts
        // agree with the actual live element delta.
        Options full;
        full.min_angle=30.0;
        full.max_angle=90.0;
        full.rounds=3;
        full.k=2;
        full.smoothing_steps=1;
        full.project=false;
        full.protect_features=false;
        full.strict_vertex_count=false;
        full.verbose=false;
        const Statistics actual=run(grid,full);
        const auto delta_v=static_cast<std::ptrdiff_t>(actual.vertices_after)-
                           static_cast<std::ptrdiff_t>(actual.vertices_before);
        if (delta_v!=static_cast<std::ptrdiff_t>(actual.insertions)-
                     static_cast<std::ptrdiff_t>(actual.collapses) ||
            !grid.is_valid() || !CGAL::is_triangle_mesh(grid))
            throw std::runtime_error("Staged pipeline violated live vertex accounting");
        std::cout << "Grid pipeline PASS, inserts=" << actual.insertions
                  << ", collapses=" << actual.collapses
                  << ", flips=" << actual.flips << '\n';

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
