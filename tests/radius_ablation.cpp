#include "rar/RadiusAblation.h"
#include <CGAL/boost/graph/generators.h>
#include <iostream>
#include <limits>

void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main() {
    try {
        for (auto values : {std::pair<double,double>{0.1,0.3}, {0.3,0.1}, {0.2,0.2}}) {
            auto h = rar::separate_radius_mechanisms(values.first, values.second);
            check(h.cap <= h.raw && h.recovery >= h.raw, "direction invariants");
            check(h.cap == (h.radius < h.raw ? h.radius : h.raw), "cap only");
            check(h.recovery == (h.radius > h.raw ? h.radius : h.raw), "recovery only");
        }
        for (double k : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
            check(rar::cgal_adaptive_radius_target_length(k, .01, .001, .5) == .5,
                  "radius invalid/flat curvature fallback");
        }
        bool rejected = false;
        try { rar::separate_radius_mechanisms(std::numeric_limits<double>::quiet_NaN(), .1); }
        catch (const std::runtime_error&) { rejected = true; }
        check(rejected, "invalid CGAL raw must be reported");
        check(rar::cgal_adaptive_radius_target_length(1e9, .01, .001, .5) == .001, "min clip");
        check(rar::cgal_adaptive_radius_target_length(.001, .01, .001, .5) == .5, "max clip");
        const double e = .183 / 2;
        check(rar::cgal_adaptive_radius_target_length(2, e, .001, 2) ==
              std::sqrt(6*e*.5 - 3*e*e), "threshold retains original branch");
        check(rar::cgal_adaptive_radius_target_length(2, std::nextafter(e, 1.0), .001, 2) == .5,
              "above threshold uses radius");

        rar::Mesh mesh;
        CGAL::make_tetrahedron(rar::Point(1,1,1), rar::Point(-1,-1,1),
                              rar::Point(-1,1,-1), rar::Point(1,-1,-1), mesh);
        rar::RemeshConfig cfg; cfg.min_edge_length = .001; cfg.max_edge_length = 10;
        auto curv = rar::make_cgal_adaptive_radius_field(mesh, cfg);
        double min_k = 1e100;
        for (auto v : vertices(mesh)) min_k = (std::min)(min_k,
            rar::cgal_adaptive_radius_curvature(curv.curvature_map[v]));
        check(min_k > 0, "tetrahedron curvature");
        cfg.epsilon = 3/min_k;
        auto fields = rar::make_radius_ablation_field(mesh, cfg);
        CGAL::Polygon_mesh_processing::Adaptive_sizing_field<rar::Mesh> direct(
            cfg.epsilon, {cfg.min_edge_length,cfg.max_edge_length}, faces(mesh), mesh);
        for (auto v : vertices(mesh)) {
            const double k = rar::cgal_adaptive_radius_curvature(fields.radius.curvature_map[v]);
            auto h = rar::separate_radius_mechanisms(fields.raw[v], fields.radius.target_map[v]);
            check(cfg.epsilon*k > 2, "t > 2 fixture");
            check(h.raw == CGAL::to_double(direct.at(v,mesh)), "raw copied from CGAL at()");
            check(h.raw == cfg.min_edge_length, "CGAL negative radicand clamps to min");
            check(h.radius > h.raw && h.cap == h.raw && h.recovery == h.radius, "t > 2 recovery");
        }
        // A flat triangle verifies the actual CGAL zero-curvature behavior.
        rar::Mesh flat;
        auto a=flat.add_vertex({0,0,0}), b=flat.add_vertex({1,0,0}), c=flat.add_vertex({0,1,0});
        flat.add_face(a,b,c);
        auto ff = rar::make_radius_ablation_field(flat,cfg);
        for (auto v : vertices(flat)) check(ff.raw[v] == cfg.max_edge_length &&
            ff.radius.target_map[v] == cfg.max_edge_length, "flat raw/radius max clip");
        cfg.min_edge_length = cfg.max_edge_length = .1;
        auto clipped = rar::make_radius_ablation_field(mesh,cfg);
        for (auto v : vertices(mesh)) {
            auto h=rar::separate_radius_mechanisms(clipped.raw[v],clipped.radius.target_map[v]);
            check(h.raw==.1 && h.radius==.1 && h.cap==.1 && h.recovery==.1, "clipping cancels mechanism");
        }
        // Actual remeshing adapter control: identical inputs and scalar targets.
        cfg = rar::RemeshConfig{}; cfg.epsilon=.05; cfg.iterations=2; cfg.export_field=false;
        rar::Mesh native=mesh, pmap=mesh;
        rar::run_cgal_adaptive_remeshing(native,cfg);
        cfg.field_type=rar::FieldType::CGALAdaptivePmap;
        rar::run_radius_ablation_remeshing(pmap,cfg);
        check(num_vertices(native)==num_vertices(pmap) && num_faces(native)==num_faces(pmap),
              "adapter topology count regression");
        native.collect_garbage(); pmap.collect_garbage();
        for (auto v : vertices(native)) check(native.point(v)==pmap.point(v), "adapter geometry regression");
        std::cout << "Radius ablation tests passed (including actual CGAL t>2 and adapter remeshing)\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
