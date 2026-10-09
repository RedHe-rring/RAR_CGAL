#include "rar/SDFThicknessField.h"
#include "rar/FeatureConstraints.h"
#include "rar/FieldExport.h"
#include "rar/PropertyMapSizingField.h"
#include "rar/RemeshBackend.h"

#include <CGAL/version.h>
#include <fstream>
#include <iomanip>
#include <locale>

namespace rar {

void write_sdf_thickness_field(const Mesh& mesh, const RemeshConfig& cfg,
                             const SDFThicknessField& field)
{
    const auto& prefix = cfg.export_field_prefix;
    if (!cfg.export_field || prefix.empty()) return;
    detail::write_scalar_field_diagnostics(mesh, prefix, {
        {"curvature", [&](auto v) { return field.curvature[v]; }},
        {"target_length", [&](auto v) { return field.target[v]; }},
        {"thickness", [&](auto v) { return field.thickness[v]; }},
        {"vertex_area", [&](auto v) { return field.vertex_area[v]; }},
        {"epsilon_local", [&](auto v) { return field.samples[v].epsilon_local; }},
        {"epsilon_kappa", [&](auto v) { return cfg.epsilon * field.curvature[v]; }},
        {"epsilon_local_kappa", [&](auto v) { return field.samples[v].epsilon_local * field.curvature[v]; }},
        {"epsilon_limited", [&](auto v) { return double(field.samples[v].epsilon_limited); }},
        {"radius_limited", [&](auto v) { return double(field.samples[v].radius_limited); }},
        {"thickness_limited", [&](auto v) { return double(field.samples[v].thickness_limited); }},
        {"min_edge_clipped", [&](auto v) { return double(field.samples[v].min_edge_clipped); }},
        {"max_edge_clipped", [&](auto v) { return double(field.samples[v].max_edge_clipped); }},
        {"zero_curvature", [&](auto v) { return double(field.curvature[v] == 0.0); }}
    });
    std::ofstream faces_out;
    faces_out.exceptions(std::ios::failbit | std::ios::badbit);
    faces_out.open(prefix + "_sdf_faces.csv");
    faces_out.imbue(std::locale::classic());
    faces_out << std::setprecision(17) << "face_id,x,y,z,area,thickness,sdf_valid\n";
    for (const auto f : mesh.faces()) {
        const auto h = mesh.halfedge(f);
        const auto p = CGAL::centroid(mesh.point(mesh.source(h)), mesh.point(mesh.target(h)),
                                      mesh.point(mesh.target(mesh.next(h))));
        faces_out << f.idx() << ',' << p.x() << ',' << p.y() << ',' << p.z() << ','
                  << field.face_area[f] << ',' << field.face_thickness[f] << ",1\n";
    }
    faces_out.close();
    std::ofstream metadata;
    metadata.exceptions(std::ios::failbit | std::ios::badbit);
    metadata.open(prefix + "_sdf.json");
    metadata.imbue(std::locale::classic());
    metadata << std::setprecision(17)
        << "{\n  \"schema_version\": 1,\n  \"field\": \"cgal-adaptive-sdf-thickness\",\n"
        << "  \"cgal_version\": \"" << CGAL_VERSION_STR << "\",\n"
        << "  \"sdf_postprocess\": false,\n  \"thickness_units\": \"input mesh length units\",\n"
        << "  \"vertex_transfer\": \"incident face area weighted mean\",\n"
        << "  \"invalid_sdf_policy\": \"error\",\n  \"chen_applied\": false,\n"
        << "  \"sdf_rays\": " << cfg.sdf_rays << ",\n"
        << "  \"sdf_cone_angle_degrees\": " << cfg.sdf_cone_angle_degrees << ",\n"
        << "  \"radius_epsilon_ratio\": 0.183,\n"
        << "  \"epsilon\": " << cfg.epsilon << ",\n"
        << "  \"min_edge_length\": " << cfg.min_edge_length << ",\n"
        << "  \"max_edge_length\": " << cfg.max_edge_length << ",\n"
        << "  \"thickness_limited_area_fraction\": " << field.stats.thickness_limited_area_fraction << "\n}\n";
    metadata.close();
}

SDFThicknessStats run_sdf_thickness_remeshing(Mesh& mesh, const RemeshConfig& cfg)
{
    const auto field = make_sdf_thickness_field(mesh, cfg);
    // Export initial values before topology/geometry changes. Only the target map
    // is maintained by the existing sizing adapter during remeshing.
    write_sdf_thickness_field(mesh, cfg, field);
    PropertyMapSizingField sizing(field.target, cfg.max_edge_length);
    const auto constraints = make_feature_constraints(mesh, cfg.preserve_features, cfg.feature_angle_degrees);
    dispatch_remeshing(faces(mesh), sizing, mesh, cfg,
        CGAL::parameters::number_of_iterations(cfg.iterations)
        .number_of_relaxation_steps(cfg.relaxation_steps)
        .edge_is_constrained_map(constraints.edge_map)
        .vertex_is_constrained_map(constraints.vertex_map)
        .collapse_constraints(false).relax_constraints(true).do_project(cfg.do_project));
    return field.stats;
}

} // namespace rar
