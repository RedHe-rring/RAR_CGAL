#pragma once

#include "rar/CGALAdaptiveRemesher.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>

namespace rar {

struct RadiusMechanisms {
    double raw, radius, cap, recovery;
};

inline RadiusMechanisms separate_radius_mechanisms(double raw, double radius)
{
    // CGAL can return NaN for invalid curvature. Do not silently reinterpret it.
    if (!std::isfinite(raw) || !(raw > 0.0) ||
        !std::isfinite(radius) || !(radius > 0.0))
        throw std::runtime_error("Non-finite or non-positive initial sizing value");
    return {raw, radius, (std::min)(raw, radius), (std::max)(raw, radius)};
}

inline const char* radius_mechanism(const RadiusMechanisms& h)
{
    // Diagnostics only: ignore arithmetic roundoff from 6*eps/k vs 6*eps*(1/k).
    const double tol = 64 * std::numeric_limits<double>::epsilon() *
                       (std::max)(h.raw, h.radius);
    if (h.radius < h.raw - tol) return "Cap";
    if (h.radius > h.raw + tol) return "Recovery";
    return "Unchanged";
}

struct RadiusAblationField {
    CGALAdaptiveRadiusField radius;
    PropertyMapSizingField::ScalarMap raw;
};

inline RadiusAblationField make_radius_ablation_field(Mesh& mesh, const RemeshConfig& cfg)
{
    validate_cgal_adaptive_config(cfg);
    PMP::Adaptive_sizing_field<Mesh> raw_field(
        cfg.epsilon, {cfg.min_edge_length, cfg.max_edge_length}, faces(mesh), mesh);
    auto raw = mesh.add_property_map<Mesh::Vertex_index, double>(
        "v:radius_ablation_raw", cfg.max_edge_length).first;
    for (auto v : vertices(mesh)) raw[v] = CGAL::to_double(raw_field.at(v, mesh));
    auto radius = make_cgal_adaptive_radius_field(mesh, cfg);
    for (auto v : vertices(mesh)) separate_radius_mechanisms(raw[v], radius.target_map[v]);
    return {radius, raw};
}

inline void write_radius_ablation_diagnostics(Mesh& mesh, const RemeshConfig& cfg)
{
    auto field = make_radius_ablation_field(mesh, cfg);
    auto area = mesh.add_property_map<Mesh::Vertex_index, double>("v:ablation_area", 0).first;
    for (auto v : vertices(mesh)) area[v] = 0;
    for (auto f : faces(mesh)) {
        auto h = halfedge(f, mesh);
        auto a = source(h, mesh), b = target(h, mesh), c = target(next(h, mesh), mesh);
        const double share = std::sqrt(CGAL::to_double(CGAL::cross_product(
            mesh.point(b) - mesh.point(a), mesh.point(c) - mesh.point(a)).squared_length())) / 6;
        area[a] += share; area[b] += share; area[c] += share;
    }
    const std::filesystem::path path(cfg.radius_diagnostics_path);
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open radius diagnostics: " + path.string());
    out << std::setprecision(17)
        << "vertex_id,x,y,z,vertex_area,kappa,epsilon_kappa,h_raw,h_radius,h_cap,h_recovery,delta,ratio,mechanism,guard_branch,clipped_unchanged,raw_at_min,raw_at_max,radius_at_min,radius_at_max\n";
    for (auto v : vertices(mesh)) {
        const double k = cgal_adaptive_radius_curvature(field.radius.curvature_map[v]);
        const auto h = separate_radius_mechanisms(field.raw[v], field.radius.target_map[v]);
        const bool branch = std::isfinite(k) && k > 1e-15 &&
            cfg.epsilon > kCGALAdaptiveRadiusThresholdRatio / k;
        const bool clipped = branch && h.raw == h.radius &&
            (h.raw == cfg.min_edge_length || h.raw == cfg.max_edge_length);
        const auto& p = mesh.point(v);
        out << v.idx() << ',' << p.x() << ',' << p.y() << ',' << p.z() << ','
            << area[v] << ',' << k << ',' << cfg.epsilon*k << ','
            << h.raw << ',' << h.radius << ',' << h.cap << ',' << h.recovery << ','
            << h.radius-h.raw << ',' << h.radius/h.raw << ',' << radius_mechanism(h) << ','
            << branch << ',' << clipped << ',' << (h.raw == cfg.min_edge_length) << ','
            << (h.raw == cfg.max_edge_length) << ',' << (h.radius == cfg.min_edge_length) << ','
            << (h.radius == cfg.max_edge_length) << '\n';
    }
    if (!out) throw std::runtime_error("Failed writing radius diagnostics: " + path.string());
}

inline void run_radius_ablation_remeshing(Mesh& mesh, const RemeshConfig& cfg)
{
    const auto field = make_radius_ablation_field(mesh, cfg);
    auto target = mesh.add_property_map<Mesh::Vertex_index, double>(
        "v:radius_ablation_target", cfg.max_edge_length).first;
    for (auto v : vertices(mesh)) {
        const auto h = separate_radius_mechanisms(field.raw[v], field.radius.target_map[v]);
        switch (cfg.field_type) {
        case FieldType::CGALAdaptivePmap: target[v] = h.raw; break;
        case FieldType::CGALAdaptiveRadiusCap: target[v] = h.cap; break;
        case FieldType::CGALAdaptiveRadiusRecovery: target[v] = h.recovery; break;
        default: throw std::invalid_argument("Unsupported radius ablation mode");
        }
    }
    if (cfg.export_field && !cfg.export_field_prefix.empty())
        write_field_diagnostics(mesh, cfg.export_field_prefix,
            [&](auto v) { return cgal_adaptive_radius_curvature(field.radius.curvature_map[v]); },
            [&](auto v) { return target[v]; });
    PropertyMapSizingField sizing_field(target, cfg.max_edge_length);
    const auto constraints = make_feature_constraints(mesh, cfg.preserve_features, cfg.feature_angle_degrees);
    PMP::isotropic_remeshing(faces(mesh), sizing_field, mesh,
        CGAL::parameters::number_of_iterations(cfg.iterations)
        .number_of_relaxation_steps(cfg.relaxation_steps)
        .edge_is_constrained_map(constraints.edge_map)
        .vertex_is_constrained_map(constraints.vertex_map)
        .collapse_constraints(false).relax_constraints(true).do_project(cfg.do_project));
}

} // namespace rar
