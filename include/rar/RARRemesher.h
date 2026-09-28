#pragma once

#include "rar/FieldExport.h"
#include "rar/RARSizingField.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <optional>
#include <utility>

namespace rar {

namespace PMP = CGAL::Polygon_mesh_processing;

inline RARFieldStats run_rar_like_field_cgal_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg,
    const std::optional<double> chen_beta)
{
    RARSizingField sizing_field(
        cfg.epsilon,
        {cfg.min_edge_length, cfg.max_edge_length},
        mesh,
        chen_beta);

    const RARFieldStats initial_stats = sizing_field.stats();

    if (cfg.export_field && !cfg.export_field_prefix.empty()) {
        if (chen_beta.has_value()) {
            write_rar_chen_field_diagnostics(
                mesh,
                cfg.export_field_prefix,
                [&](const Mesh::Vertex_index v) {
                    return sizing_field.curvature(v);
                },
                [&](const Mesh::Vertex_index v) {
                    return sizing_field.raw_target_length(v);
                },
                [&](const Mesh::Vertex_index v) {
                    return sizing_field.target_length(v);
                });
        } else {
            write_field_diagnostics(
                mesh,
                cfg.export_field_prefix,
                [&](const Mesh::Vertex_index v) {
                    return sizing_field.curvature(v);
                },
                [&](const Mesh::Vertex_index v) {
                    return sizing_field.target_length(v);
                });
        }
    }

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::number_of_iterations(cfg.iterations)
            .number_of_relaxation_steps(cfg.relaxation_steps)
            .do_project(cfg.do_project)
    );

    return initial_stats;
}

inline RARFieldStats run_rar_field_cgal_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    return run_rar_like_field_cgal_remeshing(
        mesh,
        cfg,
        std::nullopt);
}

inline RARFieldStats run_rar_chen_field_cgal_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    return run_rar_like_field_cgal_remeshing(
        mesh,
        cfg,
        cfg.chen_beta);
}

} // namespace rar
