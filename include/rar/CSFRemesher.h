#pragma once

#include "rar/CSFSizingField.h"
#include "rar/FieldExport.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>

namespace rar {

namespace PMP = CGAL::Polygon_mesh_processing;

inline CSFFieldStats run_csf_field_cgal_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    CSFSizingField sizing_field(
        cfg.csf_mesh_scale,
        mesh);

    const CSFFieldStats initial_stats =
        sizing_field.stats();

    if (cfg.export_field &&
        !cfg.export_field_prefix.empty()) {
        write_csf_field_diagnostics(
            mesh,
            cfg.export_field_prefix,
            [&](const Mesh::Vertex_index v) {
                return sizing_field.raw_curvature(v);
            },
            [&](const Mesh::Vertex_index v) {
                return sizing_field.smoothed_curvature(v);
            },
            [&](const Mesh::Vertex_index v) {
                return sizing_field.target_length(v);
            });
    }

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::
            number_of_iterations(cfg.iterations)
            .number_of_relaxation_steps(
                cfg.relaxation_steps)
            .do_project(cfg.do_project)
    );

    return initial_stats;
}

} // namespace rar
