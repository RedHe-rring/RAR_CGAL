#pragma once

#include "rar/FieldExport.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rar {

namespace PMP = CGAL::Polygon_mesh_processing;

inline void run_cgal_adaptive_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    if (!(cfg.epsilon > 0.0)) {
        throw std::invalid_argument("epsilon must be > 0");
    }
    if (!(cfg.min_edge_length > 0.0) ||
        !(cfg.max_edge_length > 0.0)) {
        throw std::invalid_argument(
            "edge length bounds must be > 0");
    }
    if (cfg.min_edge_length > cfg.max_edge_length) {
        throw std::invalid_argument(
            "min_edge_length must be <= max_edge_length");
    }

    const std::pair<double, double> edge_min_max{
        cfg.min_edge_length,
        cfg.max_edge_length
    };

    PMP::Adaptive_sizing_field<Mesh> sizing_field(
        cfg.epsilon,
        edge_min_max,
        faces(mesh),
        mesh
    );

    if (cfg.export_field && !cfg.export_field_prefix.empty()) {
        using PrincipalCurvatures =
            PMP::Principal_curvatures_and_directions<Kernel>;

        auto curvature_map =
            mesh.add_property_map<
                Mesh::Vertex_index,
                PrincipalCurvatures>(
                    "v:cgal_principal_curvatures",
                    PrincipalCurvatures{}).first;

        PMP::interpolated_corrected_curvatures(
            mesh,
            CGAL::parameters::
                vertex_principal_curvatures_and_directions_map(
                    curvature_map));

        write_field_diagnostics(
            mesh,
            cfg.export_field_prefix,
            [&](const Mesh::Vertex_index v) {
                const PrincipalCurvatures& pc =
                    curvature_map[v];
                const double kmin =
                    std::abs(CGAL::to_double(
                        pc.min_curvature));
                const double kmax =
                    std::abs(CGAL::to_double(
                        pc.max_curvature));
                return (std::max)(kmin, kmax);
            },
            [&](const Mesh::Vertex_index v) {
                return CGAL::to_double(
                    sizing_field.at(v, mesh));
            });
    }

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::number_of_iterations(
            cfg.iterations)
            .number_of_relaxation_steps(
                cfg.relaxation_steps)
            .do_project(cfg.do_project)
    );
}

} // namespace rar
