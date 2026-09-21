#pragma once

#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <stdexcept>
#include <utility>

namespace rar {

namespace PMP = CGAL::Polygon_mesh_processing;

inline void run_cgal_adaptive_remeshing(Mesh& mesh, const RemeshConfig& cfg) {
    if (!(cfg.epsilon > 0.0)) {
        throw std::invalid_argument("epsilon must be > 0");
    }
    if (!(cfg.min_edge_length > 0.0) || !(cfg.max_edge_length > 0.0)) {
        throw std::invalid_argument("edge length bounds must be > 0");
    }
    if (cfg.min_edge_length > cfg.max_edge_length) {
        throw std::invalid_argument("min_edge_length must be <= max_edge_length");
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

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::number_of_iterations(cfg.iterations)
            .number_of_relaxation_steps(cfg.relaxation_steps)
            .do_project(cfg.do_project)
    );
}

} // namespace rar
