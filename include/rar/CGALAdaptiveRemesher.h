#pragma once

#include "rar/CGALAdaptiveSizing.h"
#include "rar/ChenSizingCorrection.h"
#include "rar/FeatureConstraints.h"
#include "rar/FieldExport.h"
#include "rar/PropertyMapSizingField.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rar {

namespace PMP = CGAL::Polygon_mesh_processing;

struct CGALAdaptiveChenStats {
    std::size_t vertex_count = 0;

    double raw_sizing_min = 0.0;
    double raw_sizing_max = 0.0;
    double raw_sizing_mean = 0.0;

    double corrected_sizing_min = 0.0;
    double corrected_sizing_max = 0.0;
    double corrected_sizing_mean = 0.0;

    double chen_beta = 0.0;
    double chen_objective = 0.0;
    double chen_max_gradient_before = 0.0;
    double chen_max_gradient_after = 0.0;
    std::size_t chen_changed_vertex_count = 0;
};

inline void validate_cgal_adaptive_config(
    const RemeshConfig& cfg)
{
    if (!(cfg.epsilon > 0.0)) {
        throw std::invalid_argument(
            "epsilon must be > 0");
    }
    if (!(cfg.min_edge_length > 0.0) ||
        !(cfg.max_edge_length > 0.0)) {
        throw std::invalid_argument(
            "edge length bounds must be > 0");
    }
    if (cfg.min_edge_length >
        cfg.max_edge_length) {
        throw std::invalid_argument(
            "min_edge_length must be <= max_edge_length");
    }
}

inline void validate_cgal_adaptive_chen_config(
    const RemeshConfig& cfg)
{
    validate_cgal_adaptive_config(cfg);

    if (!(cfg.chen_beta > 1.0) ||
        !std::isfinite(cfg.chen_beta)) {
        throw std::invalid_argument(
            "beta must be > 1 for CGAL adaptive Chen modes");
    }
}

template <typename SizingField, typename TargetFn>
inline void export_cgal_adaptive_field(
    Mesh& mesh,
    const RemeshConfig& cfg,
    const SizingField&,
    TargetFn target_writer)
{
    if (!cfg.export_field ||
        cfg.export_field_prefix.empty()) {
        return;
    }

    using PrincipalCurvatures =
        PMP::Principal_curvatures_and_directions<Kernel>;

    const auto curvature_property =
        mesh.add_property_map<
            Mesh::Vertex_index,
            PrincipalCurvatures>(
                "v:cgal_principal_curvatures",
                PrincipalCurvatures{});
    const auto curvature_map =
        curvature_property.first;

    if (curvature_property.second) {
        PMP::interpolated_corrected_curvatures(
            mesh,
            CGAL::parameters::
                vertex_principal_curvatures_and_directions_map(
                    curvature_map));
    }

    target_writer(
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
        });
}

inline void run_cgal_adaptive_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    validate_cgal_adaptive_config(cfg);

    const std::pair<double, double>
        edge_min_max{
            cfg.min_edge_length,
            cfg.max_edge_length
        };

    PMP::Adaptive_sizing_field<Mesh>
        sizing_field(
            cfg.epsilon,
            edge_min_max,
            faces(mesh),
            mesh);

    export_cgal_adaptive_field(
        mesh,
        cfg,
        sizing_field,
        [&](auto curvature_fn) {
            write_field_diagnostics(
                mesh,
                cfg.export_field_prefix,
                curvature_fn,
                [&](const Mesh::Vertex_index v) {
                    return CGAL::to_double(
                        sizing_field.at(v, mesh));
                });
        });

    const FeatureConstraints
        feature_constraints =
            make_feature_constraints(
                mesh,
                cfg.preserve_features,
                cfg.feature_angle_degrees);

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::
            number_of_iterations(
                cfg.iterations)
            .number_of_relaxation_steps(
                cfg.relaxation_steps)
            .edge_is_constrained_map(
                feature_constraints.edge_map)
            .vertex_is_constrained_map(
                feature_constraints.vertex_map)
            .collapse_constraints(false)
            .relax_constraints(true)
            .do_project(cfg.do_project)
    );
}

using CGALPrincipalCurvatures =
    PMP::Principal_curvatures_and_directions<Kernel>;

struct CGALAdaptiveRadiusField {
    Mesh::Property_map<Mesh::Vertex_index, CGALPrincipalCurvatures>
        curvature_map;
    Mesh::Property_map<Mesh::Vertex_index, double> target_map;
};

inline double cgal_adaptive_radius_curvature(
    const CGALPrincipalCurvatures& pc)
{
    return (std::max)(
        std::abs(CGAL::to_double(pc.min_curvature)),
        std::abs(CGAL::to_double(pc.max_curvature)));
}

inline CGALAdaptiveRadiusField make_cgal_adaptive_radius_field(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    const auto curvature_property =
        mesh.add_property_map<
            Mesh::Vertex_index,
            CGALPrincipalCurvatures>(
                "v:cgal_principal_curvatures",
                CGALPrincipalCurvatures{});
    const auto curvature_map = curvature_property.first;

    if (curvature_property.second) {
        PMP::interpolated_corrected_curvatures(
            mesh,
            CGAL::parameters::
                vertex_principal_curvatures_and_directions_map(
                    curvature_map));
    }

    auto target_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            double>(
                "v:cgal_adaptive_radius_raw_target_length",
                cfg.max_edge_length).first;

    for (const Mesh::Vertex_index v : vertices(mesh)) {
        target_map[v] = cgal_adaptive_radius_target_length(
            cgal_adaptive_radius_curvature(curvature_map[v]),
            cfg.epsilon,
            cfg.min_edge_length,
            cfg.max_edge_length);
    }

    return {curvature_map, target_map};
}

inline void run_cgal_adaptive_radius_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    validate_cgal_adaptive_config(cfg);

    const CGALAdaptiveRadiusField field =
        make_cgal_adaptive_radius_field(mesh, cfg);

    if (cfg.export_field && !cfg.export_field_prefix.empty()) {
        write_field_diagnostics(
            mesh,
            cfg.export_field_prefix,
            [&](const Mesh::Vertex_index v) {
                return cgal_adaptive_radius_curvature(
                    field.curvature_map[v]);
            },
            [&](const Mesh::Vertex_index v) {
                return field.target_map[v];
            });
    }

    PropertyMapSizingField sizing_field(
        field.target_map,
        cfg.max_edge_length);

    const FeatureConstraints feature_constraints =
        make_feature_constraints(
            mesh,
            cfg.preserve_features,
            cfg.feature_angle_degrees);

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::
            number_of_iterations(cfg.iterations)
            .number_of_relaxation_steps(cfg.relaxation_steps)
            .edge_is_constrained_map(feature_constraints.edge_map)
            .vertex_is_constrained_map(feature_constraints.vertex_map)
            .collapse_constraints(false)
            .relax_constraints(true)
            .do_project(cfg.do_project)
    );
}

template <typename ScalarMap, typename ExportFn>
inline CGALAdaptiveChenStats
run_cgal_adaptive_chen_from_maps(
    Mesh& mesh,
    const RemeshConfig& cfg,
    const ScalarMap raw_map,
    const ScalarMap corrected_map,
    ExportFn export_field)
{
    CGALAdaptiveChenStats stats;
    stats.vertex_count =
        num_vertices(mesh);
    stats.chen_beta =
        cfg.chen_beta;

    double raw_sum = 0.0;
    stats.raw_sizing_min =
        std::numeric_limits<double>::infinity();
    stats.raw_sizing_max = 0.0;

    for (const Mesh::Vertex_index v :
         vertices(mesh)) {
        const double value = raw_map[v];

        corrected_map[v] = value;

        raw_sum += value;
        stats.raw_sizing_min =
            (std::min)(
                stats.raw_sizing_min,
                value);
        stats.raw_sizing_max =
            (std::max)(
                stats.raw_sizing_max,
                value);
    }

    if (stats.vertex_count > 0) {
        stats.raw_sizing_mean =
            raw_sum /
            static_cast<double>(
                stats.vertex_count);
    } else {
        stats.raw_sizing_min = 0.0;
    }

    const ChenCorrectionStats chen_stats =
        apply_chen_sizing_correction(
            mesh,
            corrected_map,
            cfg.chen_beta,
            cfg.min_edge_length);

    stats.chen_objective =
        chen_stats.objective;
    stats.chen_max_gradient_before =
        chen_stats.max_gradient_before;
    stats.chen_max_gradient_after =
        chen_stats.max_gradient_after;
    stats.chen_changed_vertex_count =
        chen_stats.changed_vertex_count;

    double corrected_sum = 0.0;
    stats.corrected_sizing_min =
        std::numeric_limits<double>::infinity();
    stats.corrected_sizing_max = 0.0;

    for (const Mesh::Vertex_index v :
         vertices(mesh)) {
        const double value =
            corrected_map[v];

        corrected_sum += value;
        stats.corrected_sizing_min =
            (std::min)(
                stats.corrected_sizing_min,
                value);
        stats.corrected_sizing_max =
            (std::max)(
                stats.corrected_sizing_max,
                value);
    }

    if (stats.vertex_count > 0) {
        stats.corrected_sizing_mean =
            corrected_sum /
            static_cast<double>(
                stats.vertex_count);
    } else {
        stats.corrected_sizing_min = 0.0;
    }

    export_field(raw_map, corrected_map);

    PropertyMapSizingField sizing_field(
        corrected_map,
        cfg.max_edge_length);

    const FeatureConstraints
        feature_constraints =
            make_feature_constraints(
                mesh,
                cfg.preserve_features,
                cfg.feature_angle_degrees);

    PMP::isotropic_remeshing(
        faces(mesh),
        sizing_field,
        mesh,
        CGAL::parameters::
            number_of_iterations(
                cfg.iterations)
            .number_of_relaxation_steps(
                cfg.relaxation_steps)
            .edge_is_constrained_map(
                feature_constraints.edge_map)
            .vertex_is_constrained_map(
                feature_constraints.vertex_map)
            .collapse_constraints(false)
            .relax_constraints(true)
            .do_project(cfg.do_project)
    );

    return stats;
}

inline CGALAdaptiveChenStats
run_cgal_adaptive_chen_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    validate_cgal_adaptive_chen_config(cfg);

    const std::pair<double, double>
        edge_min_max{
            cfg.min_edge_length,
            cfg.max_edge_length
        };

    PMP::Adaptive_sizing_field<Mesh>
        raw_field(
            cfg.epsilon,
            edge_min_max,
            faces(mesh),
            mesh);

    auto raw_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            double>(
                "v:cgal_adaptive_raw_target_length",
                cfg.max_edge_length).first;

    auto corrected_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            double>(
                "v:cgal_adaptive_chen_target_length",
                cfg.max_edge_length).first;

    for (const Mesh::Vertex_index v : vertices(mesh)) {
        raw_map[v] =
            CGAL::to_double(raw_field.at(v, mesh));
    }

    return run_cgal_adaptive_chen_from_maps(
        mesh,
        cfg,
        raw_map,
        corrected_map,
        [&](const auto raw_values,
            const auto corrected_values) {
            export_cgal_adaptive_field(
                mesh,
                cfg,
                raw_field,
                [&](auto curvature_fn) {
                    write_rar_chen_field_diagnostics(
                        mesh,
                        cfg.export_field_prefix,
                        curvature_fn,
                        [&](const Mesh::Vertex_index v) {
                            return raw_values[v];
                        },
                        [&](const Mesh::Vertex_index v) {
                            return corrected_values[v];
                        });
                });
        });
}

inline CGALAdaptiveChenStats
run_cgal_adaptive_radius_chen_remeshing(
    Mesh& mesh,
    const RemeshConfig& cfg)
{
    validate_cgal_adaptive_chen_config(cfg);

    const CGALAdaptiveRadiusField field =
        make_cgal_adaptive_radius_field(mesh, cfg);

    auto corrected_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            double>(
                "v:cgal_adaptive_radius_chen_target_length",
                cfg.max_edge_length).first;

    return run_cgal_adaptive_chen_from_maps(
        mesh,
        cfg,
        field.target_map,
        corrected_map,
        [&](const auto raw_values,
            const auto corrected_values) {
            if (!cfg.export_field ||
                cfg.export_field_prefix.empty()) {
                return;
            }

            write_rar_chen_field_diagnostics(
                mesh,
                cfg.export_field_prefix,
                [&](const Mesh::Vertex_index v) {
                    return cgal_adaptive_radius_curvature(
                        field.curvature_map[v]);
                },
                [&](const Mesh::Vertex_index v) {
                    return raw_values[v];
                },
                [&](const Mesh::Vertex_index v) {
                    return corrected_values[v];
                });
        });
}

} // namespace rar
