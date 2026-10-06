#pragma once

#include "rar/CGALAdaptiveSizing.h"
#include "rar/Types.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace rar {

inline constexpr double kAutoEpsilonCurvatureQuantile = 0.90;

struct AutoEpsilonStats {
    double epsilon = 1e-3;
    double weighted_p90_kappa = 0.0;
    double valid_vertex_area = 0.0;
    std::size_t valid_vertex_count = 0;
    bool used_fallback = true;
};

inline AutoEpsilonStats estimate_auto_epsilon(
    Mesh& mesh,
    const double fallback_epsilon = 1e-3)
{
    using PrincipalCurvatures =
        CGAL::Polygon_mesh_processing::
            Principal_curvatures_and_directions<Kernel>;

    auto curvature_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            PrincipalCurvatures>(
                "v:cgal_principal_curvatures",
                PrincipalCurvatures{}).first;
    auto vertex_area_map =
        mesh.add_property_map<
            Mesh::Vertex_index,
            double>(
                "v:auto_epsilon_area",
                0.0).first;

    for (const Mesh::Vertex_index v : vertices(mesh)) {
        vertex_area_map[v] = 0.0;
    }

    CGAL::Polygon_mesh_processing::
        interpolated_corrected_curvatures(
            mesh,
            CGAL::parameters::
                vertex_principal_curvatures_and_directions_map(
                    curvature_map));

    for (const Mesh::Face_index f : faces(mesh)) {
        const auto h = halfedge(f, mesh);
        const Mesh::Vertex_index v0 = source(h, mesh);
        const Mesh::Vertex_index v1 = target(h, mesh);
        const Mesh::Vertex_index v2 =
            target(next(h, mesh), mesh);

        const Kernel::Vector_3 normal =
            CGAL::cross_product(
                Kernel::Vector_3(
                    mesh.point(v0),
                    mesh.point(v1)),
                Kernel::Vector_3(
                    mesh.point(v0),
                    mesh.point(v2)));
        const double normal_sq =
            CGAL::to_double(normal.squared_length());

        if (!(normal_sq > 1e-30) ||
            !std::isfinite(normal_sq)) {
            continue;
        }

        // |cross(e01, e02)| is twice the triangle area. One third of the
        // triangle area is assigned to each incident vertex.
        const double area_share =
            std::sqrt(normal_sq) / 6.0;
        vertex_area_map[v0] += area_share;
        vertex_area_map[v1] += area_share;
        vertex_area_map[v2] += area_share;
    }

    std::vector<std::pair<double, double>> samples;
    samples.reserve(num_vertices(mesh));

    AutoEpsilonStats stats;
    stats.epsilon = fallback_epsilon;

    for (const Mesh::Vertex_index v : vertices(mesh)) {
        const PrincipalCurvatures& pc = curvature_map[v];
        const double kappa =
            (std::max)(
                std::abs(CGAL::to_double(pc.min_curvature)),
                std::abs(CGAL::to_double(pc.max_curvature)));
        const double area = vertex_area_map[v];

        if (!std::isfinite(kappa) ||
            kappa <= 1e-15 ||
            !std::isfinite(area) ||
            !(area > 0.0)) {
            continue;
        }

        samples.emplace_back(kappa, area);
        stats.valid_vertex_area += area;
    }

    stats.valid_vertex_count = samples.size();

    if (samples.empty() ||
        !(stats.valid_vertex_area > 0.0) ||
        !std::isfinite(stats.valid_vertex_area)) {
        return stats;
    }

    std::sort(
        samples.begin(),
        samples.end(),
        [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });

    const double target_area =
        kAutoEpsilonCurvatureQuantile *
        stats.valid_vertex_area;
    double cumulative_area = 0.0;
    double p90_kappa = samples.back().first;

    for (const auto& sample : samples) {
        cumulative_area += sample.second;
        if (cumulative_area >= target_area) {
            p90_kappa = sample.first;
            break;
        }
    }

    const double epsilon =
        kCGALAdaptiveRadiusThresholdRatio /
        p90_kappa;

    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        return stats;
    }

    stats.epsilon = epsilon;
    stats.weighted_p90_kappa = p90_kappa;
    stats.used_fallback = false;
    return stats;
}

} // namespace rar
