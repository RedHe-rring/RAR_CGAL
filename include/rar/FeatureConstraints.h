#pragma once

#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <cstddef>

namespace rar {

struct FeatureConstraints {
    Mesh::Property_map<Mesh::Edge_index, bool> edge_map;
    Mesh::Property_map<Mesh::Vertex_index, bool> vertex_map;
    std::size_t feature_edge_count = 0;
    std::size_t fixed_feature_vertex_count = 0;
};

inline FeatureConstraints make_feature_constraints(
    Mesh& mesh,
    const bool preserve_features,
    const double feature_angle_degrees)
{
    auto edge_map =
        mesh.add_property_map<Mesh::Edge_index, bool>(
            "e:sharp_feature", false).first;

    auto feature_degree_map =
        mesh.add_property_map<Mesh::Vertex_index, int>(
            "v:sharp_feature_degree", 0).first;

    auto vertex_map =
        mesh.add_property_map<Mesh::Vertex_index, bool>(
            "v:sharp_feature_corner", false).first;

    // Property maps can already exist when a mesh is processed more than
    // once, so explicitly reset them before running feature detection.
    for (const Mesh::Edge_index e : edges(mesh)) {
        edge_map[e] = false;
    }
    for (const Mesh::Vertex_index v : vertices(mesh)) {
        feature_degree_map[v] = 0;
        vertex_map[v] = false;
    }

    if (!preserve_features) {
        return {
            edge_map,
            vertex_map,
            0,
            0
        };
    }

    CGAL::Polygon_mesh_processing::detect_sharp_edges(
        mesh,
        feature_angle_degrees,
        edge_map,
        CGAL::parameters::vertex_feature_degree_map(
            feature_degree_map));

    std::size_t feature_edge_count = 0;
    for (const Mesh::Edge_index e : edges(mesh)) {
        if (edge_map[e]) {
            ++feature_edge_count;
        }
    }

    std::size_t fixed_feature_vertex_count = 0;
    for (const Mesh::Vertex_index v : vertices(mesh)) {
        const int degree = feature_degree_map[v];

        // Degree-2 feature vertices belong to a regular feature polyline and
        // are allowed to relax along that polyline. Endpoints and junctions
        // are fixed so that feature corners cannot drift or disappear.
        if (degree > 0 && degree != 2) {
            vertex_map[v] = true;
            ++fixed_feature_vertex_count;
        }
    }

    return {
        edge_map,
        vertex_map,
        feature_edge_count,
        fixed_feature_vertex_count
    };
}

} // namespace rar
