#pragma once

#include "rar/Types.h"

#include <optional>

namespace rar {

struct CSFFieldStats {
    double base_edge_length = 0.0;
    double raw_curvature_min = 0.0;
    double raw_curvature_max = 0.0;
    double raw_curvature_mean = 0.0;
    double smoothed_curvature_min = 0.0;
    double smoothed_curvature_max = 0.0;
    double smoothed_curvature_mean = 0.0;
    double sizing_min = 0.0;
    double sizing_max = 0.0;
    double sizing_mean = 0.0;
    std::size_t vertex_count = 0;
};

class CSFSizingField {
public:
    using vertex_descriptor =
        boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor =
        boost::graph_traits<Mesh>::halfedge_descriptor;
    using FT = Kernel::FT;
    using Point_3 = Point;
    using ScalarMap =
        Mesh::Property_map<vertex_descriptor, double>;

    explicit CSFSizingField(
        double mesh_scale,
        Mesh& mesh);

    FT at(
        vertex_descriptor v,
        const Mesh& mesh) const;

    std::optional<FT> is_too_long(
        vertex_descriptor va,
        vertex_descriptor vb,
        const Mesh& mesh) const;

    std::optional<FT> is_too_short(
        halfedge_descriptor h,
        const Mesh& mesh) const;

    Point_3 split_placement(
        halfedge_descriptor h,
        const Mesh& mesh) const;

    void register_split_vertex(
        vertex_descriptor v,
        const Mesh& mesh);

    const CSFFieldStats& stats() const;

    double raw_curvature(
        vertex_descriptor v) const;

    double smoothed_curvature(
        vertex_descriptor v) const;

    double target_length(
        vertex_descriptor v) const;

private:
    void compute_initial_field(Mesh& mesh);

    double mesh_scale_;
    ScalarMap raw_curvature_map_;
    ScalarMap smoothed_curvature_map_;
    ScalarMap sizing_map_;
    CSFFieldStats stats_;
};

} // namespace rar
