#pragma once

#include "rar/Types.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <cstddef>
#include <optional>

namespace rar {

// Lightweight PMPSizingField adapter backed by a per-vertex target-length map.
// It is used when an externally constructed field (for example
// CGAL Adaptive_sizing_field followed by Chen correction) must be consumed by
// PMP::isotropic_remeshing().
class PropertyMapSizingField {
public:
    using vertex_descriptor =
        boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor =
        boost::graph_traits<Mesh>::halfedge_descriptor;
    using FT = Kernel::FT;
    using Point_3 = Point;
    using ScalarMap =
        Mesh::Property_map<vertex_descriptor, double>;

    PropertyMapSizingField(
        const ScalarMap sizing_map,
        const double fallback_target_length)
        : sizing_map_(sizing_map),
          fallback_target_length_(fallback_target_length)
    {
    }

    FT at(
        const vertex_descriptor v,
        const Mesh&) const
    {
        return FT(sizing_map_[v]);
    }

    std::optional<FT> is_too_long(
        const vertex_descriptor va,
        const vertex_descriptor vb,
        const Mesh& mesh) const
    {
        const double edge_sq =
            CGAL::to_double(
                CGAL::squared_distance(
                    mesh.point(va),
                    mesh.point(vb)));

        const double target =
            (4.0 / 3.0) *
            (std::min)(
                sizing_map_[va],
                sizing_map_[vb]);
        const double target_sq =
            target * target;

        if (edge_sq > target_sq) {
            return FT(edge_sq / target_sq);
        }
        return std::nullopt;
    }

    std::optional<FT> is_too_short(
        const halfedge_descriptor h,
        const Mesh& mesh) const
    {
        const vertex_descriptor va =
            source(h, mesh);
        const vertex_descriptor vb =
            target(h, mesh);

        const double edge_sq =
            CGAL::to_double(
                CGAL::squared_distance(
                    mesh.point(va),
                    mesh.point(vb)));

        const double target_length =
            (4.0 / 5.0) *
            (std::min)(
                sizing_map_[va],
                sizing_map_[vb]);
        const double target_sq =
            target_length * target_length;

        if (edge_sq < target_sq) {
            return FT(edge_sq / target_sq);
        }
        return std::nullopt;
    }

    Point_3 split_placement(
        const halfedge_descriptor h,
        const Mesh& mesh) const
    {
        return CGAL::midpoint(
            mesh.point(source(h, mesh)),
            mesh.point(target(h, mesh)));
    }

    void register_split_vertex(
        const vertex_descriptor v,
        const Mesh& mesh)
    {
        double sizing_sum = 0.0;
        std::size_t count = 0;

        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(
                 v, mesh)) {
            const vertex_descriptor n =
                source(h, mesh);
            sizing_sum += sizing_map_[n];
            ++count;
        }

        if (count > 0) {
            sizing_map_[v] =
                sizing_sum /
                static_cast<double>(count);
        } else {
            sizing_map_[v] =
                fallback_target_length_;
        }
    }

private:
    ScalarMap sizing_map_;
    double fallback_target_length_;
};

} // namespace rar
