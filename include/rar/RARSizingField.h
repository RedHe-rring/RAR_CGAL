#pragma once

#include "rar/Types.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Weights/mixed_voronoi_region_weights.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace rar {

struct RARFieldStats {
    double curvature_min = 0.0;
    double curvature_max = 0.0;
    double curvature_mean = 0.0;
    double sizing_min = 0.0;
    double sizing_max = 0.0;
    double sizing_mean = 0.0;
    std::size_t vertex_count = 0;
};

inline double rar_target_length(
    const double kappa,
    const double epsilon,
    const double min_edge_length,
    const double max_edge_length)
{
    if (!(epsilon > 0.0)) {
        throw std::invalid_argument("RAR epsilon must be > 0");
    }
    if (!(min_edge_length > 0.0) || min_edge_length > max_edge_length) {
        throw std::invalid_argument("Invalid RAR edge-length bounds");
    }

    if (!std::isfinite(kappa) || kappa <= 1e-15) {
        return max_edge_length;
    }

    const double value_sq =
        6.0 * epsilon / kappa - 3.0 * epsilon * epsilon;

    if (!std::isfinite(value_sq) ||
        value_sq >= max_edge_length * max_edge_length) {
        return max_edge_length;
    }

    if (value_sq <= min_edge_length * min_edge_length) {
        return min_edge_length;
    }

    return std::sqrt(value_sq);
}

class RARSizingField {
public:
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
    using FT = Kernel::FT;
    using Point_3 = Point;
    using Vector_3 = Kernel::Vector_3;
    using ScalarMap = Mesh::Property_map<vertex_descriptor, double>;

    RARSizingField(
        const double epsilon,
        const std::pair<double, double>& edge_len_min_max,
        Mesh& mesh)
        : epsilon_(epsilon),
          min_edge_length_(edge_len_min_max.first),
          max_edge_length_(edge_len_min_max.second)
    {
        if (!(epsilon_ > 0.0)) {
            throw std::invalid_argument("RAR epsilon must be > 0");
        }
        if (!(min_edge_length_ > 0.0) ||
            min_edge_length_ > max_edge_length_) {
            throw std::invalid_argument("Invalid RAR edge-length bounds");
        }

        curvature_map_ =
            mesh.add_property_map<vertex_descriptor, double>(
                "v:rar_curvature", 0.0).first;
        sizing_map_ =
            mesh.add_property_map<vertex_descriptor, double>(
                "v:rar_target_length", max_edge_length_).first;

        compute_initial_field(mesh);
    }

    FT at(const vertex_descriptor v, const Mesh&) const {
        return FT(sizing_map_[v]);
    }

    std::optional<FT> is_too_long(
        const vertex_descriptor va,
        const vertex_descriptor vb,
        const Mesh& mesh) const
    {
        const double edge_sq =
            CGAL::to_double(CGAL::squared_distance(
                mesh.point(va), mesh.point(vb)));

        const double target =
            (4.0 / 3.0) *
            (std::min)(sizing_map_[va], sizing_map_[vb]);
        const double target_sq = target * target;

        if (edge_sq > target_sq) {
            return FT(edge_sq / target_sq);
        }
        return std::nullopt;
    }

    std::optional<FT> is_too_short(
        const halfedge_descriptor h,
        const Mesh& mesh) const
    {
        const vertex_descriptor va = source(h, mesh);
        const vertex_descriptor vb = target(h, mesh);

        const double edge_sq =
            CGAL::to_double(CGAL::squared_distance(
                mesh.point(va), mesh.point(vb)));

        const double target =
            (4.0 / 5.0) *
            (std::min)(sizing_map_[va], sizing_map_[vb]);
        const double target_sq = target * target;

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
        double sum = 0.0;
        std::size_t count = 0;

        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(v, mesh)) {
            sum += sizing_map_[source(h, mesh)];
            ++count;
        }

        sizing_map_[v] =
            count > 0 ? sum / static_cast<double>(count)
                      : max_edge_length_;
        curvature_map_[v] = 0.0;
    }

    const RARFieldStats& stats() const {
        return stats_;
    }

    double curvature(const vertex_descriptor v) const {
        return curvature_map_[v];
    }

    double target_length(const vertex_descriptor v) const {
        return sizing_map_[v];
    }

private:
    static double vector_norm(const Vector_3& v) {
        return std::sqrt(CGAL::to_double(v.squared_length()));
    }

    static double cotangent_at(
        const Point_3& p,
        const Point_3& q,
        const Point_3& r)
    {
        const Vector_3 u(q, p);
        const Vector_3 v(q, r);
        const Vector_3 cross = CGAL::cross_product(u, v);
        const double denom = vector_norm(cross);

        if (!(denom > 1e-20)) {
            return 0.0;
        }

        return CGAL::to_double(u * v) / denom;
    }

    static double angle_at(
        const Point_3& p,
        const Point_3& q,
        const Point_3& r)
    {
        const Vector_3 u(q, p);
        const Vector_3 v(q, r);
        const double nu = vector_norm(u);
        const double nv = vector_norm(v);

        if (!(nu > 1e-20) || !(nv > 1e-20)) {
            return 0.0;
        }

        double cosine = CGAL::to_double(u * v) / (nu * nv);
        cosine = (std::max)(-1.0, (std::min)(1.0, cosine));
        return std::acos(cosine);
    }

    static bool is_boundary_vertex(
        const vertex_descriptor v,
        const Mesh& mesh)
    {
        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(v, mesh)) {
            if (CGAL::is_border(h, mesh)) {
                return true;
            }
        }
        return false;
    }

    double vertex_curvature(
        const vertex_descriptor v,
        const Mesh& mesh) const
    {
        const Point_3& pi = mesh.point(v);

        double mixed_area = 0.0;
        double angle_sum = 0.0;

        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(v, mesh)) {
            if (CGAL::is_border(h, mesh)) {
                continue;
            }

            const vertex_descriptor vj = source(h, mesh);
            const vertex_descriptor vk = target(next(h, mesh), mesh);
            const Point_3& pj = mesh.point(vj);
            const Point_3& pk = mesh.point(vk);

            const double local_area = CGAL::to_double(
                CGAL::Weights::mixed_voronoi_area(pi, pj, pk));
            if (std::isfinite(local_area) && local_area > 0.0) {
                mixed_area += local_area;
            }

            angle_sum += angle_at(pj, pi, pk);
        }

        if (!(mixed_area > 1e-20) || !std::isfinite(mixed_area)) {
            return 0.0;
        }

        Vector_3 laplace(0.0, 0.0, 0.0);

        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(v, mesh)) {
            const vertex_descriptor vj = source(h, mesh);
            double cot_sum = 0.0;

            if (!CGAL::is_border(h, mesh)) {
                const vertex_descriptor vk =
                    target(next(h, mesh), mesh);
                cot_sum += cotangent_at(
                    mesh.point(v), mesh.point(vk), mesh.point(vj));
            }

            const halfedge_descriptor ho = opposite(h, mesh);
            if (!CGAL::is_border(ho, mesh)) {
                const vertex_descriptor vk =
                    target(next(ho, mesh), mesh);
                cot_sum += cotangent_at(
                    mesh.point(v), mesh.point(vk), mesh.point(vj));
            }

            laplace =
                laplace +
                Vector_3(pi, mesh.point(vj)) * cot_sum;
        }

        laplace = laplace / (2.0 * mixed_area);

        // RAR Eq. (3): H_i = 1/2 ||Delta x_i||.
        const double H = 0.5 * vector_norm(laplace);

        const double pi_const = 3.14159265358979323846;
        const double angle_budget =
            is_boundary_vertex(v, mesh) ? pi_const
                                        : 2.0 * pi_const;
        const double K = (angle_budget - angle_sum) / mixed_area;

        // RAR Eq. (5). Clamp the discriminant against numerical noise.
        const double discriminant =
            (std::max)(0.0, H * H - K);
        const double kappa =
            H + std::sqrt(discriminant);

        return std::isfinite(kappa) ? kappa : 0.0;
    }

    void compute_initial_field(Mesh& mesh) {
        stats_.vertex_count = num_vertices(mesh);

        if (stats_.vertex_count == 0) {
            return;
        }

        double curvature_sum = 0.0;
        double sizing_sum = 0.0;
        stats_.curvature_min =
            std::numeric_limits<double>::infinity();
        stats_.curvature_max = 0.0;
        stats_.sizing_min =
            std::numeric_limits<double>::infinity();
        stats_.sizing_max = 0.0;

        for (const vertex_descriptor v : vertices(mesh)) {
            const double kappa = vertex_curvature(v, mesh);
            const double target = rar_target_length(
                kappa,
                epsilon_,
                min_edge_length_,
                max_edge_length_);

            curvature_map_[v] = kappa;
            sizing_map_[v] = target;

            stats_.curvature_min =
                (std::min)(stats_.curvature_min, kappa);
            stats_.curvature_max =
                (std::max)(stats_.curvature_max, kappa);
            stats_.sizing_min =
                (std::min)(stats_.sizing_min, target);
            stats_.sizing_max =
                (std::max)(stats_.sizing_max, target);

            curvature_sum += kappa;
            sizing_sum += target;
        }

        const double count =
            static_cast<double>(stats_.vertex_count);
        stats_.curvature_mean = curvature_sum / count;
        stats_.sizing_mean = sizing_sum / count;
    }

private:
    double epsilon_;
    double min_edge_length_;
    double max_edge_length_;
    ScalarMap curvature_map_;
    ScalarMap sizing_map_;
    RARFieldStats stats_;
};

} // namespace rar
