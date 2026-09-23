#pragma once

#include "rar/Types.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

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
    using Vector_3 = Kernel::Vector_3;
    using ScalarMap =
        Mesh::Property_map<vertex_descriptor, double>;

    explicit CSFSizingField(
        const double mesh_scale,
        Mesh& mesh)
        : mesh_scale_(mesh_scale)
    {
        if (!(mesh_scale_ > 0.0)) {
            throw std::invalid_argument(
                "CSF mesh scale must be > 0");
        }

        raw_curvature_map_ =
            mesh.add_property_map<vertex_descriptor, double>(
                "v:csf_raw_curvature", 0.0).first;
        smoothed_curvature_map_ =
            mesh.add_property_map<vertex_descriptor, double>(
                "v:csf_smoothed_curvature", 0.0).first;
        sizing_map_ =
            mesh.add_property_map<vertex_descriptor, double>(
                "v:csf_target_length", 0.0).first;

        compute_initial_field(mesh);
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
        double sizing_sum = 0.0;
        double raw_sum = 0.0;
        double smooth_sum = 0.0;
        std::size_t count = 0;

        for (const halfedge_descriptor h :
             CGAL::halfedges_around_target(v, mesh)) {
            const vertex_descriptor n = source(h, mesh);
            sizing_sum += sizing_map_[n];
            raw_sum += raw_curvature_map_[n];
            smooth_sum += smoothed_curvature_map_[n];
            ++count;
        }

        if (count == 0) {
            sizing_map_[v] =
                stats_.base_edge_length * mesh_scale_;
            raw_curvature_map_[v] = 0.0;
            smoothed_curvature_map_[v] = 0.0;
            return;
        }

        const double denom = static_cast<double>(count);
        sizing_map_[v] = sizing_sum / denom;
        raw_curvature_map_[v] = raw_sum / denom;
        smoothed_curvature_map_[v] = smooth_sum / denom;
    }

    const CSFFieldStats& stats() const {
        return stats_;
    }

    double raw_curvature(
        const vertex_descriptor v) const
    {
        return raw_curvature_map_[v];
    }

    double smoothed_curvature(
        const vertex_descriptor v) const
    {
        return smoothed_curvature_map_[v];
    }

    double target_length(
        const vertex_descriptor v) const
    {
        return sizing_map_[v];
    }

private:
    struct Neighbor {
        std::size_t vertex = 0;
        double weight = 0.0;
    };

    using Adjacency =
        std::vector<std::vector<Neighbor>>;

    static constexpr double kMinWeight = 1e-4;
    static constexpr double kMinEdgeLength = 1e-12;
    static constexpr int kHistogramBins = 40;

    static double vector_norm(const Vector_3& v) {
        return std::sqrt(
            CGAL::to_double(v.squared_length()));
    }

    static Vector_3 normalized(const Vector_3& v) {
        const double n = vector_norm(v);
        if (!(n > 1e-14) || !std::isfinite(n)) {
            return Vector_3(0.0, 0.0, 0.0);
        }
        return v / n;
    }

    static double edge_length(
        const vertex_descriptor a,
        const vertex_descriptor b,
        const Mesh& mesh)
    {
        return std::sqrt(CGAL::to_double(
            CGAL::squared_distance(
                mesh.point(a), mesh.point(b))));
    }

    static double cotangent_at(
        const Point_3& center,
        const Point_3& a,
        const Point_3& b)
    {
        const Vector_3 u(center, a);
        const Vector_3 v(center, b);
        const Vector_3 cross =
            CGAL::cross_product(u, v);
        const double denom = vector_norm(cross);
        if (!(denom > 1e-14)) {
            return 0.0;
        }
        return CGAL::to_double(u * v) / denom;
    }

    static std::vector<Vector_3> compute_vertex_normals(
        const Mesh& mesh,
        const Mesh::Property_map<
            vertex_descriptor,
            std::size_t>& dense_index)
    {
        std::vector<Vector_3> normals(
            num_vertices(mesh),
            Vector_3(0.0, 0.0, 0.0));

        for (const auto f : faces(mesh)) {
            const halfedge_descriptor h =
                halfedge(f, mesh);
            const vertex_descriptor v0 =
                source(h, mesh);
            const vertex_descriptor v1 =
                target(h, mesh);
            const vertex_descriptor v2 =
                target(next(h, mesh), mesh);

            const Point_3& p0 = mesh.point(v0);
            const Point_3& p1 = mesh.point(v1);
            const Point_3& p2 = mesh.point(v2);

            const Vector_3 face_normal =
                CGAL::cross_product(
                    Vector_3(p0, p1),
                    Vector_3(p0, p2));

            normals[dense_index[v0]] =
                normals[dense_index[v0]] + face_normal;
            normals[dense_index[v1]] =
                normals[dense_index[v1]] + face_normal;
            normals[dense_index[v2]] =
                normals[dense_index[v2]] + face_normal;
        }

        for (Vector_3& n : normals) {
            n = normalized(n);
        }
        return normals;
    }

    static double normal_angle(
        const Vector_3& a,
        const Vector_3& b)
    {
        double cosine =
            CGAL::to_double(a * b);
        cosine =
            (std::max)(-1.0, (std::min)(1.0, cosine));
        return std::acos(cosine);
    }

    static std::vector<double> compute_raw_curvature(
        const Mesh& mesh,
        const Mesh::Property_map<
            vertex_descriptor,
            std::size_t>& dense_index,
        const std::vector<Vector_3>& normals)
    {
        std::vector<double> curvature(
            num_vertices(mesh), 0.0);

        for (const vertex_descriptor v :
             vertices(mesh)) {
            double sum = 0.0;
            std::size_t count = 0;

            for (const halfedge_descriptor h :
                 CGAL::halfedges_around_target(v, mesh)) {
                const vertex_descriptor n =
                    source(h, mesh);
                sum += normal_angle(
                    normals[dense_index[v]],
                    normals[dense_index[n]]);
                ++count;
            }

            if (count > 0) {
                curvature[dense_index[v]] =
                    sum / static_cast<double>(count);
            }
        }

        return curvature;
    }

    static double positive_weight(
        const double value)
    {
        return
            std::isfinite(value) && value > 0.0
                ? value
                : kMinWeight;
    }

    static Adjacency build_weight_graph(
        const Mesh& mesh,
        const Mesh::Property_map<
            vertex_descriptor,
            std::size_t>& dense_index)
    {
        Adjacency adjacency(num_vertices(mesh));

        for (const auto e : edges(mesh)) {
            const halfedge_descriptor h =
                halfedge(e, mesh);
            const vertex_descriptor a =
                source(h, mesh);
            const vertex_descriptor b =
                target(h, mesh);

            double cot_sum = 0.0;

            const auto accumulate_side =
                [&](const halfedge_descriptor side) {
                    if (CGAL::is_border(side, mesh)) {
                        return;
                    }

                    const vertex_descriptor c =
                        target(next(side, mesh), mesh);
                    cot_sum += cotangent_at(
                        mesh.point(c),
                        mesh.point(source(side, mesh)),
                        mesh.point(target(side, mesh)));
                };

            accumulate_side(h);
            accumulate_side(opposite(h, mesh));

            double weight =
                positive_weight(0.5 * cot_sum);

            // Code-faithful CSF variant:
            // cotangent weight additionally divided by edge length.
            weight /=
                (std::max)(
                    edge_length(a, b, mesh),
                    kMinEdgeLength);
            weight = positive_weight(weight);

            const std::size_t ia = dense_index[a];
            const std::size_t ib = dense_index[b];
            adjacency[ia].push_back({ib, weight});
            adjacency[ib].push_back({ia, weight});
        }

        for (auto& neighbors : adjacency) {
            std::sort(
                neighbors.begin(),
                neighbors.end(),
                [](const Neighbor& lhs,
                   const Neighbor& rhs) {
                    return lhs.vertex < rhs.vertex;
                });
        }

        return adjacency;
    }

    static std::vector<std::vector<std::size_t>>
    connected_components(
        const Adjacency& adjacency)
    {
        std::vector<std::vector<std::size_t>>
            components;
        std::vector<char> visited(
            adjacency.size(), 0);

        for (std::size_t start = 0;
             start < adjacency.size();
             ++start) {
            if (visited[start]) {
                continue;
            }

            std::vector<std::size_t> component;
            std::vector<std::size_t> stack{start};
            visited[start] = 1;

            while (!stack.empty()) {
                const std::size_t v = stack.back();
                stack.pop_back();
                component.push_back(v);

                for (const Neighbor& n :
                     adjacency[v]) {
                    if (!visited[n.vertex]) {
                        visited[n.vertex] = 1;
                        stack.push_back(n.vertex);
                    }
                }
            }

            std::sort(
                component.begin(),
                component.end());
            components.push_back(
                std::move(component));
        }

        return components;
    }

    static bool component_is_constant(
        const std::vector<std::size_t>& component,
        const std::vector<double>& values)
    {
        double min_value =
            std::numeric_limits<double>::infinity();
        double max_value =
            -std::numeric_limits<double>::infinity();

        for (const std::size_t v : component) {
            min_value =
                (std::min)(min_value, values[v]);
            max_value =
                (std::max)(max_value, values[v]);
        }

        return max_value - min_value <= 1e-14;
    }

    static std::vector<char> select_fixed_vertices(
        const Adjacency& adjacency,
        const std::vector<double>& values)
    {
        std::vector<char> fixed(
            values.size(), 0);

        for (const auto& component :
             connected_components(adjacency)) {
            if (component.empty()) {
                continue;
            }

            if (component.size() <= 2 ||
                component_is_constant(
                    component, values)) {
                for (const std::size_t v :
                     component) {
                    fixed[v] = 1;
                }
                continue;
            }

            const auto vmax_it =
                std::max_element(
                    component.begin(),
                    component.end(),
                    [&](const std::size_t lhs,
                        const std::size_t rhs) {
                        if (values[lhs] ==
                            values[rhs]) {
                            return lhs > rhs;
                        }
                        return
                            values[lhs] <
                            values[rhs];
                    });
            const std::size_t vmax = *vmax_it;

            std::vector<char> excluded(
                values.size(), 0);
            excluded[vmax] = 1;
            for (const Neighbor& n :
                 adjacency[vmax]) {
                excluded[n.vertex] = 1;
            }

            auto choose_min =
                [&](const bool respect_excluded)
                    -> std::optional<std::size_t> {
                    std::optional<std::size_t> best;
                    double best_value =
                        std::numeric_limits<double>::
                            infinity();

                    for (const std::size_t v :
                         component) {
                        if (v == vmax) {
                            continue;
                        }
                        if (respect_excluded &&
                            excluded[v]) {
                            continue;
                        }

                        if (!best ||
                            values[v] < best_value ||
                            (values[v] == best_value &&
                             v < *best)) {
                            best = v;
                            best_value = values[v];
                        }
                    }
                    return best;
                };

            std::optional<std::size_t> vmin =
                choose_min(true);
            if (!vmin) {
                vmin = choose_min(false);
            }
            if (!vmin) {
                throw std::runtime_error(
                    "CSF anchor selection failed");
            }

            fixed[vmax] = 1;
            fixed[*vmin] = 1;
        }

        return fixed;
    }

    static std::vector<double> solve_csf(
        const Adjacency& adjacency,
        const std::vector<double>& initial)
    {
        if (adjacency.size() != initial.size()) {
            throw std::runtime_error(
                "CSF solve input size mismatch");
        }
        if (initial.empty()) {
            return initial;
        }

        for (const double value : initial) {
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "CSF initial curvature is not finite");
            }
        }

        const std::vector<char> fixed =
            select_fixed_vertices(
                adjacency, initial);

        std::vector<int> free_index(
            initial.size(), -1);
        int free_count = 0;

        for (std::size_t i = 0;
             i < initial.size();
             ++i) {
            if (!fixed[i]) {
                free_index[i] = free_count++;
            }
        }

        if (free_count == 0) {
            return initial;
        }

        using SparseMatrix =
            Eigen::SparseMatrix<double>;
        using Triplet =
            Eigen::Triplet<double>;

        std::vector<Triplet> triplets;
        triplets.reserve(
            static_cast<std::size_t>(
                free_count) * 8);

        Eigen::VectorXd rhs =
            Eigen::VectorXd::Zero(free_count);

        for (std::size_t i = 0;
             i < adjacency.size();
             ++i) {
            const int row = free_index[i];
            if (row < 0) {
                continue;
            }

            double diagonal = 0.0;

            for (const Neighbor& n :
                 adjacency[i]) {
                diagonal += n.weight;
                const int col =
                    free_index[n.vertex];

                if (col >= 0) {
                    triplets.emplace_back(
                        row, col, -n.weight);
                } else {
                    rhs[row] +=
                        n.weight *
                        initial[n.vertex];
                }
            }

            if (!(diagonal > 0.0) ||
                !std::isfinite(diagonal)) {
                throw std::runtime_error(
                    "CSF vertex has no positive "
                    "Laplace diagonal");
            }

            triplets.emplace_back(
                row, row, diagonal);
        }

        SparseMatrix matrix(
            free_count, free_count);
        matrix.setFromTriplets(
            triplets.begin(),
            triplets.end());
        matrix.makeCompressed();

        Eigen::SimplicialLDLT<SparseMatrix>
            solver;
        solver.compute(matrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error(
                "CSF Laplace factorization failed");
        }

        const Eigen::VectorXd solution =
            solver.solve(rhs);
        if (solver.info() != Eigen::Success ||
            !solution.allFinite()) {
            throw std::runtime_error(
                "CSF Laplace solve failed");
        }

        std::vector<double> result = initial;
        for (std::size_t i = 0;
             i < free_index.size();
             ++i) {
            const int idx = free_index[i];
            if (idx >= 0) {
                result[i] = solution[idx];
            }
        }

        return result;
    }

    static std::size_t clamp_index(
        const long long index,
        const std::size_t size)
    {
        if (size == 0) {
            return 0;
        }

        const long long hi =
            static_cast<long long>(size - 1);
        return static_cast<std::size_t>(
            (std::max)(
                0LL,
                (std::min)(index, hi)));
    }

    static std::vector<double>
    histogram_multipliers(
        const std::vector<double>& curvature)
    {
        static constexpr std::array<
            double, 5> multipliers{
                1.8, 1.4, 1.0, 0.8, 0.6};

        std::vector<double> result(
            curvature.size(), 1.0);

        if (curvature.empty()) {
            return result;
        }

        const auto [min_it, max_it] =
            std::minmax_element(
                curvature.begin(),
                curvature.end());
        const double min_value = *min_it;
        const double max_value = *max_it;

        if (std::abs(max_value - min_value) <=
            1e-14) {
            return result;
        }

        const double width =
            (max_value - min_value) /
            static_cast<double>(kHistogramBins);

        std::array<int, kHistogramBins> counts{};
        for (const double value : curvature) {
            int bin = static_cast<int>(
                std::floor(
                    (value - min_value) / width));
            bin = (std::max)(
                0,
                (std::min)(
                    bin, kHistogramBins - 1));
            ++counts[static_cast<std::size_t>(bin)];
        }

        const int max_count =
            *std::max_element(
                counts.begin(), counts.end());

        std::vector<int> max_bins;
        for (int i = 0;
             i < kHistogramBins;
             ++i) {
            if (counts[
                    static_cast<std::size_t>(i)]
                == max_count) {
                max_bins.push_back(i);
            }
        }

        const int base_bin =
            max_bins[max_bins.size() / 2];
        const double base =
            min_value +
            (static_cast<double>(base_bin) + 0.5) *
                width;

        std::vector<double> sorted = curvature;
        std::sort(
            sorted.begin(), sorted.end());

        auto lower =
            std::lower_bound(
                sorted.begin(),
                sorted.end(),
                base);

        long long idx =
            static_cast<long long>(
                std::distance(
                    sorted.begin(), lower));
        idx = (std::max)(
            1LL,
            (std::min)(
                idx,
                static_cast<long long>(
                    sorted.size() - 1)));

        const long long unit_right =
            (std::max)(
                1LL,
                (static_cast<long long>(
                     sorted.size()) -
                 idx) /
                    5LL);

        const std::array<double, 4> cuts{
            sorted[clamp_index(
                (idx / 5LL) * 2LL,
                sorted.size())],
            sorted[clamp_index(
                (idx / 5LL) * 4LL,
                sorted.size())],
            sorted[clamp_index(
                idx + unit_right,
                sorted.size())],
            sorted[clamp_index(
                idx + unit_right * 3LL,
                sorted.size())]
        };

        for (std::size_t i = 0;
             i < curvature.size();
             ++i) {
            const double c = curvature[i];

            if (c < cuts[0]) {
                result[i] = multipliers[0];
            } else if (c < cuts[1]) {
                result[i] = multipliers[1];
            } else if (c < cuts[2]) {
                result[i] = multipliers[2];
            } else if (c < cuts[3]) {
                result[i] = multipliers[3];
            } else {
                result[i] = multipliers[4];
            }
        }

        return result;
    }

    static double mean_edge_length(
        const Mesh& mesh)
    {
        if (num_edges(mesh) == 0) {
            return 0.0;
        }

        double sum = 0.0;
        for (const auto e : edges(mesh)) {
            const halfedge_descriptor h =
                halfedge(e, mesh);
            sum += edge_length(
                source(h, mesh),
                target(h, mesh),
                mesh);
        }

        return
            sum /
            static_cast<double>(num_edges(mesh));
    }

    static void accumulate_stats(
        const std::vector<double>& values,
        double& min_value,
        double& mean_value,
        double& max_value)
    {
        if (values.empty()) {
            min_value = 0.0;
            mean_value = 0.0;
            max_value = 0.0;
            return;
        }

        min_value =
            std::numeric_limits<double>::infinity();
        max_value =
            -std::numeric_limits<double>::infinity();
        double sum = 0.0;

        for (const double value : values) {
            min_value =
                (std::min)(min_value, value);
            max_value =
                (std::max)(max_value, value);
            sum += value;
        }

        mean_value =
            sum /
            static_cast<double>(values.size());
    }

    void compute_initial_field(Mesh& mesh) {
        stats_.vertex_count = num_vertices(mesh);
        stats_.base_edge_length =
            mean_edge_length(mesh);

        if (stats_.vertex_count == 0) {
            return;
        }

        if (!(stats_.base_edge_length > 0.0)) {
            throw std::runtime_error(
                "CSF requires a mesh with "
                "positive edge lengths");
        }

        auto dense_index =
            mesh.add_property_map<
                vertex_descriptor,
                std::size_t>(
                    "v:csf_dense_index", 0).first;

        std::vector<vertex_descriptor>
            dense_vertices;
        dense_vertices.reserve(
            stats_.vertex_count);

        std::size_t next_id = 0;
        for (const vertex_descriptor v :
             vertices(mesh)) {
            dense_index[v] = next_id++;
            dense_vertices.push_back(v);
        }

        const std::vector<Vector_3> normals =
            compute_vertex_normals(
                mesh, dense_index);
        const std::vector<double> raw =
            compute_raw_curvature(
                mesh, dense_index, normals);
        const Adjacency adjacency =
            build_weight_graph(
                mesh, dense_index);
        const std::vector<double> smooth =
            solve_csf(adjacency, raw);
        const std::vector<double> multipliers =
            histogram_multipliers(smooth);

        std::vector<double> sizing(
            stats_.vertex_count, 0.0);

        for (std::size_t i = 0;
             i < dense_vertices.size();
             ++i) {
            const vertex_descriptor v =
                dense_vertices[i];

            raw_curvature_map_[v] = raw[i];
            smoothed_curvature_map_[v] =
                smooth[i];

            const double target =
                (std::max)(
                    kMinEdgeLength,
                    stats_.base_edge_length *
                        mesh_scale_ *
                        multipliers[i]);

            sizing_map_[v] = target;
            sizing[i] = target;
        }

        accumulate_stats(
            raw,
            stats_.raw_curvature_min,
            stats_.raw_curvature_mean,
            stats_.raw_curvature_max);
        accumulate_stats(
            smooth,
            stats_.smoothed_curvature_min,
            stats_.smoothed_curvature_mean,
            stats_.smoothed_curvature_max);
        accumulate_stats(
            sizing,
            stats_.sizing_min,
            stats_.sizing_mean,
            stats_.sizing_max);
    }

private:
    double mesh_scale_;
    ScalarMap raw_curvature_map_;
    ScalarMap smoothed_curvature_map_;
    ScalarMap sizing_map_;
    CSFFieldStats stats_;
};

} // namespace rar
