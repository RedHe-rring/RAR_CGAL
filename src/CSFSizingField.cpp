#include "rar/CSFSizingField.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rar {
namespace {

using Vertex = CSFSizingField::vertex_descriptor;
using Halfedge = CSFSizingField::halfedge_descriptor;
using Vector = Kernel::Vector_3;
using DenseMap = Mesh::Property_map<Vertex, std::size_t>;

constexpr double kMinEdgeLength = 1e-12;
constexpr int kPrimarySmoothSteps = 3;
constexpr int kNeighborAverageSteps = 3;
constexpr double kLambda = 0.5;
constexpr int kHistogramBins = 40;
constexpr std::array<double, 5> kMultipliers{
    1.8, 1.4, 1.0, 0.8, 0.6};

struct AuthorTopology {
    std::vector<std::array<std::size_t, 2>> face_pairs_flat;
    std::vector<std::vector<std::size_t>> neighbor_entries;
};

double vector_norm(const Vector& v) {
    return std::sqrt(CGAL::to_double(v.squared_length()));
}

Vector unit_vector(const Vector& v) {
    const double n = vector_norm(v);
    if (!(n > 0.0) || !std::isfinite(n)) {
        return Vector(0.0, 0.0, 0.0);
    }
    return v / n;
}

bool is_zero_vector(const Vector& v) {
    return CGAL::to_double(v.squared_length()) == 0.0;
}

double edge_length(Vertex a, Vertex b, const Mesh& mesh) {
    return std::sqrt(CGAL::to_double(
        CGAL::squared_distance(mesh.point(a), mesh.point(b))));
}

double triangle_area(
    const Point& p0,
    const Point& p1,
    const Point& p2)
{
    return 0.5 * vector_norm(
        CGAL::cross_product(Vector(p0, p1), Vector(p0, p2)));
}

AuthorTopology build_author_topology(
    const Mesh& mesh,
    const DenseMap& dense_index)
{
    AuthorTopology topology;
    topology.neighbor_entries.resize(num_vertices(mesh));

    for (const auto f : faces(mesh)) {
        const Halfedge h = halfedge(f, mesh);
        const Vertex v0 = source(h, mesh);
        const Vertex v1 = target(h, mesh);
        const Vertex v2 = target(next(h, mesh), mesh);

        const std::size_t i0 = dense_index[v0];
        const std::size_t i1 = dense_index[v1];
        const std::size_t i2 = dense_index[v2];

        // Match the author's pointNeighbor construction exactly:
        // v0 -> (v1,v2), v1 -> (v2,v0), v2 -> (v0,v1).
        topology.neighbor_entries[i0].push_back(i1);
        topology.neighbor_entries[i0].push_back(i2);

        topology.neighbor_entries[i1].push_back(i2);
        topology.neighbor_entries[i1].push_back(i0);

        topology.neighbor_entries[i2].push_back(i0);
        topology.neighbor_entries[i2].push_back(i1);
    }

    return topology;
}

std::vector<std::size_t> unique_neighbors_in_author_order(
    const std::vector<std::size_t>& entries)
{
    std::vector<std::size_t> result;
    result.reserve(entries.size());

    for (const std::size_t v : entries) {
        if (std::find(result.begin(), result.end(), v) == result.end()) {
            result.push_back(v);
        }
    }
    return result;
}

std::vector<Vector> compute_author_vertex_normals(
    const Mesh& mesh,
    const DenseMap& dense_index)
{
    std::vector<Vector> weighted_sum(
        num_vertices(mesh), Vector(0.0, 0.0, 0.0));
    std::vector<double> area_sum(num_vertices(mesh), 0.0);

    for (const auto f : faces(mesh)) {
        const Halfedge h = halfedge(f, mesh);
        const Vertex v0 = source(h, mesh);
        const Vertex v1 = target(h, mesh);
        const Vertex v2 = target(next(h, mesh), mesh);

        const Point& p0 = mesh.point(v0);
        const Point& p1 = mesh.point(v1);
        const Point& p2 = mesh.point(v2);

        const Vector face_normal = unit_vector(
            CGAL::cross_product(Vector(p0, p1), Vector(p0, p2)));
        const double area = triangle_area(p0, p1, p2);

        for (const Vertex v : {v0, v1, v2}) {
            const std::size_t i = dense_index[v];
            weighted_sum[i] = weighted_sum[i] + face_normal * area;
            area_sum[i] += area;
        }
    }

    // Important: the author's MeshGeometric_Normal_Point() returns
    // the area-weighted average directly. It does NOT renormalize
    // the final vertex normal.
    for (std::size_t i = 0; i < weighted_sum.size(); ++i) {
        if (area_sum[i] > 0.0) {
            weighted_sum[i] = weighted_sum[i] / area_sum[i];
        }
    }

    return weighted_sum;
}

double author_normal_angle(const Vector& a, const Vector& b) {
    double cosine = CGAL::to_double(a * b);
    cosine = (std::max)(-1.0, (std::min)(1.0, cosine));
    return std::acos(cosine);
}

double author_inner_angle(
    std::size_t center,
    std::size_t p2,
    std::size_t p3,
    const std::vector<Vertex>& vertices_dense,
    const Mesh& mesh)
{
    const Point& pc = mesh.point(vertices_dense[center]);
    const Point& pa = mesh.point(vertices_dense[p2]);
    const Point& pb = mesh.point(vertices_dense[p3]);

    Vector a(pc, pa);
    Vector b(pc, pb);

    const double an = vector_norm(a);
    const double bn = vector_norm(b);
    if (!(an > 0.0) || !(bn > 0.0)) {
        return 0.0;
    }

    a = a / an;
    b = b / bn;

    double cosine = CGAL::to_double(a * b);
    cosine = (std::max)(-1.0, (std::min)(1.0, cosine));
    return std::acos(cosine);
}

std::vector<double> compute_author_raw_curvature(
    const AuthorTopology& topology,
    const std::vector<Vector>& vertex_normals)
{
    std::vector<double> raw(vertex_normals.size(), 0.0);

    for (std::size_t i = 0; i < raw.size(); ++i) {
        const auto& entries = topology.neighbor_entries[i];
        if (entries.empty()) {
            continue;
        }

        double angle_sum = 0.0;

        for (std::size_t j = 0; j + 1 < entries.size(); j += 2) {
            const std::size_t p2 = entries[j];
            const std::size_t p3 = entries[j + 1];

            Vector p2n = vertex_normals[p2];
            Vector p3n = vertex_normals[p3];

            // Match the author's zero-normal fallback.
            if (is_zero_vector(p2n)) {
                p2n = p3n;
            }
            if (is_zero_vector(p3n)) {
                p3n = p2n;
            }

            angle_sum += author_normal_angle(vertex_normals[i], p2n);
            angle_sum += author_normal_angle(vertex_normals[i], p3n);
        }

        // The original code divides by pointNeighbor[i].size(),
        // i.e. by the flattened face-pair entry count.
        raw[i] = angle_sum / static_cast<double>(entries.size());
    }

    return raw;
}

struct AuthorWeights {
    std::vector<std::vector<std::size_t>> neighbors;
    std::vector<std::vector<double>> normalized;
};

AuthorWeights build_author_csf_weights(
    const AuthorTopology& topology,
    const std::vector<Vertex>& vertices_dense,
    const Mesh& mesh)
{
    AuthorWeights result;
    result.neighbors.resize(vertices_dense.size());
    result.normalized.resize(vertices_dense.size());

    for (std::size_t i = 0; i < vertices_dense.size(); ++i) {
        const auto& entries = topology.neighbor_entries[i];
        auto& neighbors = result.neighbors[i];
        neighbors = unique_neighbors_in_author_order(entries);

        std::vector<double> cotangent_weights(neighbors.size(), 0.0);
        std::vector<double> distances(neighbors.size(), 0.0);

        for (std::size_t j = 0; j < neighbors.size(); ++j) {
            const std::size_t p2 = neighbors[j];

            distances[j] = std::sqrt(CGAL::to_double(
                CGAL::squared_distance(
                    mesh.point(vertices_dense[i]),
                    mesh.point(vertices_dense[p2]))));

            for (std::size_t k = 0; k + 1 < entries.size(); k += 2) {
                const std::size_t p21 = entries[k];
                const std::size_t p22 = entries[k + 1];

                if (p21 != p2 && p22 != p2) {
                    continue;
                }

                const std::size_t opposite =
                    (p21 == p2) ? p22 : p21;

                const double angle = author_inner_angle(
                    opposite,
                    i,
                    p2,
                    vertices_dense,
                    mesh);

                double tan_value = std::tan(std::abs(angle));
                if (tan_value < 0.1) {
                    tan_value = 0.1;
                }
                if (tan_value > 10.0) {
                    tan_value = 10.0;
                }

                cotangent_weights[j] += 1.0 / tan_value;
            }
        }

        double sum_weight = 0.0;
        for (std::size_t j = 0; j < neighbors.size(); ++j) {
            const double distance =
                (std::max)(distances[j], kMinEdgeLength);
            sum_weight += cotangent_weights[j] / distance;
        }

        auto& normalized = result.normalized[i];
        normalized.assign(neighbors.size(), 0.0);

        if (sum_weight > 0.0 && std::isfinite(sum_weight)) {
            for (std::size_t j = 0; j < neighbors.size(); ++j) {
                const double distance =
                    (std::max)(distances[j], kMinEdgeLength);
                normalized[j] =
                    (cotangent_weights[j] / distance) / sum_weight;
            }
        }
    }

    return result;
}

std::vector<double> author_primary_csf_smoothing(
    const std::vector<double>& raw,
    const AuthorWeights& weights)
{
    std::vector<double> smooth = raw;

    // Match MeshGeometric_Harmonic_N_Value(): three in-place
    // Gauss-Seidel-like sweeps with lambda = 0.5.
    for (int step = 0; step < kPrimarySmoothSteps; ++step) {
        for (std::size_t i = 0; i < smooth.size(); ++i) {
            const auto& neighbors = weights.neighbors[i];
            const auto& normalized = weights.normalized[i];

            if (neighbors.empty()) {
                continue;
            }

            double neighbor_value = 0.0;
            for (std::size_t j = 0; j < neighbors.size(); ++j) {
                neighbor_value +=
                    normalized[j] * smooth[neighbors[j]];
            }

            smooth[i] =
                smooth[i] * kLambda +
                neighbor_value * (1.0 - kLambda);
        }
    }

    return smooth;
}

std::vector<double> author_secondary_neighbor_average(
    const std::vector<double>& input,
    const AuthorTopology& topology,
    const std::vector<Vertex>& vertices_dense,
    const Mesh& mesh)
{
    std::vector<double> values = input;

    // Match MeshOptimization_Apt_L_init(): three synchronous
    // neighbor averages, weighted by edge length.
    for (int step = 0; step < kNeighborAverageSteps; ++step) {
        const std::vector<double> old = values;

        for (std::size_t i = 0; i < values.size(); ++i) {
            const std::vector<std::size_t> neighbors =
                unique_neighbors_in_author_order(
                    topology.neighbor_entries[i]);

            if (neighbors.empty()) {
                continue;
            }

            double weight_sum = 0.0;
            std::vector<double> weights(neighbors.size(), 0.0);

            for (std::size_t j = 0; j < neighbors.size(); ++j) {
                weights[j] = std::sqrt(CGAL::to_double(
                    CGAL::squared_distance(
                        mesh.point(vertices_dense[i]),
                        mesh.point(vertices_dense[neighbors[j]]))));
                weight_sum += weights[j];
            }

            if (!(weight_sum > 0.0)) {
                continue;
            }

            double averaged = 0.0;
            for (std::size_t j = 0; j < neighbors.size(); ++j) {
                averaged +=
                    (weights[j] / weight_sum) *
                    old[neighbors[j]];
            }

            values[i] = averaged;
        }
    }

    return values;
}

std::vector<double> author_histogram_multipliers(
    const std::vector<double>& curvature)
{
    std::vector<double> result(curvature.size(), 1.0);

    if (curvature.empty()) {
        return result;
    }

    std::vector<double> sorted = curvature;
    std::sort(sorted.begin(), sorted.end());

    const double min_value = sorted.front();
    const double max_value = sorted.back();
    const double unit_step =
        (max_value - min_value) /
        static_cast<double>(kHistogramBins);

    std::array<double, kHistogramBins> histogram{};

    // Match the author's inclusive [lower, upper] tests.
    for (int i = 0; i < kHistogramBins; ++i) {
        const double lower =
            min_value + unit_step * static_cast<double>(i);
        const double upper =
            min_value + unit_step * static_cast<double>(i + 1);

        for (const double value : sorted) {
            if (value <= upper && value >= lower) {
                histogram[static_cast<std::size_t>(i)] += 1.0;
            }
        }
    }

    int his_max = 0;
    double index_sum = -1.0;
    for (int i = 0; i < kHistogramBins; ++i) {
        if (histogram[static_cast<std::size_t>(i)] > index_sum) {
            his_max = i;
            index_sum = histogram[static_cast<std::size_t>(i)];
        }
    }

    const double lower_mid =
        min_value + unit_step * static_cast<double>(his_max);
    const double upper_mid =
        min_value + unit_step * static_cast<double>(his_max + 1);
    const double c_base = (lower_mid + upper_mid) / 2.0;

    std::size_t index_c_base = 0;
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i - 1] <= c_base && c_base <= sorted[i]) {
            index_c_base = i;
            break;
        }
    }

    const std::size_t max_cu = sorted.size();
    const std::size_t unit_right =
        (max_cu - index_c_base) / 5;

    const auto safe_index =
        [max_cu](std::size_t i) {
            return (std::min)(i, max_cu - 1);
        };

    const double li1 = sorted[safe_index((index_c_base / 5) * 2)];
    const double li2 = sorted[safe_index((index_c_base / 5) * 4)];
    const double li3 = sorted[safe_index(index_c_base + unit_right)];
    const double li4 = sorted[safe_index(index_c_base + unit_right * 3)];

    for (std::size_t i = 0; i < curvature.size(); ++i) {
        const double c = curvature[i];

        if (c < li1) {
            result[i] = kMultipliers[0];
        } else if (c < li2 && c >= li1) {
            result[i] = kMultipliers[1];
        } else if (c < li3 && c >= li2) {
            result[i] = kMultipliers[2];
        } else if (c < li4 && c >= li3) {
            result[i] = kMultipliers[3];
        } else {
            result[i] = kMultipliers[4];
        }
    }

    return result;
}

double mean_edge_length(const Mesh& mesh) {
    if (num_edges(mesh) == 0) {
        return 0.0;
    }

    double sum = 0.0;
    for (const auto e : edges(mesh)) {
        const Halfedge h = halfedge(e, mesh);
        sum += edge_length(source(h, mesh), target(h, mesh), mesh);
    }

    return sum / static_cast<double>(num_edges(mesh));
}

void accumulate_stats(
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

    min_value = std::numeric_limits<double>::infinity();
    max_value = -std::numeric_limits<double>::infinity();
    double sum = 0.0;

    for (const double value : values) {
        min_value = (std::min)(min_value, value);
        max_value = (std::max)(max_value, value);
        sum += value;
    }

    mean_value = sum / static_cast<double>(values.size());
}

} // namespace

CSFSizingField::CSFSizingField(
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

CSFSizingField::FT CSFSizingField::at(
    const vertex_descriptor v,
    const Mesh&) const
{
    return FT(sizing_map_[v]);
}

std::optional<CSFSizingField::FT>
CSFSizingField::is_too_long(
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

std::optional<CSFSizingField::FT>
CSFSizingField::is_too_short(
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

CSFSizingField::Point_3 CSFSizingField::split_placement(
    const halfedge_descriptor h,
    const Mesh& mesh) const
{
    return CGAL::midpoint(
        mesh.point(source(h, mesh)),
        mesh.point(target(h, mesh)));
}

void CSFSizingField::register_split_vertex(
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

const CSFFieldStats& CSFSizingField::stats() const {
    return stats_;
}

double CSFSizingField::raw_curvature(
    const vertex_descriptor v) const
{
    return raw_curvature_map_[v];
}

double CSFSizingField::smoothed_curvature(
    const vertex_descriptor v) const
{
    return smoothed_curvature_map_[v];
}

double CSFSizingField::target_length(
    const vertex_descriptor v) const
{
    return sizing_map_[v];
}

void CSFSizingField::compute_initial_field(Mesh& mesh) {
    stats_.vertex_count = num_vertices(mesh);
    stats_.base_edge_length = mean_edge_length(mesh);

    if (stats_.vertex_count == 0) {
        return;
    }

    if (!(stats_.base_edge_length > 0.0)) {
        throw std::runtime_error(
            "CSF requires a mesh with positive edge lengths");
    }

    auto dense_index =
        mesh.add_property_map<vertex_descriptor, std::size_t>(
            "v:csf_dense_index", 0).first;

    std::vector<vertex_descriptor> dense_vertices;
    dense_vertices.reserve(stats_.vertex_count);

    std::size_t next_id = 0;
    for (const vertex_descriptor v : vertices(mesh)) {
        dense_index[v] = next_id++;
        dense_vertices.push_back(v);
    }

    const AuthorTopology topology =
        build_author_topology(mesh, dense_index);

    const std::vector<Vector> vertex_normals =
        compute_author_vertex_normals(mesh, dense_index);

    const std::vector<double> raw =
        compute_author_raw_curvature(topology, vertex_normals);

    const AuthorWeights weights =
        build_author_csf_weights(
            topology, dense_vertices, mesh);

    const std::vector<double> primary_smooth =
        author_primary_csf_smoothing(raw, weights);

    const std::vector<double> smooth =
        author_secondary_neighbor_average(
            primary_smooth,
            topology,
            dense_vertices,
            mesh);

    const std::vector<double> multipliers =
        author_histogram_multipliers(smooth);

    std::vector<double> sizing(stats_.vertex_count, 0.0);

    for (std::size_t i = 0; i < dense_vertices.size(); ++i) {
        const vertex_descriptor v = dense_vertices[i];

        raw_curvature_map_[v] = raw[i];
        smoothed_curvature_map_[v] = smooth[i];

        // The original code first assigns multiplier * L_ave,
        // then multiplies all local targets by meshScale.
        // This is algebraically identical.
        const double target = (std::max)(
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

} // namespace rar
