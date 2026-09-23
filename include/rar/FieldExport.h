#pragma once

#include "rar/Types.h"

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rar {

namespace detail {

using VertexScalarFn =
    std::function<double(Mesh::Vertex_index)>;

inline void write_scalar_field_diagnostics(
    const Mesh& mesh,
    const std::string& prefix,
    const std::vector<
        std::pair<std::string, VertexScalarFn>>&
        scalar_fields)
{
    if (prefix.empty()) {
        return;
    }

    const std::filesystem::path prefix_path(prefix);
    if (prefix_path.has_parent_path()) {
        std::filesystem::create_directories(
            prefix_path.parent_path());
    }

    const std::string csv_path = prefix + ".csv";
    const std::string ply_path = prefix + ".ply";

    std::map<Mesh::Vertex_index, std::size_t>
        vertex_ids;
    std::size_t next_id = 0;
    for (const auto v : vertices(mesh)) {
        vertex_ids.emplace(v, next_id++);
    }

    std::ofstream csv(csv_path);
    if (!csv) {
        throw std::runtime_error(
            "Failed to open field CSV: " +
            csv_path);
    }

    csv << std::setprecision(17);
    csv << "vertex_id,x,y,z";
    for (const auto& field : scalar_fields) {
        csv << ',' << field.first;
    }
    csv << '\n';

    for (const auto v : vertices(mesh)) {
        const auto& p = mesh.point(v);
        csv
            << vertex_ids.at(v) << ','
            << CGAL::to_double(p.x()) << ','
            << CGAL::to_double(p.y()) << ','
            << CGAL::to_double(p.z());

        for (const auto& field : scalar_fields) {
            csv << ',' << field.second(v);
        }
        csv << '\n';
    }

    std::ofstream ply(ply_path);
    if (!ply) {
        throw std::runtime_error(
            "Failed to open field PLY: " +
            ply_path);
    }

    ply << std::setprecision(17);
    ply
        << "ply\n"
        << "format ascii 1.0\n"
        << "comment RAR_CGAL initial sizing-field diagnostics\n"
        << "element vertex "
        << num_vertices(mesh) << "\n"
        << "property double x\n"
        << "property double y\n"
        << "property double z\n";

    for (const auto& field : scalar_fields) {
        ply
            << "property double "
            << field.first << '\n';
    }

    ply
        << "element face "
        << num_faces(mesh) << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (const auto v : vertices(mesh)) {
        const auto& p = mesh.point(v);
        ply
            << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z());

        for (const auto& field : scalar_fields) {
            ply << ' ' << field.second(v);
        }
        ply << '\n';
    }

    for (const auto f : faces(mesh)) {
        std::size_t count = 0;
        for (const auto v :
             CGAL::vertices_around_face(
                 halfedge(f, mesh), mesh)) {
            (void)v;
            ++count;
        }

        ply << count;
        for (const auto v :
             CGAL::vertices_around_face(
                 halfedge(f, mesh), mesh)) {
            ply << ' ' << vertex_ids.at(v);
        }
        ply << '\n';
    }
}

} // namespace detail

template <typename CurvatureFn, typename TargetFn>
void write_field_diagnostics(
    const Mesh& mesh,
    const std::string& prefix,
    CurvatureFn curvature_fn,
    TargetFn target_fn)
{
    detail::write_scalar_field_diagnostics(
        mesh,
        prefix,
        {
            {
                "curvature",
                detail::VertexScalarFn(
                    curvature_fn)
            },
            {
                "target_length",
                detail::VertexScalarFn(
                    target_fn)
            }
        });
}

template <
    typename RawCurvatureFn,
    typename CurvatureFn,
    typename TargetFn>
void write_csf_field_diagnostics(
    const Mesh& mesh,
    const std::string& prefix,
    RawCurvatureFn raw_curvature_fn,
    CurvatureFn curvature_fn,
    TargetFn target_fn)
{
    detail::write_scalar_field_diagnostics(
        mesh,
        prefix,
        {
            {
                "raw_curvature",
                detail::VertexScalarFn(
                    raw_curvature_fn)
            },
            {
                "curvature",
                detail::VertexScalarFn(
                    curvature_fn)
            },
            {
                "target_length",
                detail::VertexScalarFn(
                    target_fn)
            }
        });
}

} // namespace rar
