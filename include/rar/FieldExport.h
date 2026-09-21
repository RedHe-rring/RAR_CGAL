#pragma once

#include "rar/Types.h"

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>

namespace rar {

template <typename CurvatureFn, typename TargetFn>
void write_field_diagnostics(
    const Mesh& mesh,
    const std::string& prefix,
    CurvatureFn curvature_fn,
    TargetFn target_fn)
{
    if (prefix.empty()) {
        return;
    }

    const std::filesystem::path prefix_path(prefix);
    if (prefix_path.has_parent_path()) {
        std::filesystem::create_directories(
            prefix_path.parent_path());
    }

    const std::string csv_path = prefix + "_field.csv";
    const std::string ply_path = prefix + "_field.ply";

    std::map<Mesh::Vertex_index, std::size_t> vertex_ids;
    std::size_t next_id = 0;
    for (const auto v : vertices(mesh)) {
        vertex_ids.emplace(v, next_id++);
    }

    std::ofstream csv(csv_path);
    if (!csv) {
        throw std::runtime_error(
            "Failed to open field CSV: " + csv_path);
    }

    csv << std::setprecision(17);
    csv << "vertex_id,x,y,z,curvature,target_length\n";
    for (const auto v : vertices(mesh)) {
        const auto& p = mesh.point(v);
        csv
            << vertex_ids.at(v) << ','
            << CGAL::to_double(p.x()) << ','
            << CGAL::to_double(p.y()) << ','
            << CGAL::to_double(p.z()) << ','
            << curvature_fn(v) << ','
            << target_fn(v) << '\n';
    }

    std::ofstream ply(ply_path);
    if (!ply) {
        throw std::runtime_error(
            "Failed to open field PLY: " + ply_path);
    }

    ply << std::setprecision(17);
    ply
        << "ply\n"
        << "format ascii 1.0\n"
        << "comment RAR_CGAL initial sizing-field diagnostics\n"
        << "element vertex " << num_vertices(mesh) << "\n"
        << "property double x\n"
        << "property double y\n"
        << "property double z\n"
        << "property double curvature\n"
        << "property double target_length\n"
        << "element face " << num_faces(mesh) << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (const auto v : vertices(mesh)) {
        const auto& p = mesh.point(v);
        ply
            << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z()) << ' '
            << curvature_fn(v) << ' '
            << target_fn(v) << '\n';
    }

    for (const auto f : faces(mesh)) {
        std::size_t count = 0;
        for (const auto v :
             CGAL::vertices_around_face(halfedge(f, mesh), mesh)) {
            (void)v;
            ++count;
        }

        ply << count;
        for (const auto v :
             CGAL::vertices_around_face(halfedge(f, mesh), mesh)) {
            ply << ' ' << vertex_ids.at(v);
        }
        ply << '\n';
    }
}

} // namespace rar
