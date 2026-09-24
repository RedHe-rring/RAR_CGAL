#pragma once

#include "rar/Types.h"

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

namespace rar {

inline bool write_ascii_ply_compact(
    const std::string& path,
    const Mesh& mesh,
    const int precision = 17)
{
    std::ofstream out(path);
    if (!out) {
        return false;
    }

    std::map<Mesh::Vertex_index, std::size_t> vertex_ids;
    std::size_t next_id = 0;

    for (const auto v : vertices(mesh)) {
        vertex_ids.emplace(v, next_id++);
    }

    out << std::setprecision(precision);
    out
        << "ply\n"
        << "format ascii 1.0\n"
        << "comment RAR_CGAL compact final mesh export\n"
        << "element vertex " << vertex_ids.size() << "\n"
        << "property double x\n"
        << "property double y\n"
        << "property double z\n"
        << "element face " << num_faces(mesh) << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (const auto v : vertices(mesh)) {
        const auto& p = mesh.point(v);
        out
            << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z()) << '\n';
    }

    for (const auto f : faces(mesh)) {
        std::vector<std::size_t> ids;

        for (const auto v :
             CGAL::vertices_around_face(
                 halfedge(f, mesh), mesh)) {
            const auto it = vertex_ids.find(v);
            if (it == vertex_ids.end()) {
                return false;
            }
            ids.push_back(it->second);
        }

        if (ids.size() < 3 || ids.size() > 255) {
            return false;
        }

        out << ids.size();
        for (const std::size_t id : ids) {
            if (id >= vertex_ids.size()) {
                return false;
            }
            out << ' ' << id;
        }
        out << '\n';
    }

    return static_cast<bool>(out);
}

inline bool has_ply_extension(
    const std::string& path)
{
    std::string extension =
        std::filesystem::path(path).extension().string();

    std::transform(
        extension.begin(),
        extension.end(),
        extension.begin(),
        [](const unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });

    return extension == ".ply";
}

} // namespace rar
