#pragma once
// C1 external sizing adapter for the existing CGAL remeshing backend.
// CSV is tied to the original mesh by ordered vertex ID and coordinates.
// This does not certify a Hausdorff error bound; inserted vertices inherit
// an interpolated sizing value via PropertyMapSizingField.

#include "rar/FeatureConstraints.h"
#include "rar/PropertyMapSizingField.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/number_utils.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rar {

inline void run_c1_external_remeshing(Mesh& mesh, const RemeshConfig& cfg)
{
    if (cfg.c1_field_csv.empty()) {
        throw std::invalid_argument("c1-external requires --sizing-csv");
    }
    std::ifstream csv(cfg.c1_field_csv);
    if (!csv) {
        throw std::runtime_error("Cannot open C1 CSV: " + cfg.c1_field_csv);
    }
    std::string line;
    if (!std::getline(csv, line) ||
        line.find("vertex_id,x,y,z,target_length,") != 0) {
        throw std::runtime_error("Invalid C1 CSV header");
    }

    std::vector<Mesh::Vertex_index> original;
    original.reserve(num_vertices(mesh));
    for (const auto v : vertices(mesh)) original.push_back(v);
    const auto targets =
        mesh.add_property_map<Mesh::Vertex_index,double>(
            "v:c1_external_target_length",0.0).first;

    double scale = 1.0;
    for (const auto v : original) {
        const auto p = mesh.point(v);
        scale = (std::max)(scale,std::abs(CGAL::to_double(p.x())));
        scale = (std::max)(scale,std::abs(CGAL::to_double(p.y())));
        scale = (std::max)(scale,std::abs(CGAL::to_double(p.z())));
    }
    const double coord_tolerance = 1e-6 * scale;

    std::size_t row = 0;
    double sum = 0.0;
    while (std::getline(csv,line)) {
        if (line.empty()) continue;
        if (row >= original.size()) {
            throw std::runtime_error("C1 CSV contains too many vertices");
        }
        std::istringstream stream(line);
        std::string item;
        auto next = [&]() {
            if (!std::getline(stream,item,',')) {
                throw std::runtime_error(
                    "Malformed C1 CSV at vertex " + std::to_string(row));
            }
            return item;
        };
        const std::size_t id = static_cast<std::size_t>(std::stoull(next()));
        const double x = std::stod(next());
        const double y = std::stod(next());
        const double z = std::stod(next());
        const double h = std::stod(next());
        if (id != row || !std::isfinite(h) || h <= 0 ||
            !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            throw std::runtime_error("Invalid C1 CSV ID, coordinate, or length");
        }
        const auto p = mesh.point(original[row]);
        const double error = (std::max)({
            std::abs(CGAL::to_double(p.x())-x),
            std::abs(CGAL::to_double(p.y())-y),
            std::abs(CGAL::to_double(p.z())-z)
        });
        if (error > coord_tolerance) {
            throw std::runtime_error(
                "C1 CSV does not match input mesh at row " +
                std::to_string(row));
        }
        targets[original[row]] = h;
        sum += h;
        ++row;
    }
    if (!row || row != original.size()) {
        throw std::runtime_error("C1 CSV vertex count mismatch");
    }

    PropertyMapSizingField field(targets,sum/static_cast<double>(row));
    const FeatureConstraints features = make_feature_constraints(
        mesh,cfg.preserve_features,cfg.feature_angle_degrees);
    namespace PMP = CGAL::Polygon_mesh_processing;
    PMP::isotropic_remeshing(
        faces(mesh),field,mesh,
        CGAL::parameters::number_of_iterations(cfg.iterations)
            .number_of_relaxation_steps(cfg.relaxation_steps)
            .edge_is_constrained_map(features.edge_map)
            .vertex_is_constrained_map(features.vertex_map)
            .collapse_constraints(false)
            .relax_constraints(true)
            .do_project(cfg.do_project)
    );
}

} // namespace rar
