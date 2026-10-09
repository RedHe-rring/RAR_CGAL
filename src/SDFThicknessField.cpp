#include "rar/SDFThicknessField.h"

#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/mesh_segmentation.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace rar {
namespace PMP = CGAL::Polygon_mesh_processing;

void validate_sdf_thickness_config(const RemeshConfig& cfg)
{
    // Validate common parameters using the same contract as the formula.
    sdf_thickness_target_length(0.0, 1.0, cfg.epsilon,
                                cfg.min_edge_length, cfg.max_edge_length);
    if (cfg.sdf_rays == 0)
        throw std::invalid_argument("sdf-rays must be > 0");
    if (!std::isfinite(cfg.sdf_cone_angle_degrees) ||
        !(cfg.sdf_cone_angle_degrees > 0.0) || !(cfg.sdf_cone_angle_degrees < 180.0))
        throw std::invalid_argument("sdf-cone-angle must be finite and in (0, 180) degrees");
}

bool validate_sdf_thickness_mesh(const Mesh& mesh)
{
    if (mesh.number_of_faces() == 0 || !CGAL::is_valid_polygon_mesh(mesh) ||
        !CGAL::is_triangle_mesh(mesh))
        throw std::invalid_argument("SDF thickness requires a nonempty valid triangle mesh");
    if (!CGAL::is_closed(mesh))
        throw std::invalid_argument("SDF thickness requires a closed mesh (no boundary edges)");
    for (const auto v : mesh.vertices()) {
        const auto& p = mesh.point(v);
        if (!std::isfinite(CGAL::to_double(p.x())) ||
            !std::isfinite(CGAL::to_double(p.y())) ||
            !std::isfinite(CGAL::to_double(p.z())) || mesh.is_isolated(v))
            throw std::invalid_argument("SDF thickness requires finite, non-isolated vertices");
    }
    std::vector<Mesh::Halfedge_index> nonmanifold;
    PMP::non_manifold_vertices(mesh, std::back_inserter(nonmanifold));
    if (!nonmanifold.empty())
        throw std::invalid_argument("SDF thickness requires a 2-manifold mesh");
    for (const auto f : mesh.faces()) {
        const double area = CGAL::to_double(PMP::face_area(f, mesh));
        if (!std::isfinite(area) || !(area > 0.0))
            throw std::invalid_argument("SDF thickness rejects degenerate or non-finite face " +
                                        std::to_string(f.idx()));
    }
    // The volume/orientation routines assume an intersection-free mesh.
    // SDF can still be attempted, but its ray lengths may be unreliable.
    if (PMP::does_self_intersect(mesh))
        return true;
    if (!PMP::does_bound_a_volume(mesh) || !(CGAL::to_double(PMP::volume(mesh)) > 0.0))
        throw std::invalid_argument("SDF thickness requires consistently oriented boundaries of a bounded solid");
    return false;
}

SDFThicknessField make_sdf_thickness_field(Mesh& mesh, const RemeshConfig& cfg)
{
    validate_sdf_thickness_config(cfg);
    const bool input_self_intersections = validate_sdf_thickness_mesh(mesh);
    if (input_self_intersections)
        std::cerr << "Warning: SDF thickness input has self-intersections; "
                     "ray-based thickness may be unreliable. "
                     "Solid orientation/volume validation was skipped.\n";
    auto face_t = mesh.add_property_map<Mesh::Face_index, double>("f:sdf_raw_thickness", -1.0).first;
    auto face_area = mesh.add_property_map<Mesh::Face_index, double>("f:sdf_area", 0.0).first;
    // CGAL already removes outlier ray lengths in this raw computation.
    // postprocess=true would also smooth and normalize to [0,1], losing units.
    CGAL::sdf_values(mesh, face_t, cfg.sdf_cone_angle_degrees * CGAL_PI / 180.0,
                     cfg.sdf_rays, false);
    std::size_t invalid_faces = 0;
    for (const auto f : mesh.faces())
        if (!std::isfinite(face_t[f]) || !(face_t[f] > 0.0)) ++invalid_faces;
    if (invalid_faces != 0)
        throw std::runtime_error("CGAL SDF returned missing/non-positive/non-finite thickness on " +
            std::to_string(invalid_faces) + " faces; no thickness fallback is enabled");

    auto thickness = mesh.add_property_map<Mesh::Vertex_index, double>("v:sdf_thickness", 0.0).first;
    auto curvature = mesh.add_property_map<Mesh::Vertex_index, double>("v:sdf_curvature", 0.0).first;
    auto vertex_area = mesh.add_property_map<Mesh::Vertex_index, double>("v:sdf_vertex_area", 0.0).first;
    auto target = mesh.add_property_map<Mesh::Vertex_index, double>("v:sdf_h0", cfg.max_edge_length).first;
    auto samples = mesh.add_property_map<Mesh::Vertex_index, SDFThicknessSizingSample>(
        "v:sdf_sizing_sample", SDFThicknessSizingSample{}).first;
    // Reinitialize existing properties too: repeated construction must not accumulate.
    for (const auto v : mesh.vertices()) { thickness[v] = 0.0; vertex_area[v] = 0.0; }
    SDFThicknessStats stats;
    stats.input_self_intersections = input_self_intersections;
    stats.face_count = mesh.number_of_faces();
    stats.vertex_count = mesh.number_of_vertices();
    stats.face_thickness_min = stats.target_min = std::numeric_limits<double>::infinity();
    for (const auto f : mesh.faces()) {
        face_area[f] = CGAL::to_double(PMP::face_area(f, mesh));
        stats.face_thickness_min = (std::min)(stats.face_thickness_min, face_t[f]);
        stats.face_thickness_max = (std::max)(stats.face_thickness_max, face_t[f]);
        for (const auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            vertex_area[v] += face_area[f] / 3.0;
            thickness[v] += (face_area[f] / 3.0) * face_t[f];
        }
    }
    using PrincipalCurvatures = PMP::Principal_curvatures_and_directions<Kernel>;
    auto principal = mesh.add_property_map<Mesh::Vertex_index, PrincipalCurvatures>(
        "v:sdf_principal_curvatures", PrincipalCurvatures{}).first;
    PMP::interpolated_corrected_curvatures(mesh,
        CGAL::parameters::vertex_principal_curvatures_and_directions_map(principal));
    double total_area = 0.0, limited_area = 0.0;
    for (const auto v : mesh.vertices()) {
        if (!std::isfinite(vertex_area[v]) || !(vertex_area[v] > 0.0))
            throw std::runtime_error("Invalid incident area for SDF vertex " + std::to_string(v.idx()));
        thickness[v] /= vertex_area[v];
        const double k1 = CGAL::to_double(principal[v].min_curvature);
        const double k2 = CGAL::to_double(principal[v].max_curvature);
        if (!std::isfinite(k1) || !std::isfinite(k2))
            throw std::runtime_error("Invalid curvature for SDF vertex " + std::to_string(v.idx()));
        curvature[v] = (std::max)(std::abs(k1), std::abs(k2));
        samples[v] = sdf_thickness_target_length(curvature[v], thickness[v], cfg.epsilon,
                                                 cfg.min_edge_length, cfg.max_edge_length);
        target[v] = samples[v].target_length;
        stats.target_min = (std::min)(stats.target_min, target[v]);
        stats.target_max = (std::max)(stats.target_max, target[v]);
        total_area += vertex_area[v];
        if (samples[v].thickness_limited) {
            ++stats.thickness_limited_vertices;
            limited_area += vertex_area[v];
        }
        stats.min_edge_clipped_vertices += samples[v].min_edge_clipped;
        stats.max_edge_clipped_vertices += samples[v].max_edge_clipped;
    }
    stats.thickness_limited_area_fraction = limited_area / total_area;
    return {face_t, face_area, thickness, curvature, vertex_area, target, samples, stats};
}

} // namespace rar
