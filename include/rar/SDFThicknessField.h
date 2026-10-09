#pragma once

#include "rar/RemeshConfig.h"
#include "rar/ChenSizingCorrection.h"
#include "rar/SDFThicknessSizing.h"
#include "rar/Types.h"

#include <cstddef>

namespace rar {

struct SDFThicknessStats {
    std::size_t face_count = 0;
    std::size_t vertex_count = 0;
    std::size_t thickness_limited_vertices = 0;
    std::size_t min_edge_clipped_vertices = 0;
    std::size_t max_edge_clipped_vertices = 0;
    double thickness_limited_area_fraction = 0.0;
    double face_thickness_min = 0.0;
    double face_thickness_max = 0.0;
    double target_min = 0.0;
    double target_max = 0.0;
};

struct SDFThicknessChenStats {
    SDFThicknessStats initial;
    ChenCorrectionStats correction;
    double corrected_min = 0.0;
    double corrected_max = 0.0;
    double corrected_mean = 0.0;
};

struct SDFThicknessField {
    Mesh::Property_map<Mesh::Face_index, double> face_thickness;
    Mesh::Property_map<Mesh::Face_index, double> face_area;
    Mesh::Property_map<Mesh::Vertex_index, double> thickness;
    Mesh::Property_map<Mesh::Vertex_index, double> curvature;
    Mesh::Property_map<Mesh::Vertex_index, double> vertex_area;
    Mesh::Property_map<Mesh::Vertex_index, double> target;
    Mesh::Property_map<Mesh::Vertex_index, SDFThicknessSizingSample> samples;
    SDFThicknessStats stats;
};

void validate_sdf_thickness_config(const RemeshConfig& cfg);
// Checks the input before any curvature estimation or ray intersection.
void validate_sdf_thickness_mesh(const Mesh& mesh);
// Raw, unnormalized CGAL SDF, then area-weighted incident-face averaging.
// Missing SDF, invalid curvature, or invalid geometry is an explicit error.
SDFThicknessField make_sdf_thickness_field(Mesh& mesh, const RemeshConfig& cfg);
void write_sdf_thickness_field(const Mesh& mesh, const RemeshConfig& cfg,
                             const SDFThicknessField& field,
                             const Mesh::Property_map<Mesh::Vertex_index, double>* corrected = nullptr,
                             const ChenCorrectionStats* correction = nullptr);
SDFThicknessStats run_sdf_thickness_remeshing(Mesh& mesh, const RemeshConfig& cfg);
SDFThicknessChenStats run_sdf_thickness_chen_remeshing(Mesh& mesh, const RemeshConfig& cfg);

} // namespace rar
