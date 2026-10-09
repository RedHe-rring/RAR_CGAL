#pragma once

#include <string>

namespace rar {

enum class RemeshBackend { CGAL, Custom };
enum class DiagnosticLevel { Stage, Operation };

inline const char* remesh_backend_name(const RemeshBackend backend)
{
    return backend == RemeshBackend::Custom ? "custom" : "cgal";
}

enum class FieldType {
    CGALAdaptive,
    CGALAdaptiveRadius,
    CGALAdaptiveSDFThickness,
    CGALAdaptiveSDFThicknessChen,
    CGALAdaptiveChen,
    CGALAdaptiveRadiusChen,
    RAR,
    RARChen,
    CSF
};

inline const char* field_type_name(
    const FieldType type)
{
    switch (type) {
    case FieldType::CGALAdaptive:
        return "cgal-adaptive";
    case FieldType::CGALAdaptiveRadius:
        return "cgal-adaptive-radius";
    case FieldType::CGALAdaptiveSDFThickness:
        return "cgal-adaptive-sdf-thickness";
    case FieldType::CGALAdaptiveSDFThicknessChen:
        return "cgal-adaptive-sdf-thickness-chen";
    case FieldType::CGALAdaptiveChen:
        return "cgal-adaptive-chen";
    case FieldType::CGALAdaptiveRadiusChen:
        return "cgal-adaptive-radius-chen";
    case FieldType::RAR:
        return "rar";
    case FieldType::RARChen:
        return "rar-chen";
    case FieldType::CSF:
        return "csf";
    }
    return "unknown";
}

struct RemeshConfig {
    std::string input_path;
    std::string output_path;
    std::string export_field_prefix;
    bool output_path_auto = false;
    bool export_field = true;

    FieldType field_type =
        FieldType::CGALAdaptive;

    // RAR / RAR+Chen / CGAL-Adaptive / CGAL-Adaptive+Chen parameters.
    double epsilon = 1e-3;
    bool epsilon_is_explicit = false;
    double min_edge_length = 1e-3;
    double max_edge_length = 5e-1;

    // Raw CGAL Shape Diameter Function parameters (no normalization).
    unsigned int sdf_rays = 25;
    double sdf_cone_angle_degrees = 120.0;

    // Chen sizing-field gradation parameter.
    double chen_beta = 1.2;

    // CSF field parameter.
    double csf_mesh_scale = 1.0;

    // Common remeshing backend parameters.
    RemeshBackend remesh_backend = RemeshBackend::CGAL;
    unsigned int iterations = 5;
    unsigned int relaxation_steps = 3;
    bool do_project = true;
    bool remesh_diagnostics = false;
    DiagnosticLevel diagnostic_level = DiagnosticLevel::Stage;

    // Sharp-feature detection and preservation. CGAL still treats mesh
    // boundary edges as constrained independently of this setting.
    bool preserve_features = false;
    double feature_angle_degrees = 50.0;
};

} // namespace rar
