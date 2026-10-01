#pragma once

#include <string>

namespace rar {

enum class FieldType {
    CGALAdaptive,
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
    double min_edge_length = 1e-3;
    double max_edge_length = 5e-1;

    // Chen sizing-field gradation parameter.
    double chen_beta = 1.2;

    // CSF field parameter.
    double csf_mesh_scale = 1.0;

    // Common CGAL remeshing backend parameters.
    unsigned int iterations = 5;
    unsigned int relaxation_steps = 3;
    bool do_project = true;
};

} // namespace rar
