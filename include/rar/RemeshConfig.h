#pragma once

#include <string>

namespace rar {

enum class FieldType {
    CGALAdaptive,
    RAR
};

inline const char* field_type_name(const FieldType type) {
    switch (type) {
    case FieldType::CGALAdaptive:
        return "cgal-adaptive";
    case FieldType::RAR:
        return "rar";
    }
    return "unknown";
}

struct RemeshConfig {
    std::string input_path;
    std::string output_path;
    std::string export_field_prefix;
    bool output_path_auto = false;

    FieldType field_type = FieldType::CGALAdaptive;

    double epsilon = 1e-3;
    double min_edge_length = 1e-3;
    double max_edge_length = 5e-1;

    unsigned int iterations = 5;
    unsigned int relaxation_steps = 3;
    bool do_project = true;
};

} // namespace rar
