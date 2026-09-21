#pragma once

#include <string>

namespace rar {

struct RemeshConfig {
    std::string input_path;
    std::string output_path;

    double epsilon = 1e-3;
    double min_edge_length = 1e-3;
    double max_edge_length = 5e-1;

    unsigned int iterations = 5;
    unsigned int relaxation_steps = 3;
    bool do_project = true;
};

} // namespace rar
