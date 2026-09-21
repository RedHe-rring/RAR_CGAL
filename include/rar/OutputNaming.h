#pragma once

#include "rar/RemeshConfig.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>

namespace rar {

inline std::string filename_number(const double value) {
    std::ostringstream oss;
    oss << std::setprecision(10) << std::defaultfloat << value;

    std::string s = oss.str();
    std::replace(s.begin(), s.end(), '.', 'p');

    // Keep scientific notation readable and filename-safe:
    // 1e-05 -> 1em05, 1e+05 -> 1ep05.
    std::string out;
    out.reserve(s.size());
    for (const char ch : s) {
        if (ch == '-') {
            out += 'm';
        } else if (ch == '+') {
            out += 'p';
        } else {
            out += ch;
        }
    }
    return out;
}

inline std::string make_auto_output_path(
    const RemeshConfig& cfg)
{
    const std::filesystem::path input(cfg.input_path);

    std::string extension = input.extension().string();
    if (extension.empty()) {
        extension = ".ply";
    }

    const std::string stem =
        input.stem().empty()
            ? std::string("remesh")
            : input.stem().string();

    std::ostringstream name;
    name
        << stem
        << "__field-" << field_type_name(cfg.field_type)
        << "__eps-" << filename_number(cfg.epsilon)
        << "__lmin-" << filename_number(cfg.min_edge_length)
        << "__lmax-" << filename_number(cfg.max_edge_length)
        << "__it-" << cfg.iterations
        << "__relax-" << cfg.relaxation_steps
        << "__proj-" << (cfg.do_project ? "on" : "off")
        << extension;

    return (input.parent_path() / name.str()).string();
}

} // namespace rar
