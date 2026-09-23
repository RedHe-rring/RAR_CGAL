#pragma once

#include "rar/RemeshConfig.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>

namespace rar {

inline std::string filename_number(
    const double value)
{
    std::ostringstream oss;
    oss
        << std::setprecision(10)
        << std::defaultfloat
        << value;

    std::string s = oss.str();
    std::replace(
        s.begin(), s.end(), '.', 'p');

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

inline std::string make_parameterized_stem(
    const RemeshConfig& cfg)
{
    const std::filesystem::path input(
        cfg.input_path);

    const std::string stem =
        input.stem().empty()
            ? std::string("remesh")
            : input.stem().string();

    std::ostringstream name;
    name << stem;

    switch (cfg.field_type) {
    case FieldType::CGALAdaptive:
    case FieldType::RAR:
        name
            << "__eps-"
            << filename_number(cfg.epsilon)
            << "__lmin-"
            << filename_number(
                   cfg.min_edge_length)
            << "__lmax-"
            << filename_number(
                   cfg.max_edge_length);
        break;

    case FieldType::CSF:
        name
            << "__scale-"
            << filename_number(
                   cfg.csf_mesh_scale);
        break;
    }

    name
        << "__it-" << cfg.iterations
        << "__relax-"
        << cfg.relaxation_steps
        << "__proj-"
        << (cfg.do_project ? "on" : "off")
        << "__field-"
        << field_type_name(cfg.field_type);

    return name.str();
}

inline std::string make_auto_output_path(
    const RemeshConfig& cfg)
{
    const std::filesystem::path input(
        cfg.input_path);

    std::string extension =
        input.extension().string();

    if (extension.empty()) {
        extension = ".ply";
    }

    return (
        input.parent_path()
        /
        (make_parameterized_stem(cfg) +
         extension)
    ).string();
}

inline std::string make_auto_field_prefix(
    const RemeshConfig& cfg)
{
    const std::filesystem::path output(
        cfg.output_path);

    const std::string output_stem =
        output.stem().empty()
            ? std::string("remesh")
            : output.stem().string();

    // Keep diagnostic artifacts inside a folder
    // named after the output mesh, e.g.:
    // result.obj
    // result/field.ply
    // result/field.csv
    return (
        output.parent_path()
        / output_stem
        / "field"
    ).string();
}

} // namespace rar
