#include "rar/CGALAdaptiveRemesher.h"
#include "rar/CSFRemesher.h"
#include "rar/MeshOutput.h"
#include "rar/OutputNaming.h"
#include "rar/RARRemesher.h"
#include "rar/RemeshConfig.h"
#include "rar/Types.h"

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace {

void print_usage(const char* exe) {
    std::cerr
        << "Usage:\n  " << exe
        << " <input_mesh> [output_mesh] [options]\n\n"
        << "If output_mesh is omitted, the output filename is generated from the input\n"
        << "name and the parameters that actually affect the selected field.\n\n"
        << "Options:\n"
        << "  --field <name>          cgal-adaptive | cgal-adaptive-chen | cgal-adaptive-radius-chen | rar | rar-chen | csf\n"
        << "  --epsilon <value>       RAR/CGAL approximation tolerance (default: 0.001)\n"
        << "  --min-edge <value>      RAR/CGAL minimum target edge length (default: 0.001)\n"
        << "  --max-edge <value>      RAR/CGAL maximum target edge length (default: 0.5)\n"
        << "  --beta <value>          Chen progressive factor for *-chen modes (default: 1.2)\n"
        << "  --mesh-scale <value>    CSF global target-length scale (default: 1.0)\n"
        << "  --iterations <n>        Remeshing iterations (default: 5)\n"
        << "  --relax-steps <n>       Relaxation steps per iteration (default: 3)\n"
        << "  --export-field <prefix> Override field-output stem (.ply/.csv added)\n"
        << "  --no-export-field       Disable default field export\n"
        << "  --no-project            Disable projection to the input surface\n\n"
        << "Field notes:\n"
        << "  cgal-adaptive       CGAL 6.1.x Adaptive_sizing_field.\n"
        << "  cgal-adaptive-chen  CGAL adaptive field followed by Chen correction.\n"
        << "  cgal-adaptive-radius-chen  Curvature-radius-guarded CGAL field followed by Chen correction.\n"
        << "  rar                 RAR cotangent-curvature sizing field with CGAL local operators.\n"
        << "  rar-chen            RAR field followed by Chen gradient-constrained correction.\n"
        << "  csf                 CSF code-oriented smoothed-curvature field with CGAL local operators.\n";
}

rar::FieldType parse_field_type(
    const std::string& value)
{
    if (value == "cgal-adaptive" ||
        value == "cgal") {
        return rar::FieldType::CGALAdaptive;
    }

    if (value == "cgal-adaptive-chen" ||
        value == "cgal-chen") {
        return rar::FieldType::CGALAdaptiveChen;
    }

    if (value == "cgal-adaptive-radius-chen" ||
        value == "cgal-radius-chen") {
        return rar::FieldType::CGALAdaptiveRadiusChen;
    }

    if (value == "rar") {
        return rar::FieldType::RAR;
    }

    if (value == "rar-chen") {
        return rar::FieldType::RARChen;
    }

    if (value == "csf") {
        return rar::FieldType::CSF;
    }

    throw std::invalid_argument(
        "Unknown field '" + value +
        "'. Expected cgal-adaptive, cgal-adaptive-chen, "
        "cgal-adaptive-radius-chen, rar, rar-chen, or csf.");
}

rar::RemeshConfig parse_args(
    int argc,
    char** argv)
{
    if (argc == 2) {
        const std::string first = argv[1];
        if (first == "--help" ||
            first == "-h") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }
    }

    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(EXIT_FAILURE);
    }

    rar::RemeshConfig cfg;
    cfg.input_path = argv[1];

    int option_start = 2;

    if (argc > 2) {
        const std::string second = argv[2];
        if (!second.empty() &&
            second[0] != '-') {
            cfg.output_path = second;
            option_start = 3;
        }
    }

    for (int i = option_start;
         i < argc;
         ++i) {
        const std::string arg = argv[i];

        auto require_value =
            [&](const std::string& name)
                -> std::string {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "Missing value for " + name);
                }
                return argv[++i];
            };

        if (arg == "--field") {
            cfg.field_type =
                parse_field_type(
                    require_value(arg));
        } else if (arg == "--epsilon") {
            cfg.epsilon =
                std::stod(require_value(arg));
        } else if (arg == "--min-edge") {
            cfg.min_edge_length =
                std::stod(require_value(arg));
        } else if (arg == "--max-edge") {
            cfg.max_edge_length =
                std::stod(require_value(arg));
        } else if (arg == "--beta") {
            cfg.chen_beta =
                std::stod(require_value(arg));
        } else if (arg == "--mesh-scale") {
            cfg.csf_mesh_scale =
                std::stod(require_value(arg));
        } else if (arg == "--iterations") {
            cfg.iterations =
                static_cast<unsigned int>(
                    std::stoul(
                        require_value(arg)));
        } else if (arg == "--relax-steps") {
            cfg.relaxation_steps =
                static_cast<unsigned int>(
                    std::stoul(
                        require_value(arg)));
        } else if (arg == "--export-field") {
            cfg.export_field_prefix =
                require_value(arg);
            cfg.export_field = true;
        } else if (arg ==
                   "--no-export-field") {
            cfg.export_field = false;
        } else if (arg == "--no-project") {
            cfg.do_project = false;
        } else if (arg == "--help" ||
                   arg == "-h") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::invalid_argument(
                "Unknown option: " + arg);
        }
    }

    if (cfg.output_path.empty()) {
        cfg.output_path =
            rar::make_auto_output_path(cfg);
        cfg.output_path_auto = true;
    }

    if (cfg.export_field &&
        cfg.export_field_prefix.empty()) {
        cfg.export_field_prefix =
            rar::make_auto_field_prefix(cfg);
    }

    return cfg;
}

void validate_config(
    const rar::RemeshConfig& cfg)
{
    if (cfg.field_type ==
        rar::FieldType::CSF) {
        if (!(cfg.csf_mesh_scale > 0.0)) {
            throw std::invalid_argument(
                "mesh_scale must be > 0");
        }
        return;
    }

    if (!(cfg.epsilon > 0.0)) {
        throw std::invalid_argument(
            "epsilon must be > 0");
    }

    if (!(cfg.min_edge_length > 0.0) ||
        !(cfg.max_edge_length > 0.0)) {
        throw std::invalid_argument(
            "edge length bounds must be > 0");
    }

    if (cfg.min_edge_length >
        cfg.max_edge_length) {
        throw std::invalid_argument(
            "min_edge_length must be <= "
            "max_edge_length");
    }

    if ((cfg.field_type ==
             rar::FieldType::RARChen ||
         cfg.field_type ==
             rar::FieldType::CGALAdaptiveChen ||
         cfg.field_type ==
             rar::FieldType::CGALAdaptiveRadiusChen) &&
        !(cfg.chen_beta > 1.0)) {
        throw std::invalid_argument(
            "beta must be > 1 for Chen-corrected fields");
    }
}

void print_mesh_stats(
    const rar::Mesh& mesh,
    const char* label)
{
    std::cout
        << label
        << ": vertices=" << num_vertices(mesh)
        << ", edges=" << num_edges(mesh)
        << ", faces=" << num_faces(mesh)
        << '\n';
}

void print_rar_field_stats(
    const rar::RARFieldStats& s)
{
    std::cout
        << "RAR initial field:\n"
        << "  curvature min/mean/max = "
        << s.curvature_min << " / "
        << s.curvature_mean << " / "
        << s.curvature_max << '\n'
        << "  raw sizing min/mean/max = "
        << s.sizing_min << " / "
        << s.sizing_mean << " / "
        << s.sizing_max << '\n';

    if (s.chen_applied) {
        std::cout
            << "Chen correction:\n"
            << "  beta                    = "
            << s.chen_beta << '\n'
            << "  corrected min/mean/max  = "
            << s.corrected_sizing_min << " / "
            << s.corrected_sizing_mean << " / "
            << s.corrected_sizing_max << '\n'
            << "  max |grad h| before/after = "
            << s.chen_max_gradient_before << " / "
            << s.chen_max_gradient_after << '\n'
            << "  changed vertices        = "
            << s.chen_changed_vertex_count << " / "
            << s.vertex_count << '\n'
            << "  projection objective    = "
            << s.chen_objective << '\n';
    }
}

void print_cgal_adaptive_chen_stats(
    const rar::CGALAdaptiveChenStats& s)
{
    std::cout
        << "CGAL-Adaptive initial field:\n"
        << "  raw sizing min/mean/max = "
        << s.raw_sizing_min << " / "
        << s.raw_sizing_mean << " / "
        << s.raw_sizing_max << '\n'
        << "Chen correction:\n"
        << "  beta                    = "
        << s.chen_beta << '\n'
        << "  corrected min/mean/max  = "
        << s.corrected_sizing_min << " / "
        << s.corrected_sizing_mean << " / "
        << s.corrected_sizing_max << '\n'
        << "  max |grad h| before/after = "
        << s.chen_max_gradient_before << " / "
        << s.chen_max_gradient_after << '\n'
        << "  changed vertices        = "
        << s.chen_changed_vertex_count << " / "
        << s.vertex_count << '\n'
        << "  projection objective    = "
        << s.chen_objective << '\n';
}

void print_csf_field_stats(
    const rar::CSFFieldStats& s)
{
    std::cout
        << "CSF initial field:\n"
        << "  base edge length       = "
        << s.base_edge_length << '\n'
        << "  raw curvature min/mean/max = "
        << s.raw_curvature_min << " / "
        << s.raw_curvature_mean << " / "
        << s.raw_curvature_max << '\n'
        << "  smooth curvature min/mean/max = "
        << s.smoothed_curvature_min << " / "
        << s.smoothed_curvature_mean << " / "
        << s.smoothed_curvature_max << '\n'
        << "  sizing min/mean/max    = "
        << s.sizing_min << " / "
        << s.sizing_mean << " / "
        << s.sizing_max << '\n';
}

void print_config(
    const rar::RemeshConfig& cfg)
{
    std::cout
        << std::setprecision(10)
        << "field="
        << rar::field_type_name(
               cfg.field_type);

    if (cfg.field_type ==
        rar::FieldType::CSF) {
        std::cout
            << ", mesh_scale="
            << cfg.csf_mesh_scale;
    } else {
        std::cout
            << ", epsilon=" << cfg.epsilon
            << ", min_edge="
            << cfg.min_edge_length
            << ", max_edge="
            << cfg.max_edge_length;

        if (cfg.field_type ==
                rar::FieldType::RARChen ||
            cfg.field_type ==
                rar::FieldType::CGALAdaptiveChen ||
            cfg.field_type ==
                rar::FieldType::CGALAdaptiveRadiusChen) {
            std::cout
                << ", beta="
                << cfg.chen_beta;
        }
    }

    std::cout
        << ", iterations=" << cfg.iterations
        << ", relax_steps="
        << cfg.relaxation_steps
        << ", project="
        << (cfg.do_project
                ? "true"
                : "false")
        << ", export_field="
        << (cfg.export_field
                ? "true"
                : "false");

    if (cfg.export_field) {
        std::cout
            << ", field_prefix="
            << cfg.export_field_prefix;
    }

    std::cout << '\n';
}

} // namespace

int main(
    int argc,
    char** argv)
{
    try {
        const rar::RemeshConfig cfg =
            parse_args(argc, argv);
        validate_config(cfg);

        rar::Mesh mesh;

        if (!PMP::IO::read_polygon_mesh(
                cfg.input_path,
                mesh) ||
            mesh.is_empty()) {
            std::cerr
                << "Failed to read a non-empty "
                   "polygon mesh: "
                << cfg.input_path << '\n';
            return EXIT_FAILURE;
        }

        if (!CGAL::is_triangle_mesh(mesh)) {
            std::cout
                << "Input is not triangular; "
                   "triangulating faces...\n";

            if (!PMP::triangulate_faces(mesh)) {
                std::cerr
                    << "Triangulation failed.\n";
                return EXIT_FAILURE;
            }
        }

        if (!CGAL::is_triangle_mesh(mesh)) {
            std::cerr
                << "Input could not be converted "
                   "to a triangle mesh.\n";
            return EXIT_FAILURE;
        }

        print_mesh_stats(mesh, "Input");
        print_config(cfg);

        std::cout
            << (cfg.output_path_auto
                    ? "Auto output: "
                    : "Output: ")
            << cfg.output_path << '\n';

        if (cfg.export_field) {
            std::cout
                << "Field outputs: "
                << cfg.export_field_prefix
                << ".ply, "
                << cfg.export_field_prefix
                << ".csv\n";
        }

        const auto begin =
            std::chrono::steady_clock::now();

        if (cfg.field_type ==
            rar::FieldType::CGALAdaptive) {
            rar::run_cgal_adaptive_remeshing(
                mesh, cfg);
        } else if (cfg.field_type ==
                   rar::FieldType::CGALAdaptiveChen) {
            const rar::CGALAdaptiveChenStats field_stats =
                rar::run_cgal_adaptive_chen_remeshing(
                    mesh, cfg);
            print_cgal_adaptive_chen_stats(
                field_stats);
        } else if (cfg.field_type ==
                   rar::FieldType::CGALAdaptiveRadiusChen) {
            const rar::CGALAdaptiveChenStats field_stats =
                rar::run_cgal_adaptive_radius_chen_remeshing(
                    mesh, cfg);
            print_cgal_adaptive_chen_stats(
                field_stats);
        } else if (cfg.field_type ==
                   rar::FieldType::RAR) {
            const rar::RARFieldStats field_stats =
                rar::run_rar_field_cgal_remeshing(
                    mesh, cfg);
            print_rar_field_stats(field_stats);
        } else if (cfg.field_type ==
                   rar::FieldType::RARChen) {
            const rar::RARFieldStats field_stats =
                rar::run_rar_chen_field_cgal_remeshing(
                    mesh, cfg);
            print_rar_field_stats(field_stats);
        } else {
            const rar::CSFFieldStats field_stats =
                rar::run_csf_field_cgal_remeshing(
                    mesh, cfg);
            print_csf_field_stats(field_stats);
        }

        const auto end =
            std::chrono::steady_clock::now();

        // CGAL edge collapses can leave Surface_mesh descriptors
        // sparse until garbage is collected. Some mesh writers
        // (notably PLY readers downstream such as MeshLab) expect
        // face indices to refer to a compact 0..N-1 vertex array.
        // Compact the mesh before exporting the final result.
        if (mesh.has_garbage()) {
            mesh.collect_garbage();
        }

        print_mesh_stats(mesh, "Output");

        const double seconds =
            std::chrono::duration<double>(
                end - begin)
                .count();

        std::cout
            << "Remeshing runtime: "
            << seconds << " s\n";

        bool write_ok = false;

        if (rar::has_ply_extension(
                cfg.output_path)) {
            // Write PLY through an explicit descriptor->dense-ID map.
            // This avoids relying on Surface_mesh internal descriptor
            // values after split/collapse operations.
            write_ok =
                rar::write_ascii_ply_compact(
                    cfg.output_path,
                    mesh,
                    17);
        } else {
            write_ok =
                CGAL::IO::write_polygon_mesh(
                    cfg.output_path,
                    mesh,
                    CGAL::parameters::
                        stream_precision(17));
        }

        if (!write_ok) {
            std::cerr
                << "Failed to write output mesh: "
                << cfg.output_path << '\n';
            return EXIT_FAILURE;
        }

        std::cout
            << "Wrote: "
            << cfg.output_path << '\n';

        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr
            << "Error: "
            << e.what()
            << '\n';
        return EXIT_FAILURE;
    }
}
