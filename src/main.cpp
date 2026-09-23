#include "rar/CGALAdaptiveRemesher.h"
#include "rar/CSFRemesher.h"
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
        << "  --field <name>          cgal-adaptive | rar | csf (default: cgal-adaptive)\n"
        << "  --epsilon <value>       RAR/CGAL approximation tolerance (default: 0.001)\n"
        << "  --min-edge <value>      RAR/CGAL minimum target edge length (default: 0.001)\n"
        << "  --max-edge <value>      RAR/CGAL maximum target edge length (default: 0.5)\n"
        << "  --mesh-scale <value>    CSF global target-length scale (default: 1.0)\n"
        << "  --iterations <n>        Remeshing iterations (default: 5)\n"
        << "  --relax-steps <n>       Relaxation steps per iteration (default: 3)\n"
        << "  --export-field <prefix> Override field-output stem (.ply/.csv added)\n"
        << "  --no-export-field       Disable default field export\n"
        << "  --no-project            Disable projection to the input surface\n\n"
        << "Field notes:\n"
        << "  cgal-adaptive  CGAL 6.1.x Adaptive_sizing_field.\n"
        << "  rar            RAR cotangent-curvature sizing field with CGAL local operators.\n"
        << "  csf            CSF code-oriented smoothed-curvature field with CGAL local operators.\n";
}

rar::FieldType parse_field_type(
    const std::string& value)
{
    if (value == "cgal-adaptive" ||
        value == "cgal") {
        return rar::FieldType::CGALAdaptive;
    }

    if (value == "rar") {
        return rar::FieldType::RAR;
    }

    if (value == "csf") {
        return rar::FieldType::CSF;
    }

    throw std::invalid_argument(
        "Unknown field '" + value +
        "'. Expected cgal-adaptive, rar, or csf.");
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
        << "  sizing min/mean/max    = "
        << s.sizing_min << " / "
        << s.sizing_mean << " / "
        << s.sizing_max << '\n';
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
                   rar::FieldType::RAR) {
            const rar::RARFieldStats field_stats =
                rar::run_rar_field_cgal_remeshing(
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

        print_mesh_stats(mesh, "Output");

        const double seconds =
            std::chrono::duration<double>(
                end - begin)
                .count();

        std::cout
            << "Remeshing runtime: "
            << seconds << " s\n";

        if (!CGAL::IO::write_polygon_mesh(
                cfg.output_path,
                mesh,
                CGAL::parameters::
                    stream_precision(17))) {
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
