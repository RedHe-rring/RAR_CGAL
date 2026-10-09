#include "rar/Wang2019Remesher.h"
#if __has_include(<CGAL/IO/polygon_mesh_io.h>)
#  include <CGAL/IO/polygon_mesh_io.h>
#else
#  include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#endif
#include <CGAL/boost/graph/helpers.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void usage() {
    std::cout <<
        "wang2019_remesh input.ply output.ply [options]\n"
        "  --min-angle D      small-angle threshold in degrees (default 30)\n"
        "  --max-angle D      large-angle threshold in degrees (default 90)\n"
        "  --rounds N         maximal alternating rounds (default 10)\n"
        "  --budget F         attempted operations / bad triangle count (default 0.02)\n"
        "  --k N              override automatic number of operations per stage\n"
        "  --allow-drift      allow vertex count changes instead of rollback\n"
        "  --smooth N         tangential smoothing sweeps (default 3)\n"
        "  --feature-angle D  sharp edge dihedral threshold (default 50)\n"
        "  --no-features      preserve boundary only; allow changes near creases\n"
        "  --no-project       disable closest-point projection to input surface\n"
        "  --quiet            disable per-round diagnostics\n"
        "\n"
        "Experimental Wang et al. (TVCG 2019)-inspired postprocessor.\n"
        "Not an exact reproduction of the author's complete pipeline.\n";
}
}
int main(int argc, char** argv) {
    try {
        if (argc < 3 || std::string(argv[1]) == "--help") {
            usage();
            return argc < 3 ? 2 : 0;
        }
        const std::filesystem::path input(argv[1]), output(argv[2]);
        rar::wang2019::Options opt;
        for (int i = 3; i < argc; ++i) {
            const std::string arg(argv[i]);
            auto value = [&]() -> std::string {
                if (++i >= argc) throw std::invalid_argument("Missing value for " + arg);
                return argv[i];
            };
            if (arg == "--min-angle") opt.min_angle = std::stod(value());
            else if (arg == "--max-angle") opt.max_angle = std::stod(value());
            else if (arg == "--rounds") opt.rounds = static_cast<unsigned>(std::stoul(value()));
            else if (arg == "--budget") opt.budget_fraction = std::stod(value());
            else if (arg == "--k") opt.k = static_cast<unsigned>(std::stoul(value()));
            else if (arg == "--allow-drift") opt.strict_vertex_count = false;
            else if (arg == "--smooth") opt.smoothing_steps = static_cast<unsigned>(std::stoul(value()));
            else if (arg == "--feature-angle") opt.feature_angle = std::stod(value());
            else if (arg == "--no-project") opt.project = false;
            else if (arg == "--no-features") opt.protect_features = false;
            else if (arg == "--quiet") opt.verbose = false;
            else throw std::invalid_argument("Unknown argument: " + arg);
        }
        rar::Mesh mesh;
        if (!CGAL::IO::read_polygon_mesh(input.string(), mesh) ||
            !CGAL::is_triangle_mesh(mesh))
            throw std::runtime_error("Could not read a triangular surface mesh: " + input.string());
        const rar::wang2019::Statistics s = rar::wang2019::run(mesh, opt);
        if (!output.parent_path().empty())
            std::filesystem::create_directories(output.parent_path());
        if (!CGAL::IO::write_polygon_mesh(output.string(), mesh))
            throw std::runtime_error("Cannot write output mesh: " + output.string());
        std::cout << "Wang2019 prototype: " << input << " -> " << output
                  << "\n  V: " << s.vertices_before << " -> " << s.vertices_after
                  << "  F: " << s.faces_before << " -> " << s.faces_after
                  << "\n  rolled-back rounds: " << s.rolled_back_rounds
                  << "\n  inserts: " << s.insertions
                  << " collapses: " << s.collapses << " flips: " << s.flips
                  << "\n  min/max angle: " << s.angles.min << " / " << s.angles.max
                  << "\n  below/above bound: " << s.min_angle_violations
                  << " / " << s.max_angle_violations << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "wang2019_remesh: " << e.what() << '\n';
        return 1;
    }
}
