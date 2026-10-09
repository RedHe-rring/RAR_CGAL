#include "rar/RemeshDiagnostics.h"
#include "rar/OutputNaming.h"
#include "rar/BuildInfo.h"

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/version.h>
#include <boost/uuid/detail/sha1.hpp>

#include <array>
#include <ctime>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace rar {
namespace {
std::string json_string(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    for (const unsigned char ch : value) {
        switch (ch) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (ch < 0x20) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(ch);
            else out << static_cast<char>(ch);
        }
    }
    out << '"';
    return out.str();
}

std::string padded(std::size_t n, int width)
{
    std::ostringstream out;
    out << std::setw(width) << std::setfill('0') << n;
    return out.str();
}

std::string input_sha1(const std::string& path)
{
    if (path.empty()) return "null"; // An in-memory library call has no source file.
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read diagnostic input checksum: " + path);
    boost::uuids::detail::sha1 hash;
    std::array<char, 65536> buffer;
    while (input) {
        input.read(buffer.data(), buffer.size());
        hash.process_bytes(buffer.data(), static_cast<std::size_t>(input.gcount()));
    }
    if (!input.eof()) throw std::runtime_error("Failed to checksum diagnostic input: " + path);
    boost::uuids::detail::sha1::digest_type digest;
    hash.get_digest(digest);
    std::ostringstream out;
    for (const auto word : digest)
        out << std::hex << std::setw(sizeof(word) * 2) << std::setfill('0') << static_cast<unsigned long>(word);
    return json_string(out.str());
}

std::ofstream open_output(const std::filesystem::path& path)
{
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(path);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    return out;
}

const char* stage_file(const std::string& stage)
{
    if (stage == "split") return "01_split.off";
    if (stage == "collapse") return "02_collapse.off";
    if (stage == "flip") return "03_flip.off";
    if (stage == "relax") return "04_relax.off";
    if (stage == "project") return "05_project.off";
    throw std::runtime_error("Unknown diagnostic stage: " + stage);
}
} // namespace

RemeshDiagnostics::RemeshDiagnostics(const RemeshConfig& cfg)
{
    const auto root = std::filesystem::absolute(remesh_diagnostics_root(cfg));
    std::filesystem::create_directories(root);
    // create_directory is the exclusive claim; existing runs are never overwritten.
    for (std::size_t run = 1; ; ++run) {
        directory_ = root / ("run_" + padded(run, 4));
        if (std::filesystem::create_directory(directory_)) break;
    }
    events_ = open_output(directory_ / "events.jsonl");
    write_status("incomplete");
    const auto checksum = input_sha1(cfg.input_path);
    auto manifest = open_output(directory_ / "manifest.json");
    const auto now = std::time(nullptr);
    std::ostringstream utc;
    utc << std::put_time(std::gmtime(&now), "%Y-%m-%dT%H:%M:%SZ");
    manifest << "{\n"
        << "  \"schema_version\": 1,\n"
        << "  \"started_utc\": " << json_string(utc.str()) << ",\n"
        << "  \"source_commit\": " << json_string(build_info::commit) << ",\n"
        << "  \"source_worktree_dirty_at_build\": " << build_info::dirty_json << ",\n"
        << "  \"compiler\": " << json_string(build_info::compiler) << ",\n"
        << "  \"build_configuration\": " << json_string(build_info::configuration) << ",\n"
        << "  \"cgal_version\": " << json_string(CGAL_VERSION_STR) << ",\n"
        << "  \"input_path\": " << json_string(cfg.input_path) << ",\n"
        << "  \"input_sha1\": " << checksum << ",\n"
        << "  \"output_path\": " << json_string(cfg.output_path) << ",\n"
        << "  \"field_prefix\": " << json_string(cfg.export_field_prefix.empty() ? make_auto_field_prefix(cfg) : cfg.export_field_prefix) << ",\n"
        << "  \"directory\": " << json_string(directory_.generic_string()) << ",\n"
        << "  \"backend\": \"custom\",\n"
        << "  \"diagnostic_level\": " << json_string(cfg.diagnostic_level == DiagnosticLevel::Operation ? "operation" : "stage") << ",\n"
        << "  \"field\": " << json_string(field_type_name(cfg.field_type)) << ",\n"
        << "  \"epsilon\": " << cfg.epsilon << ",\n"
        << "  \"epsilon_is_explicit\": " << (cfg.epsilon_is_explicit ? "true" : "false") << ",\n"
        << "  \"min_edge_length\": " << cfg.min_edge_length << ",\n"
        << "  \"max_edge_length\": " << cfg.max_edge_length << ",\n"
        << "  \"chen_beta\": " << cfg.chen_beta << ",\n"
        << "  \"csf_mesh_scale\": " << cfg.csf_mesh_scale << ",\n"
        << "  \"iterations\": " << cfg.iterations << ",\n"
        << "  \"relaxation_steps\": " << cfg.relaxation_steps << ",\n"
        << "  \"do_project\": " << (cfg.do_project ? "true" : "false") << ",\n"
        << "  \"preserve_features\": " << (cfg.preserve_features ? "true" : "false") << ",\n"
        << "  \"feature_angle_degrees\": " << cfg.feature_angle_degrees << ",\n"
        << "  \"collapse_constraints\": false,\n"
        << "  \"relax_constraints\": true,\n"
        << "  \"descriptor_scope\": \"live mesh at this event; not compact OFF indices\",\n"
        << "  \"collapse_scope\": \"primary collapse including its degeneracy repairs\"\n"
        << "}\n";
    manifest.close();
    observer_.operations = cfg.diagnostic_level == DiagnosticLevel::Operation;
    observer_.callback = [this](const Mesh& mesh, const custom_remesh::RemeshEvent& event) {
        record(mesh, event);
    };
}

void RemeshDiagnostics::record(const Mesh& mesh, const custom_remesh::RemeshEvent& event)
{
    std::filesystem::path relative;
    const auto iteration = std::filesystem::path("iter_" + padded(event.iteration, 3));
    if (event.kind == "input") relative = "input.off";
    else if (event.kind == "iteration") relative = iteration / "result.off";
    else if (event.kind == "stage") relative = iteration / stage_file(event.stage);
    else if (event.kind == "sweep") relative = iteration / ("04_relax_" + padded(event.step, 3) + ".off");
    else if (event.kind == "operation") relative = iteration / event.stage / ("op_" + padded(event.step, 6) + ".off");
    else if (event.kind != "skipped") throw std::runtime_error("Unknown diagnostic event: " + event.kind);

    if (!relative.empty()) {
        const auto path = directory_ / relative;
        std::filesystem::create_directories(path.parent_path());
        // Never collect garbage on the live mesh: queued descriptors must stay valid.
        Mesh snapshot(mesh);
        snapshot.collect_garbage();
        auto out = open_output(path);
        if (!CGAL::IO::write_OFF(out, snapshot, CGAL::parameters::stream_precision(17)))
            throw std::runtime_error("Failed to write diagnostic snapshot: " + path.string());
        out.close();
    }

    events_ << "{\"event\":" << ++event_count_
        << ",\"kind\":" << json_string(event.kind)
        << ",\"iteration\":" << event.iteration
        << ",\"stage\":" << json_string(event.stage)
        << ",\"step\":" << event.step
        << ",\"attempted\":" << event.attempted
        << ",\"accepted\":" << event.accepted
        << ",\"rejected\":" << (event.attempted - event.accepted)
        << ",\"vertices\":" << mesh.number_of_vertices()
        << ",\"edges\":" << mesh.number_of_edges()
        << ",\"faces\":" << mesh.number_of_faces()
        << ",\"file\":" << (relative.empty() ? "null" : json_string(relative.generic_string()))
        << ",\"edge_before\":" << json_string(event.edge_before)
        << ",\"vertices_before\":[";
    for (std::size_t i = 0; i < event.vertices_before.size(); ++i) {
        if (i) events_ << ',';
        events_ << json_string(event.vertices_before[i]);
    }
    events_ << "],\"vertex_created\":" << json_string(event.vertex_created) << "}\n";
    events_.flush();
}

void RemeshDiagnostics::write_status(const char* state)
{
    auto out = open_output(directory_ / "status.json");
    out << "{\"state\":" << json_string(state) << ",\"events\":" << event_count_ << "}\n";
    out.close();
}

void RemeshDiagnostics::finish()
{
    events_.close();
    write_status("complete");
}
} // namespace rar
