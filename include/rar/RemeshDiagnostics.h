#pragma once

#include "rar/RemeshConfig.h"
#include "rar/Types.h"
#include "rar/custom_remesh/RemeshObserver.h"

#include <filesystem>
#include <fstream>
#include <string>

namespace rar {

// One recorder owns one exclusively created run directory.
class RemeshDiagnostics {
public:
    explicit RemeshDiagnostics(const RemeshConfig& cfg);
    RemeshDiagnostics(const RemeshDiagnostics&) = delete;
    RemeshDiagnostics& operator=(const RemeshDiagnostics&) = delete;

    custom_remesh::RemeshObserver<Mesh>& observer() { return observer_; }
    const std::filesystem::path& directory() const { return directory_; }
    void finish();

private:
    void record(const Mesh& mesh, const custom_remesh::RemeshEvent& event);
    void write_status(const char* state);

    std::filesystem::path directory_;
    std::ofstream events_;
    std::size_t event_count_ = 0;
    custom_remesh::RemeshObserver<Mesh> observer_;
};

} // namespace rar
