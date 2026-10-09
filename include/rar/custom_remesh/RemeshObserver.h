#pragma once

#include <cstddef>
#include <functional>
#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>

namespace rar { namespace custom_remesh {

struct RemeshEvent {
    std::string kind;
    std::string stage;
    unsigned int iteration = 0;
    std::size_t step = 0;
    std::size_t attempted = 0;
    std::size_t accepted = 0;
    std::string edge_before;
    std::vector<std::string> vertices_before;
    std::string vertex_created;
};

// Optional, synchronous observation of completed transactions. No mesh mutation.
template<class PolygonMesh>
struct RemeshObserver {
    bool operations = false;
    unsigned int iteration = 0;
    std::string stage;
    std::size_t attempted = 0;
    std::size_t accepted = 0;
    std::function<void(const PolygonMesh&, const RemeshEvent&)> callback;

    void begin_stage(const char* name)
    {
        stage = name;
        attempted = accepted = 0;
    }

    void emit(const PolygonMesh& mesh, const char* kind, std::size_t step = 0) const
    {
        if (callback) callback(mesh, {kind, stage, iteration, step, attempted, accepted, {}, {}, {}});
    }

    template<class Descriptor>
    static std::string describe(const Descriptor& descriptor)
    {
        std::ostringstream out;
        out << descriptor;
        return out.str();
    }

    template<class Edge, class Vertex>
    void operation(const PolygonMesh& mesh, const Edge& edge_before,
                   std::initializer_list<Vertex> vertices_before,
                   const std::string& created = {})
    {
        ++accepted;
        if (!operations || !callback) return;
        RemeshEvent event{"operation", stage, iteration, accepted, attempted, accepted, {}, {}, created};
        event.edge_before = describe(edge_before);
        for (const auto& vertex : vertices_before) event.vertices_before.push_back(describe(vertex));
        callback(mesh, event);
    }
};

} } // namespace rar::custom_remesh
