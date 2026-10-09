#pragma once

#include "rar/RemeshConfig.h"
#include "rar/RemeshDiagnostics.h"
#include "rar/custom_remesh/Remesh.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <iostream>
#include <stdexcept>

namespace rar {

template<class FaceRange, class SizingField,
         class PolygonMesh, class NamedParameters>
void dispatch_remeshing(
    const FaceRange& face_range,
    SizingField& sizing,
    PolygonMesh& mesh,
    const RemeshBackend backend,
    const NamedParameters& np,
    custom_remesh::RemeshObserver<PolygonMesh>* observer = nullptr)
{
    if (backend == RemeshBackend::Custom) {
        custom_remesh::isotropic_remeshing(face_range, sizing, mesh, np, observer);
    } else {
        if (observer) throw std::invalid_argument("Remesh diagnostics require --remesh-backend custom");
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            face_range, sizing, mesh, np);
    }
}

template<class FaceRange, class SizingField, class NamedParameters>
void dispatch_remeshing(const FaceRange& face_range, SizingField& sizing,
                        Mesh& mesh, const RemeshConfig& cfg, const NamedParameters& np)
{
    if (!cfg.remesh_diagnostics) {
        dispatch_remeshing(face_range, sizing, mesh, cfg.remesh_backend, np);
        return;
    }
    if (cfg.remesh_backend != RemeshBackend::Custom)
        throw std::invalid_argument("Remesh diagnostics require --remesh-backend custom");
    RemeshDiagnostics diagnostics(cfg);
    std::cout << "Remesh diagnostics: " << diagnostics.directory().string() << '\n';
    dispatch_remeshing(face_range, sizing, mesh, cfg.remesh_backend, np, &diagnostics.observer());
    diagnostics.finish();
}

} // namespace rar
