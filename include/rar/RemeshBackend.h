#pragma once

#include "rar/RemeshConfig.h"
#include "rar/custom_remesh/Remesh.h"

#include <CGAL/Polygon_mesh_processing/remesh.h>

namespace rar {

template<class FaceRange, class SizingField,
         class PolygonMesh, class NamedParameters>
void dispatch_remeshing(
    const FaceRange& face_range,
    SizingField& sizing,
    PolygonMesh& mesh,
    const RemeshBackend backend,
    const NamedParameters& np)
{
    if (backend == RemeshBackend::Custom) {
        custom_remesh::isotropic_remeshing(face_range, sizing, mesh, np);
    } else {
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            face_range, sizing, mesh, np);
    }
}

} // namespace rar
