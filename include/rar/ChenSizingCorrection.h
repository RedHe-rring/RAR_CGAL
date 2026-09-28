#pragma once

#include "rar/Types.h"

#include <cstddef>

namespace rar {

struct ChenCorrectionStats {
    std::size_t constraint_count = 0;
    std::size_t changed_vertex_count = 0;
    double objective = 0.0;
    double max_gradient_before = 0.0;
    double max_gradient_after = 0.0;
};

ChenCorrectionStats apply_chen_sizing_correction(
    const Mesh& mesh,
    Mesh::Property_map<Mesh::Vertex_index, double> sizing_map,
    double beta,
    double hmin);

} // namespace rar
