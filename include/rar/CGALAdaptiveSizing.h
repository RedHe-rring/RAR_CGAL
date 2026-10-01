#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rar {

// Curvature-radius-guarded variant of CGAL's adaptive target-length formula.
//
// Let r = 1 / kappa_max. For epsilon <= r / 2, retain the CGAL-style formula
// without imposing r as a lower bound. Once epsilon exceeds half the local
// radius, the tolerance is considered too large for the local geometric scale
// and the target falls back to r:
//
//   h = sqrt(6 epsilon r - 3 epsilon^2),  epsilon <= r / 2
//   h = r,                               epsilon >  r / 2
//
// This is intentionally discontinuous at epsilon = r / 2: the value at the
// threshold is 3 r / 2, while values immediately above it fall back to r. The
// configured global edge-length bounds are applied last.
inline double cgal_adaptive_radius_target_length(
    const double kappa_max,
    const double epsilon,
    const double min_edge_length,
    const double max_edge_length)
{
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        throw std::invalid_argument(
            "CGAL adaptive epsilon must be finite and > 0");
    }
    if (!(min_edge_length > 0.0) ||
        !std::isfinite(min_edge_length) ||
        !std::isfinite(max_edge_length) ||
        min_edge_length > max_edge_length) {
        throw std::invalid_argument(
            "Invalid CGAL adaptive edge-length bounds");
    }

    if (!std::isfinite(kappa_max) || kappa_max <= 1e-15) {
        return max_edge_length;
    }

    const double radius = 1.0 / kappa_max;
    double target_length = radius;

    if (epsilon <= 0.5 * radius) {
        const double value_sq =
            6.0 * epsilon * radius -
            3.0 * epsilon * epsilon;
        target_length =
            std::sqrt((std::max)(0.0, value_sq));
    }

    return (std::max)(
        min_edge_length,
        (std::min)(max_edge_length, target_length));
}

} // namespace rar
