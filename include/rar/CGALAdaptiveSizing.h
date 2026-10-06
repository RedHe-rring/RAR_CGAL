#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rar {

inline constexpr double kCGALAdaptiveRadiusThresholdRatio = 0.183;

// Curvature-radius-guarded variant of CGAL's adaptive target-length formula.
//
// Let r = 1 / kappa_max. For epsilon <= 0.183 r, retain the CGAL-style
// formula without imposing r as a lower bound. Once epsilon exceeds that
// threshold, the target falls back to r:
//
//   h = sqrt(6 epsilon r - 3 epsilon^2),  epsilon <= 0.183 r
//   h = r,                               epsilon >  0.183 r
//
// The literal 0.183 ratio is a rounded version of the point where the formula
// equals r, so only a small jump remains at the threshold. The configured
// global edge-length bounds are applied last.
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

    if (epsilon <=
        kCGALAdaptiveRadiusThresholdRatio * radius) {
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
