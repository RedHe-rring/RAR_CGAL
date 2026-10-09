#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace rar {

inline constexpr double kSDFThicknessEpsilonRatio = 0.25;

struct SDFThicknessSizingSample {
    double epsilon_local = 0.0;
    double target_length = 0.0;
    bool epsilon_limited = false;
    bool radius_limited = false;
    bool thickness_limited = false;
    bool min_edge_clipped = false;
    bool max_edge_clipped = false;
};

// kappa = max(abs(k1), abs(k2)); thickness and epsilon use mesh length units.
// This defines an initial sizing field, not a guarantee against self-intersection.
inline SDFThicknessSizingSample sdf_thickness_target_length(
    double kappa, double thickness, double epsilon, double min_edge, double max_edge)
{
    if (!std::isfinite(kappa) || kappa < 0.0)
        throw std::invalid_argument("SDF sizing curvature must be finite and >= 0");
    if (!std::isfinite(thickness) || !(thickness > 0.0))
        throw std::invalid_argument("SDF thickness must be finite and > 0");
    if (!std::isfinite(epsilon) || !(epsilon > 0.0))
        throw std::invalid_argument("SDF epsilon must be finite and > 0");
    if (!std::isfinite(min_edge) || !std::isfinite(max_edge) ||
        !(min_edge > 0.0) || min_edge > max_edge)
        throw std::invalid_argument("Invalid SDF edge-length bounds");

    const double radius_epsilon = kappa == 0.0
        ? std::numeric_limits<double>::infinity() : 0.183 / kappa;
    const double thickness_epsilon = thickness * kSDFThicknessEpsilonRatio;
    SDFThicknessSizingSample result;
    result.epsilon_local = (std::min)({epsilon, radius_epsilon, thickness_epsilon});
    if (!(result.epsilon_local > 0.0))
        throw std::runtime_error("SDF local epsilon underflow; rescale the input mesh");
    // Ties intentionally activate all equal limiting terms.
    result.epsilon_limited = result.epsilon_local == epsilon;
    result.radius_limited = result.epsilon_local == radius_epsilon;
    result.thickness_limited = result.epsilon_local == thickness_epsilon;

    // Algebraically sqrt(6*epsilon_local/kappa - 3*epsilon_local^2).
    // Factoring avoids overflow from squaring large lengths. The radius cap
    // keeps kappa*epsilon_local <= 0.183, so the radicand is positive.
    const double h = kappa == 0.0 ? std::numeric_limits<double>::infinity()
        : std::sqrt(result.epsilon_local) / std::sqrt(kappa)
            * std::sqrt(6.0 - 3.0 * kappa * result.epsilon_local);
    result.min_edge_clipped = h < min_edge;
    result.max_edge_clipped = h > max_edge;
    result.target_length = (std::max)(min_edge, (std::min)(max_edge, h));
    return result;
}

} // namespace rar
