#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rar::detail {

struct ChenStageResult {
    bool accepted;
    int status;
    std::vector<double> solution;
};

struct ChenContinuationResult {
    std::vector<double> solution;
    double completed_fraction = 0.0;
    double failed_fraction = 0.0;
    int failed_status = 0;
    bool used_fallback = false;
};

// Own the last validated solution: a failed solve may overwrite all solver
// buffers. Validation may normalize roundoff at bounds, but must reject an
// infeasible field before it can become a checkpoint.
template <typename Solve, typename Validate>
ChenContinuationResult run_chen_continuation(
    double fraction,
    const double growth,
    const std::vector<double>& initial,
    Solve solve,
    Validate validate)
{
    ChenContinuationResult result;
    bool first = true;
    while (true) {
        ChenStageResult stage = solve(
            fraction, first ? initial : result.solution, first);
        if (!stage.accepted) {
            if (first) {
                throw std::runtime_error(
                    "IPOPT failed in the first Chen continuation stage "
                    "(status " + std::to_string(stage.status) +
                    "); no successful stage is available for fallback");
            }
            result.used_fallback = true;
            result.failed_fraction = fraction;
            result.failed_status = stage.status;
            break;
        }

        validate(stage.solution, fraction);
        result.solution = std::move(stage.solution);
        result.completed_fraction = fraction;
        if (fraction == 1.0) {
            break;
        }
        fraction = (std::min)(1.0, fraction * growth);
        first = false;
    }
    validate(result.solution, result.completed_fraction);
    return result;
}

} // namespace rar::detail
