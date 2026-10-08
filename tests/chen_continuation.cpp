#include "rar/detail/ChenContinuation.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(const bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void validate(std::vector<double>& values, const double fraction) {
    require(values.size() == 1, "invalid solution size");
    require(std::isfinite(values[0]) && values[0] <= fraction,
            "invalid solution must be rejected");
}

void test_complete() {
    std::vector<double> fractions;
    const auto result = rar::detail::run_chen_continuation(
        0.02, 10.0, {0.0},
        [&](const double f, const std::vector<double>& start, const bool first) {
            require(first == fractions.empty(), "incorrect first-stage flag");
            require(start[0] == (first ? 0.0 : fractions.back()),
                    "stage must start from previous success");
            fractions.push_back(f);
            return rar::detail::ChenStageResult{true, 0, {f}};
        }, validate);
    require(fractions.size() == 3 && fractions.back() == 1.0,
            "continuation must reach exactly f=1");
    require(!result.used_fallback && result.completed_fraction == 1.0 &&
            result.solution[0] == 1.0 && result.failed_fraction == 0.0,
            "successful continuation reported incorrectly");
}

void test_fallback() {
    int calls = 0;
    std::vector<double> solver_buffer;
    const auto result = rar::detail::run_chen_continuation(
        0.01, 10.0, {0.0},
        [&](const double f, const std::vector<double>& start, const bool) {
            ++calls;
            require(start[0] == (calls == 1 ? 0.0 : 0.01),
                    "fallback test warm start is wrong");
            solver_buffer = calls == 1
                ? std::vector<double>{f}
                : std::vector<double>{std::numeric_limits<double>::quiet_NaN()};
            return rar::detail::ChenStageResult{
                calls == 1, calls == 1 ? 0 : -1, solver_buffer};
        }, validate);
    require(calls == 2, "failed stage must stop continuation");
    require(result.used_fallback && result.completed_fraction == 0.01 &&
            result.failed_fraction == 0.1 && result.failed_status == -1,
            "fallback metadata is wrong");
    require(result.solution[0] == 0.01,
            "failed solver iterate overwrote successful checkpoint");
}

void test_first_stage_failure() {
    bool threw = false;
    try {
        rar::detail::run_chen_continuation(
            0.01, 10.0, {0.0},
            [](double, const std::vector<double>&, bool) {
                return rar::detail::ChenStageResult{false, -1, {0.0}};
            }, validate);
    } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).find("no successful stage") != std::string::npos;
    }
    require(threw, "first failure must not turn the seed into a fallback");
}

void test_reject_invalid_success() {
    bool threw = false;
    try {
        rar::detail::run_chen_continuation(
            0.01, 10.0, {0.0},
            [](double, const std::vector<double>&, bool) {
                return rar::detail::ChenStageResult{true, 0, {10.0}};
            }, validate);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    require(threw, "accepted solver status must still pass validation");
}

} // namespace

int main() {
    try {
        test_complete();
        test_fallback();
        test_first_stage_failure();
        test_reject_invalid_success();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
