#include "rar/OutputNaming.h"
#include "rar/RARSizingField.h"
#include "rar/RemeshConfig.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <string>

int main() {
    rar::RemeshConfig cfg;
    assert(cfg.field_type == rar::FieldType::CGALAdaptive);
    assert(cfg.epsilon > 0.0);
    assert(cfg.min_edge_length > 0.0);
    assert(cfg.max_edge_length >= cfg.min_edge_length);
    assert(cfg.iterations > 0);

    const double kappa = 2.0;
    const double epsilon = 0.01;
    const double expected =
        std::sqrt(6.0 * epsilon / kappa -
                  3.0 * epsilon * epsilon);

    const double actual =
        rar::rar_target_length(kappa, epsilon, 1e-6, 10.0);

    assert(std::abs(actual - expected) < 1e-12);

    const double flat =
        rar::rar_target_length(0.0, epsilon, 0.1, 2.0);
    assert(std::abs(flat - 2.0) < 1e-12);

    const double clamped_small =
        rar::rar_target_length(1e12, epsilon, 0.1, 2.0);
    assert(std::abs(clamped_small - 0.1) < 1e-12);

    rar::RemeshConfig naming;
    naming.input_path = "models/sample.obj";
    naming.field_type = rar::FieldType::RAR;
    naming.epsilon = 0.001;
    naming.min_edge_length = 0.001;
    naming.max_edge_length = 0.5;
    naming.iterations = 5;
    naming.relaxation_steps = 3;
    naming.do_project = true;

    const std::string auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(naming))
            .filename()
            .string();

    assert(
        auto_name ==
        "sample__field-rar__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on.obj");

    return 0;
}
