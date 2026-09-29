#include "rar/OutputNaming.h"
#include "rar/RARSizingField.h"
#include "rar/RemeshConfig.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <string>

int main() {
    rar::RemeshConfig cfg;
    assert(
        cfg.field_type ==
        rar::FieldType::CGALAdaptive);
    assert(cfg.epsilon > 0.0);
    assert(cfg.min_edge_length > 0.0);
    assert(
        cfg.max_edge_length >=
        cfg.min_edge_length);
    assert(cfg.chen_beta > 1.0);
    assert(cfg.csf_mesh_scale > 0.0);
    assert(cfg.iterations > 0);

    const double kappa = 2.0;
    const double epsilon = 0.01;
    const double expected =
        std::sqrt(
            6.0 * epsilon / kappa -
            3.0 * epsilon * epsilon);

    const double actual =
        rar::rar_target_length(
            kappa,
            epsilon,
            1e-6,
            10.0);

    assert(
        std::abs(actual - expected) <
        1e-12);

    const double eta = 0.01;
    const double kappa_a = 2.0;
    const double kappa_b = 5.0;
    const double relative_a =
        rar::rar_curvature_normalized_target_length(
            kappa_a,
            eta,
            1e-9,
            10.0);
    const double relative_b =
        rar::rar_curvature_normalized_target_length(
            kappa_b,
            eta,
            1e-9,
            10.0);

    const double expected_relative_factor =
        std::sqrt(
            6.0 * eta -
            3.0 * eta * eta);

    assert(
        std::abs(
            relative_a * kappa_a -
            expected_relative_factor) <
        1e-12);
    assert(
        std::abs(
            relative_b * kappa_b -
            expected_relative_factor) <
        1e-12);

    const double flat_relative =
        rar::rar_curvature_normalized_target_length(
            0.0,
            eta,
            0.1,
            2.0);
    assert(
        std::abs(flat_relative - 2.0) <
        1e-12);

    const double flat =
        rar::rar_target_length(
            0.0,
            epsilon,
            0.1,
            2.0);

    assert(
        std::abs(flat - 2.0) <
        1e-12);

    const double clamped_small =
        rar::rar_target_length(
            1e12,
            epsilon,
            0.1,
            2.0);

    assert(
        std::abs(clamped_small - 0.1) <
        1e-12);

    rar::RemeshConfig naming;
    naming.input_path =
        "models/sample.obj";
    naming.field_type =
        rar::FieldType::RAR;
    naming.epsilon = 0.001;
    naming.min_edge_length = 0.001;
    naming.max_edge_length = 0.5;
    naming.chen_beta = 1.2;
    naming.iterations = 5;
    naming.relaxation_steps = 3;
    naming.do_project = true;

    const std::string rar_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        rar_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar.obj");

    naming.output_path =
        rar::make_auto_output_path(naming);

    const std::filesystem::path
        auto_field_prefix =
            rar::make_auto_field_prefix(
                naming);

    assert(
        auto_field_prefix
            .filename()
            .string() ==
        "field");

    assert(
        auto_field_prefix
            .parent_path()
            .filename()
            .string() ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar");

    naming.field_type =
        rar::FieldType::RARRelative;
    naming.relative_error = 0.01;

    const std::string rar_relative_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        rar_relative_auto_name ==
        "sample__eta-0p01__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar-relative.obj");

    naming.field_type =
        rar::FieldType::CGALAdaptiveChen;

    const std::string cgal_chen_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        cgal_chen_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__beta-1p2__it-5__relax-3__proj-on__field-cgal-adaptive-chen.obj");

    naming.field_type =
        rar::FieldType::RARChen;

    const std::string rar_chen_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        rar_chen_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__beta-1p2__it-5__relax-3__proj-on__field-rar-chen.obj");

    naming.field_type =
        rar::FieldType::CSF;
    naming.csf_mesh_scale = 1.2;

    const std::string csf_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        csf_auto_name ==
        "sample__scale-1p2__it-5__relax-3__proj-on__field-csf.obj");

    assert(naming.export_field);

    return 0;
}
