#include "rar/CGALAdaptiveSizing.h"
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
    assert(!cfg.epsilon_is_explicit);
    assert(cfg.min_edge_length > 0.0);
    assert(
        cfg.max_edge_length >=
        cfg.min_edge_length);
    assert(cfg.chen_beta > 1.0);
    assert(cfg.csf_mesh_scale > 0.0);
    assert(cfg.iterations > 0);
    assert(!cfg.preserve_features);
    assert(cfg.feature_angle_degrees == 50.0);

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

    const double below_radius =
        rar::cgal_adaptive_radius_target_length(
            2.0,
            0.01,
            1e-6,
            10.0);
    const double expected_below_radius =
        std::sqrt(6.0 * 0.01 * 0.5 - 3.0 * 0.01 * 0.01);
    assert(
        std::abs(below_radius - expected_below_radius) <
        1e-12);
    assert(below_radius < 0.5);

    const double radius = 0.5;
    const double radius_threshold =
        rar::kCGALAdaptiveRadiusThresholdRatio * radius;

    const double below_radius_threshold =
        rar::cgal_adaptive_radius_target_length(
            2.0,
            0.09,
            1e-6,
            10.0);
    const double expected_below_radius_threshold =
        std::sqrt(6.0 * 0.09 * radius - 3.0 * 0.09 * 0.09);
    assert(
        std::abs(
            below_radius_threshold -
            expected_below_radius_threshold) <
        1e-12);

    const double at_radius_threshold =
        rar::cgal_adaptive_radius_target_length(
            2.0,
            radius_threshold,
            1e-6,
            10.0);
    const double expected_at_radius_threshold =
        std::sqrt(
            6.0 * radius_threshold * radius -
            3.0 * radius_threshold * radius_threshold);
    assert(
        std::abs(
            at_radius_threshold -
            expected_at_radius_threshold) <
        1e-12);

    const double above_radius_threshold =
        rar::cgal_adaptive_radius_target_length(
            2.0,
            0.1,
            1e-6,
            10.0);
    assert(
        std::abs(above_radius_threshold - radius) <
        1e-12);

    const double flat_radius_field =
        rar::cgal_adaptive_radius_target_length(
            0.0,
            0.01,
            0.1,
            2.0);
    assert(std::abs(flat_radius_field - 2.0) < 1e-12);

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
    naming.preserve_features = false;
    naming.feature_angle_degrees = 50.0;

    const std::string rar_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        rar_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar.obj");

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
        "sample__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar");

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
        "sample__eps-0p001__lmin-0p001__lmax-0p5__beta-1p2__it-5__relax-3__proj-on__features-off__field-cgal-adaptive-chen.obj");

    naming.field_type =
        rar::FieldType::CGALAdaptiveRadius;

    const std::string cgal_radius_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(naming))
            .filename()
            .string();

    assert(
        cgal_radius_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-cgal-adaptive-radius.obj");

    naming.field_type =
        rar::FieldType::CGALAdaptiveRadiusChen;

    const std::string cgal_radius_chen_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        cgal_radius_chen_auto_name ==
        "sample__eps-0p001__lmin-0p001__lmax-0p5__beta-1p2__it-5__relax-3__proj-on__features-off__field-cgal-adaptive-radius-chen.obj");

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
        "sample__eps-0p001__lmin-0p001__lmax-0p5__beta-1p2__it-5__relax-3__proj-on__features-off__field-rar-chen.obj");

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
        "sample__scale-1p2__it-5__relax-3__proj-on__features-off__field-csf.obj");

    naming.preserve_features = true;
    naming.feature_angle_degrees = 50.0;

    const std::string features_auto_name =
        std::filesystem::path(
            rar::make_auto_output_path(
                naming))
            .filename()
            .string();

    assert(
        features_auto_name ==
        "sample__scale-1p2__it-5__relax-3__proj-on__features-on__feature-angle-50__field-csf.obj");

    assert(naming.export_field);

    return 0;
}
