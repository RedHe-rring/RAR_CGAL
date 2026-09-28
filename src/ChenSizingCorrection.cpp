#include "rar/ChenSizingCorrection.h"

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(RAR_HAS_IPOPT)
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#endif

namespace rar {

#if defined(RAR_HAS_IPOPT)

namespace {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 to_vec3(const Kernel::Vector_3& v) {
    return {
        CGAL::to_double(v.x()),
        CGAL::to_double(v.y()),
        CGAL::to_double(v.z())
    };
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct TriangleConstraint {
    std::array<std::size_t, 3> vertex_ids{};
    std::array<Vec3, 3> basis_gradients{};
};

inline double constraint_value(
    const TriangleConstraint& c,
    const double* x)
{
    Vec3 grad;

    for (std::size_t local = 0; local < 3; ++local) {
        const double value = x[c.vertex_ids[local]];
        grad.x += value * c.basis_gradients[local].x;
        grad.y += value * c.basis_gradients[local].y;
        grad.z += value * c.basis_gradients[local].z;
    }

    return dot(grad, grad);
}

class ChenProjectionNLP final : public Ipopt::TNLP {
public:
    ChenProjectionNLP(
        std::vector<double> h0,
        std::vector<TriangleConstraint> constraints,
        const double beta,
        const double hmin)
        : h0_(std::move(h0)),
          constraints_(std::move(constraints)),
          hmin_(hmin),
          gradient_limit_sq_(
              std::log(beta) * std::log(beta)),
          start_(h0_),
          solution_(h0_)
    {
        // The feasible set is never empty because the constant field h=hmin
        // satisfies every gradient constraint and the variable bounds.
        //
        // Starting IPOPT directly from h0 can be strongly infeasible. Build a
        // feasible warm start by uniformly shrinking the variation of h0
        // around the constant hmin field. Since grad(hmin)=0,
        //
        //   grad(hmin + s * (h0 - hmin)) = s * grad(h0).
        //
        // Choosing s from the maximum initial gradient therefore guarantees
        // the triangle gradient bounds up to a small safety margin.
        double max_gradient_sq = 0.0;
        for (const TriangleConstraint& c : constraints_) {
            max_gradient_sq =
                (std::max)(
                    max_gradient_sq,
                    constraint_value(c, h0_.data()));
        }

        const double max_gradient =
            std::sqrt((std::max)(0.0, max_gradient_sq));
        const double gradient_limit =
            std::sqrt(gradient_limit_sq_);

        if (max_gradient > gradient_limit &&
            max_gradient > 0.0) {
            const double scale =
                0.95 * gradient_limit / max_gradient;

            for (std::size_t i = 0; i < start_.size(); ++i) {
                start_[i] =
                    hmin_ +
                    scale * (h0_[i] - hmin_);
            }
        }
    }

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style) override
    {
        n = static_cast<Ipopt::Index>(h0_.size());
        m = static_cast<Ipopt::Index>(constraints_.size());
        nnz_jac_g = 3 * m;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u) override
    {
        if (static_cast<std::size_t>(n) != h0_.size() ||
            static_cast<std::size_t>(m) != constraints_.size()) {
            return false;
        }

        for (Ipopt::Index i = 0; i < n; ++i) {
            x_l[i] = hmin_;
            x_u[i] = h0_[static_cast<std::size_t>(i)];
        }

        for (Ipopt::Index i = 0; i < m; ++i) {
            // This constraint has only an upper bound:
            // ||grad h||^2 <= log(beta)^2.
            // IPOPT treats bounds <= -1e19 as unbounded by default.
            g_l[i] = -1e19;
            g_u[i] = gradient_limit_sq_;
        }

        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index,
        bool init_lambda,
        Ipopt::Number*) override
    {
        if (!init_x || init_z || init_lambda ||
            static_cast<std::size_t>(n) != h0_.size()) {
            return false;
        }

        for (Ipopt::Index i = 0; i < n; ++i) {
            x[i] = start_[static_cast<std::size_t>(i)];
        }

        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& obj_value) override
    {
        if (static_cast<std::size_t>(n) != h0_.size()) {
            return false;
        }

        obj_value = 0.0;

        for (Ipopt::Index i = 0; i < n; ++i) {
            const double d =
                x[i] - h0_[static_cast<std::size_t>(i)];
            obj_value += d * d;
        }

        return true;
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* grad_f) override
    {
        if (static_cast<std::size_t>(n) != h0_.size()) {
            return false;
        }

        for (Ipopt::Index i = 0; i < n; ++i) {
            grad_f[i] =
                2.0 *
                (x[i] - h0_[static_cast<std::size_t>(i)]);
        }

        return true;
    }

    bool eval_g(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* g) override
    {
        if (static_cast<std::size_t>(m) != constraints_.size()) {
            return false;
        }

        for (Ipopt::Index i = 0; i < m; ++i) {
            g[i] = constraint_value(
                constraints_[static_cast<std::size_t>(i)],
                x);
        }

        return true;
    }

    bool eval_jac_g(
        Ipopt::Index,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nele_jac,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values) override
    {
        if (static_cast<std::size_t>(m) != constraints_.size() ||
            nele_jac != 3 * m) {
            return false;
        }

        Ipopt::Index cursor = 0;

        if (values == nullptr) {
            for (Ipopt::Index i = 0; i < m; ++i) {
                const auto& c =
                    constraints_[static_cast<std::size_t>(i)];

                for (std::size_t local = 0; local < 3; ++local) {
                    iRow[cursor] = i;
                    jCol[cursor] =
                        static_cast<Ipopt::Index>(
                            c.vertex_ids[local]);
                    ++cursor;
                }
            }

            return true;
        }

        for (Ipopt::Index i = 0; i < m; ++i) {
            const auto& c =
                constraints_[static_cast<std::size_t>(i)];

            Vec3 grad;
            for (std::size_t local = 0; local < 3; ++local) {
                const double value = x[c.vertex_ids[local]];
                grad.x += value * c.basis_gradients[local].x;
                grad.y += value * c.basis_gradients[local].y;
                grad.z += value * c.basis_gradients[local].z;
            }

            for (std::size_t local = 0; local < 3; ++local) {
                values[cursor++] =
                    2.0 * dot(grad, c.basis_gradients[local]);
            }
        }

        return true;
    }

    bool eval_h(
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index,
        Ipopt::Index*,
        Ipopt::Index*,
        Ipopt::Number*) override
    {
        // IPOPT uses its limited-memory Hessian approximation.
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Index,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Number,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*) override
    {
        solution_.assign(
            x,
            x + static_cast<std::size_t>(n));
    }

    const std::vector<double>& solution() const {
        return solution_;
    }

private:
    std::vector<double> h0_;
    std::vector<TriangleConstraint> constraints_;
    double hmin_;
    double gradient_limit_sq_;
    std::vector<double> start_;
    std::vector<double> solution_;

};

std::vector<TriangleConstraint> build_constraints(
    const Mesh& mesh,
    const std::map<Mesh::Vertex_index, std::size_t>& dense_index)
{
    std::vector<TriangleConstraint> result;
    result.reserve(num_faces(mesh));

    for (const auto f : faces(mesh)) {
        const auto h = halfedge(f, mesh);
        const auto v0 = source(h, mesh);
        const auto v1 = target(h, mesh);
        const auto v2 = target(next(h, mesh), mesh);

        const Point& p0 = mesh.point(v0);
        const Point& p1 = mesh.point(v1);
        const Point& p2 = mesh.point(v2);

        const Kernel::Vector_3 normal =
            CGAL::cross_product(
                Kernel::Vector_3(p0, p1),
                Kernel::Vector_3(p0, p2));

        const double normal_sq =
            CGAL::to_double(normal.squared_length());

        if (!(normal_sq > 1e-30) ||
            !std::isfinite(normal_sq)) {
            continue;
        }

        const Kernel::Vector_3 grad0 =
            CGAL::cross_product(
                normal,
                Kernel::Vector_3(p1, p2)) /
            normal_sq;

        const Kernel::Vector_3 grad1 =
            CGAL::cross_product(
                normal,
                Kernel::Vector_3(p2, p0)) /
            normal_sq;

        const Kernel::Vector_3 grad2 =
            CGAL::cross_product(
                normal,
                Kernel::Vector_3(p0, p1)) /
            normal_sq;

        TriangleConstraint c;
        c.vertex_ids = {
            dense_index.at(v0),
            dense_index.at(v1),
            dense_index.at(v2)
        };
        c.basis_gradients = {
            to_vec3(grad0),
            to_vec3(grad1),
            to_vec3(grad2)
        };

        result.push_back(c);
    }

    return result;
}

double max_gradient(
    const std::vector<TriangleConstraint>& constraints,
    const std::vector<double>& h)
{
    double max_sq = 0.0;

    for (const auto& c : constraints) {
        max_sq =
            (std::max)(
                max_sq,
                constraint_value(c, h.data()));
    }

    return std::sqrt((std::max)(0.0, max_sq));
}

} // namespace

ChenCorrectionStats apply_chen_sizing_correction(
    const Mesh& mesh,
    Mesh::Property_map<Mesh::Vertex_index, double> sizing_map,
    const double beta,
    const double hmin)
{
    if (!(beta > 1.0) || !std::isfinite(beta)) {
        throw std::invalid_argument(
            "Chen progressive factor beta must be > 1");
    }

    if (!(hmin > 0.0) || !std::isfinite(hmin)) {
        throw std::invalid_argument(
            "Chen hmin must be finite and > 0");
    }

    std::map<Mesh::Vertex_index, std::size_t> dense_index;
    std::vector<Mesh::Vertex_index> dense_vertices;
    dense_vertices.reserve(num_vertices(mesh));

    std::vector<double> h0;
    h0.reserve(num_vertices(mesh));

    std::size_t next_id = 0;
    for (const auto v : vertices(mesh)) {
        const double value = sizing_map[v];

        if (!std::isfinite(value) || value < hmin) {
            throw std::runtime_error(
                "RAR sizing field contains an invalid value for Chen correction");
        }

        dense_index.emplace(v, next_id++);
        dense_vertices.push_back(v);
        h0.push_back(value);
    }

    const std::vector<TriangleConstraint> constraints =
        build_constraints(mesh, dense_index);

    ChenCorrectionStats stats;
    stats.constraint_count = constraints.size();
    stats.max_gradient_before =
        max_gradient(constraints, h0);

    if (h0.empty() || constraints.empty()) {
        stats.max_gradient_after =
            stats.max_gradient_before;
        return stats;
    }

    Ipopt::SmartPtr<ChenProjectionNLP> problem =
        new ChenProjectionNLP(
            h0,
            constraints,
            beta,
            hmin);

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app =
        IpoptApplicationFactory();

    if (!app->Options()->SetStringValue("linear_solver", "mumps") ||
        !app->Options()->SetIntegerValue("print_level", 0) ||
        !app->Options()->SetStringValue("sb", "yes") ||
        !app->Options()->SetStringValue(
            "hessian_approximation",
            "limited-memory") ||
        !app->Options()->SetNumericValue("tol", 1e-8) ||
        !app->Options()->SetNumericValue("constr_viol_tol", 1e-8) ||
        !app->Options()->SetIntegerValue("max_iter", 1000)) {
        throw std::runtime_error(
            "Failed to configure IPOPT for Chen correction");
    }

    const Ipopt::ApplicationReturnStatus init_status =
        app->Initialize();

    if (init_status != Ipopt::Solve_Succeeded) {
        throw std::runtime_error(
            "Failed to initialize IPOPT for Chen correction");
    }

    const Ipopt::ApplicationReturnStatus solve_status =
        app->OptimizeTNLP(problem);

    if (solve_status != Ipopt::Solve_Succeeded &&
        solve_status != Ipopt::Solved_To_Acceptable_Level &&
        solve_status != Ipopt::Feasible_Point_Found) {
        throw std::runtime_error(
            "IPOPT failed while solving the Chen sizing correction (status " +
            std::to_string(static_cast<int>(solve_status)) + ")");
    }

    const std::vector<double>& solution =
        problem->solution();

    if (solution.size() != dense_vertices.size()) {
        throw std::runtime_error(
            "IPOPT returned a Chen solution with an unexpected size");
    }

    for (std::size_t i = 0; i < solution.size(); ++i) {
        const double value =
            (std::max)(
                hmin,
                (std::min)(h0[i], solution[i]));

        sizing_map[dense_vertices[i]] = value;

        const double scale =
            (std::max)(1.0, std::abs(h0[i]));
        if (std::abs(value - h0[i]) > 1e-10 * scale) {
            ++stats.changed_vertex_count;
        }

        const double d = value - h0[i];
        stats.objective += d * d;
    }

    std::vector<double> corrected(solution.size(), 0.0);
    for (std::size_t i = 0; i < dense_vertices.size(); ++i) {
        corrected[i] = sizing_map[dense_vertices[i]];
    }

    stats.max_gradient_after =
        max_gradient(constraints, corrected);

    return stats;
}

#else

ChenCorrectionStats apply_chen_sizing_correction(
    const Mesh&,
    Mesh::Property_map<Mesh::Vertex_index, double>,
    double,
    double)
{
    throw std::runtime_error(
        "RAR+Chen requires IPOPT. Install coin-or-ipopt and reconfigure CMake.");
}

#endif

} // namespace rar
