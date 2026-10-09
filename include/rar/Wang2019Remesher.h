#pragma once
// Wang et al. (TVCG 2019) angle-driven remeshing: independent, experimental
// post-processing backend. This is an implementation study, not author code.
#include "rar/Types.h"
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#if __has_include(<CGAL/AABB_traits_3.h>)
#  include <CGAL/AABB_traits_3.h>
#  define RAR_WANG2019_CGAL_AABB_3 1
#else
#  include <CGAL/AABB_traits.h>
#  define RAR_WANG2019_CGAL_AABB_3 0
#endif
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/number_utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rar {
namespace wang2019 {
namespace PMP = CGAL::Polygon_mesh_processing;

using Vertex = Mesh::Vertex_index;
using Edge = Mesh::Edge_index;
using Face = Mesh::Face_index;
using Halfedge = Mesh::Halfedge_index;
using Vector = Kernel::Vector_3;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
#if RAR_WANG2019_CGAL_AABB_3
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
#else
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
#endif
using Tree = CGAL::AABB_tree<Traits>;

struct Options {
    double min_angle = 30.0;        // degrees, user-controllable
    double max_angle = 90.0;        // degrees, user-controllable
    double feature_angle = 50.0;    // dihedral degrees
    unsigned rounds = 10;
    unsigned smoothing_steps = 3;
    double budget_fraction = 0.02;  // at most this fraction of initial N / pass
    bool protect_features = true;
    bool project = true;
    bool verbose = true;
};

struct Angles {
    double min = 180.0;
    double max = 0.0;
};

struct Statistics {
    std::size_t insertions = 0;
    std::size_t collapses = 0;
    std::size_t flips = 0;
    std::size_t vertices_before = 0;
    std::size_t vertices_after = 0;
    std::size_t faces_before = 0;
    std::size_t faces_after = 0;
    std::size_t min_angle_violations = 0;
    std::size_t max_angle_violations = 0;
    Angles angles;
};

inline double norm2(const Vector& v) {
    return CGAL::to_double(v.squared_length());
}
inline double dot(const Vector& a, const Vector& b) {
    return CGAL::to_double(a * b);
}
inline Vector cross(const Vector& a, const Vector& b) {
    return CGAL::cross_product(a, b);
}
inline double degrees(double radians) {
    return radians * (180.0 / 3.14159265358979323846);
}
inline double angle(const Point& a, const Point& b, const Point& c) {
    const Vector u = a - b, v = c - b;
    const double denom = std::sqrt(norm2(u) * norm2(v));
    if (denom <= 1e-30) return 0.0;
    const double cosine = (std::max)(-1.0, (std::min)(1.0, dot(u, v) / denom));
    return degrees(std::acos(cosine));
}
inline Angles triangle_angles(const Point& a, const Point& b, const Point& c) {
    const double x = angle(b, a, c), y = angle(a, b, c), z = angle(a, c, b);
    return {(std::min)({x, y, z}), (std::max)({x, y, z})};
}
inline std::array<Vertex, 3> vertices_of(const Mesh& m, Face f) {
    const Halfedge h = halfedge(f, m);
    return {source(h, m), target(h, m), target(next(h, m), m)};
}
inline Angles angles_of(const Mesh& m, Face f) {
    const auto v = vertices_of(m, f);
    return triangle_angles(m.point(v[0]), m.point(v[1]), m.point(v[2]));
}
inline Vector normal(const Point& a, const Point& b, const Point& c) {
    return cross(b - a, c - a);
}
// Exact face topology is checked by Euler; this test additionally rejects local
// degenerate/inverted triangles (but it is NOT a global self-intersection check).
inline bool acceptable_triangle(const Point& a, const Point& b,
                                const Point& c, const Vector& old_n) {
    const Vector n = normal(a, b, c);
    return norm2(n) > 1e-28 &&
           dot(n, old_n) > 1e-6 * std::sqrt(norm2(n) * norm2(old_n));
}

inline std::set<Vertex> protected_vertices(const Mesh& m, double feature_angle) {
    std::set<Vertex> result;
    const double cos_limit = std::cos(feature_angle * 3.14159265358979323846 / 180.0);
    for (const Edge e : edges(m)) {
        const Halfedge h = halfedge(e, m);
        if (is_border(h, m) || is_border(opposite(h, m), m)) {
            result.insert(source(h, m));
            result.insert(target(h, m));
            continue;
        }
        const auto a = vertices_of(m, face(h, m));
        const auto b = vertices_of(m, face(opposite(h, m), m));
        const Vector na = normal(m.point(a[0]), m.point(a[1]), m.point(a[2]));
        const Vector nb = normal(m.point(b[0]), m.point(b[1]), m.point(b[2]));
        const double d = std::sqrt(norm2(na) * norm2(nb));
        if (d <= 1e-30 || dot(na, nb) / d < cos_limit) {
            result.insert(source(h, m));
            result.insert(target(h, m));
        }
    }
    return result;
}

inline bool is_protected(const std::set<Vertex>& fixed, Vertex v) {
    return fixed.find(v) != fixed.end();
}

inline bool split_quad(Mesh& m, Face f, Vertex vm) {
    // Euler::split_edge gives a quad, NOT two triangles. Connect the new
    // vertex to the opposite quad vertex with Euler::split_face.
    if (m.is_removed(f)) return false;
    const Halfedge start = halfedge(f, m);
    for (const Halfedge h : CGAL::halfedges_around_face(start, m)) {
        if (target(h, m) != vm) continue;
        const Halfedge opposite_corner = next(next(h, m), m);
        CGAL::Euler::split_face(h, opposite_corner, m);
        return true;
    }
    return false;
}

inline bool insert_at_large_angle(Mesh& m, Face f, const Tree& reference,
                                  const Options& opt, const std::set<Vertex>& fixed) {
    if (m.is_removed(f) || angles_of(m, f).max <= opt.max_angle) return false;
    // The edge opposite a triangle's largest angle is its longest edge.
    Halfedge chosen = halfedge(f, m);
    double max_len2 = -1.0;
    for (const Halfedge h : CGAL::halfedges_around_face(halfedge(f, m), m)) {
        const double d2 = norm2(m.point(source(h, m)) - m.point(target(h, m)));
        if (d2 > max_len2) { max_len2 = d2; chosen = h; }
    }
    const Halfedge ho = opposite(chosen, m);
    if (is_border(chosen, m) || is_border(ho, m)) return false;
    const Vertex va = source(chosen, m), vb = target(chosen, m);
    if (is_protected(fixed, va) || is_protected(fixed, vb)) return false;
    const Vertex vc = target(next(chosen, m), m);
    const Vertex vd = target(next(ho, m), m);
    // Moving an edge next to a feature corner is deliberately disallowed in
    // this first prototype; it is stricter than the paper's feature rules.
    if (is_protected(fixed, vc) || is_protected(fixed, vd)) return false;
    const Point a = m.point(va), b = m.point(vb);
    const Point c = m.point(vc), d = m.point(vd);
    const Point midpoint = CGAL::midpoint(a, b);
    const Point p = opt.project ? reference.closest_point(midpoint) : midpoint;
    const Vector n0 = normal(a, b, c), n1 = normal(b, a, d);
    if (!acceptable_triangle(a, p, c, n0) ||
        !acceptable_triangle(p, b, c, n0) ||
        !acceptable_triangle(b, p, d, n1) ||
        !acceptable_triangle(p, a, d, n1)) return false;
    const double before = (std::max)(triangle_angles(a, b, c).max,
                                      triangle_angles(b, a, d).max);
    const double after = (std::max)({
        triangle_angles(a, p, c).max, triangle_angles(p, b, c).max,
        triangle_angles(b, p, d).max, triangle_angles(p, a, d).max
    });
    if (after >= before - 1e-7) return false;

    const Face f0 = face(chosen, m), f1 = face(ho, m);
    const Halfedge new_h = CGAL::Euler::split_edge(chosen, m);
    const Vertex vm = target(new_h, m);
    m.point(vm) = p;
    if (!split_quad(m, f0, vm) || !split_quad(m, f1, vm))
        throw std::runtime_error("Euler edge split left an unexpected non-quad face");
    return true;
}

inline void incident_faces(const Mesh& m, Vertex v, std::set<Face>& dst) {
    const Halfedge h = halfedge(v, m);
    if (h == Mesh::null_halfedge()) return;
    for (const Face f : CGAL::faces_around_target(h, m)) {
        if (f != Mesh::null_face() && !m.is_removed(f)) dst.insert(f);
    }
}

inline bool collapse_at_small_angle(Mesh& m, Face f, const Tree& reference,
                                    const Options& opt, const std::set<Vertex>& fixed) {
    if (m.is_removed(f) || angles_of(m, f).min >= opt.min_angle) return false;
    Halfedge shortest = halfedge(f, m);
    double min_len2 = (std::numeric_limits<double>::max)();
    for (const Halfedge h : CGAL::halfedges_around_face(halfedge(f, m), m)) {
        const double d2 = norm2(m.point(source(h, m)) - m.point(target(h, m)));
        if (d2 < min_len2) { min_len2 = d2; shortest = h; }
    }
    const Edge e = edge(shortest, m);
    if (is_border(e, m) || !CGAL::Euler::does_satisfy_link_condition(e, m)) return false;
    const Vertex va = source(shortest, m), vb = target(shortest, m);
    if (is_protected(fixed, va) || is_protected(fixed, vb)) return false;
    std::set<Face> affected;
    incident_faces(m, va, affected);
    incident_faces(m, vb, affected);
    const Point p0 = m.point(va), p1 = m.point(vb);
    const Point midpoint = CGAL::midpoint(p0, p1);
    const Point p = opt.project ? reference.closest_point(midpoint) : midpoint;
    const Face removed0 = face(shortest, m);
    const Face removed1 = face(opposite(shortest, m), m);
    double min_before = 180.0, min_after = 180.0;
    double max_before = 0.0, max_after = 0.0;
    bool any_remaining = false;
    for (const Face q : affected) {
        const auto vv = vertices_of(m, q);
        const Point a = m.point(vv[0]), b = m.point(vv[1]), c = m.point(vv[2]);
        const Angles old = triangle_angles(a, b, c);
        min_before = (std::min)(min_before, old.min);
        max_before = (std::max)(max_before, old.max);
        if (q == removed0 || q == removed1) continue;
        const Point na = (vv[0] == va || vv[0] == vb) ? p : a;
        const Point nb = (vv[1] == va || vv[1] == vb) ? p : b;
        const Point nc = (vv[2] == va || vv[2] == vb) ? p : c;
        if (!acceptable_triangle(na, nb, nc, normal(a, b, c))) return false;
        const Angles newer = triangle_angles(na, nb, nc);
        min_after = (std::min)(min_after, newer.min);
        max_after = (std::max)(max_after, newer.max);
        any_remaining = true;
    }
    if (!any_remaining || min_after <= min_before + 1e-7 ||
        max_after > (std::max)(max_before, opt.max_angle + 5.0)) return false;
    const Vertex kept = CGAL::Euler::collapse_edge(e, m);
    m.point(kept) = p;
    return true;
}

inline double valence_cost(unsigned valence, bool boundary) {
    const double ideal = boundary ? 4.0 : 6.0;
    const double d = static_cast<double>(valence) - ideal;
    return d * d;
}

inline bool flip_if_better(Mesh& m, Edge e, const std::set<Vertex>& fixed) {
    if (m.is_removed(e) || is_border(e, m)) return false;
    const Halfedge h = halfedge(e, m), ho = opposite(h, m);
    const Vertex a = source(h, m), b = target(h, m);
    const Vertex c = target(next(h, m), m);
    const Vertex d = target(next(ho, m), m);
    if (a == d || b == c || c == d ||
        is_protected(fixed, a) || is_protected(fixed, b) ||
        is_protected(fixed, c) || is_protected(fixed, d)) return false;
    if (m.halfedge(c, d).second) return false;
    if (m.degree(a) < 4 || m.degree(b) < 4) return false;
    const double old_energy =
        valence_cost(static_cast<unsigned>(m.degree(a)), false) +
        valence_cost(static_cast<unsigned>(m.degree(b)), false) +
        valence_cost(static_cast<unsigned>(m.degree(c)), false) +
        valence_cost(static_cast<unsigned>(m.degree(d)), false);
    const double new_energy =
        valence_cost(static_cast<unsigned>(m.degree(a)-1), false) +
        valence_cost(static_cast<unsigned>(m.degree(b)-1), false) +
        valence_cost(static_cast<unsigned>(m.degree(c)+1), false) +
        valence_cost(static_cast<unsigned>(m.degree(d)+1), false);
    if (new_energy >= old_energy) return false;
    const Point pa = m.point(a), pb = m.point(b), pc = m.point(c), pd = m.point(d);
    const Vector n0 = normal(pa, pb, pc), n1 = normal(pb, pa, pd);
    if (!acceptable_triangle(pc, pd, pb, n0) ||
        !acceptable_triangle(pd, pc, pa, n1)) return false;
    const double old_min = (std::min)(triangle_angles(pa, pb, pc).min,
                                       triangle_angles(pb, pa, pd).min);
    const double new_min = (std::min)(triangle_angles(pc, pd, pb).min,
                                       triangle_angles(pd, pc, pa).min);
    if (new_min < (std::min)(old_min, 10.0)) return false;
    CGAL::Euler::flip_edge(h, m);
    return true;
}

inline bool vertex_move_valid(const Mesh& m, Vertex v, const Point& candidate) {
    std::set<Face> adjacent;
    incident_faces(m, v, adjacent);
    for (const Face f : adjacent) {
        const auto vv = vertices_of(m, f);
        const Point a = m.point(vv[0]), b = m.point(vv[1]), c = m.point(vv[2]);
        if (!acceptable_triangle(vv[0] == v ? candidate : a,
                                 vv[1] == v ? candidate : b,
                                 vv[2] == v ? candidate : c,
                                 normal(a, b, c))) return false;
    }
    return true;
}

inline void tangential_smooth(Mesh& m, const Tree& reference,
                              const Options& opt, const std::set<Vertex>& fixed) {
    for (unsigned step = 0; step < opt.smoothing_steps; ++step) {
        std::vector<std::pair<Vertex, Point>> proposed;
        for (const Vertex v : vertices(m)) {
            if (is_protected(fixed, v) || m.is_border(v)) continue;
            const Halfedge h = halfedge(v, m);
            if (h == Mesh::null_halfedge()) continue;
            double x = 0, y = 0, z = 0;
            std::size_t count = 0;
            for (const Vertex neighbor : CGAL::vertices_around_target(h, m)) {
                const Point& p = m.point(neighbor);
                x += p.x(); y += p.y(); z += p.z(); ++count;
            }
            if (count < 3) continue;
            const Point old = m.point(v);
            const Vector lap = Point(x/count, y/count, z/count) - old;
            Vector n = CGAL::NULL_VECTOR;
            std::set<Face> adjacent;
            incident_faces(m, v, adjacent);
            for (const Face f : adjacent) {
                const auto t = vertices_of(m, f);
                n = n + normal(m.point(t[0]), m.point(t[1]), m.point(t[2]));
            }
            const double nn = norm2(n);
            if (nn < 1e-28) continue;
            const Vector tangent = lap - n * (dot(lap, n) / nn);
            const Point trial = old + tangent * 0.5;
            proposed.emplace_back(v, opt.project ? reference.closest_point(trial) : trial);
        }
        for (const auto& move : proposed) {
            if (vertex_move_valid(m, move.first, move.second))
                m.point(move.first) = move.second;
        }
    }
}

inline Statistics analyze(const Mesh& m, const Options& opt) {
    Statistics s;
    s.vertices_after = num_vertices(m);
    s.faces_after = num_faces(m);
    for (const Face f : faces(m)) {
        const Angles a = angles_of(m, f);
        s.angles.min = (std::min)(s.angles.min, a.min);
        s.angles.max = (std::max)(s.angles.max, a.max);
        if (a.min < opt.min_angle) ++s.min_angle_violations;
        if (a.max > opt.max_angle) ++s.max_angle_violations;
    }
    return s;
}

inline Statistics run(Mesh& mesh, const Options& opt) {
    if (!CGAL::is_triangle_mesh(mesh) || num_faces(mesh) == 0)
        throw std::invalid_argument("Wang2019 requires a nonempty triangle mesh");
    if (!(opt.min_angle > 0.0 && opt.min_angle < 60.0 &&
          opt.max_angle > 60.0 && opt.max_angle < 180.0 &&
          opt.budget_fraction > 0.0 && opt.budget_fraction <= 1.0 &&
          opt.feature_angle > 0.0 && opt.feature_angle < 180.0))
        throw std::invalid_argument("Invalid angle, feature, or operation-budget bounds");

    // Projection must use the original input, not the dynamically edited mesh.
    const Mesh original(mesh);
    Tree reference(faces(original).first, faces(original).second, original);
    reference.accelerate_distance_queries();

    Statistics total;
    total.vertices_before = num_vertices(mesh);
    total.faces_before = num_faces(mesh);
    // Recomputed each round; new vertices are unrestricted, while original
    // boundary and sharp vertices remain frozen (conservative approximation).
    for (unsigned round = 0; round < opt.rounds; ++round) {
        std::set<Vertex> fixed = protected_vertices(mesh,
            opt.protect_features ? opt.feature_angle : 179.999);
        if (!opt.protect_features) {
            fixed.clear();
            for (const Vertex v : vertices(mesh))
                if (mesh.is_border(v)) fixed.insert(v);
        }
        const std::size_t budget = (std::max)(std::size_t(1),
            static_cast<std::size_t>(std::ceil(opt.budget_fraction * num_vertices(mesh))));
        std::size_t added = 0, removed = 0, flipped = 0;

        std::vector<Face> descending;
        for (const Face f : faces(mesh))
            if (angles_of(mesh, f).max > opt.max_angle) descending.push_back(f);
        std::sort(descending.begin(), descending.end(), [&](Face a, Face b) {
            return angles_of(mesh, a).max > angles_of(mesh, b).max;
        });
        for (const Face f : descending) {
            if (added >= budget) break;
            if (insert_at_large_angle(mesh, f, reference, opt, fixed)) ++added;
        }

        // Refresh sharp-feature classification after insertions.
        fixed = protected_vertices(mesh,
            opt.protect_features ? opt.feature_angle : 179.999);
        if (!opt.protect_features) {
            fixed.clear();
            for (const Vertex v : vertices(mesh))
                if (mesh.is_border(v)) fixed.insert(v);
        }

        std::vector<Face> ascending;
        for (const Face f : faces(mesh))
            if (angles_of(mesh, f).min < opt.min_angle) ascending.push_back(f);
        std::sort(ascending.begin(), ascending.end(), [&](Face a, Face b) {
            return angles_of(mesh, a).min < angles_of(mesh, b).min;
        });
        for (const Face f : ascending) {
            if (removed >= budget) break;
            if (collapse_at_small_angle(mesh, f, reference, opt, fixed)) ++removed;
        }

        fixed = protected_vertices(mesh,
            opt.protect_features ? opt.feature_angle : 179.999);
        if (!opt.protect_features) {
            fixed.clear();
            for (const Vertex v : vertices(mesh))
                if (mesh.is_border(v)) fixed.insert(v);
        }
        std::vector<Edge> edge_snapshot;
        for (const Edge e : edges(mesh)) edge_snapshot.push_back(e);
        for (const Edge e : edge_snapshot)
            if (flip_if_better(mesh, e, fixed)) ++flipped;

        tangential_smooth(mesh, reference, opt, fixed);
        total.insertions += added;
        total.collapses += removed;
        total.flips += flipped;
        if (opt.verbose) {
            const Statistics a = analyze(mesh, opt);
            std::cout << "[Wang2019] round=" << (round+1)
                      << " insert=" << added << " collapse=" << removed
                      << " flip=" << flipped << " V=" << a.vertices_after
                      << " min_angle=" << a.angles.min
                      << " max_angle=" << a.angles.max
                      << " below=" << a.min_angle_violations
                      << " above=" << a.max_angle_violations << '\n';
        }
        if (added == 0 && removed == 0 && flipped == 0) break;
    }
    const Statistics last = analyze(mesh, opt);
    total.vertices_after = last.vertices_after;
    total.faces_after = last.faces_after;
    total.min_angle_violations = last.min_angle_violations;
    total.max_angle_violations = last.max_angle_violations;
    total.angles = last.angles;
    return total;
}
} // namespace wang2019
} // namespace rar
