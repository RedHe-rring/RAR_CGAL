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
#include <string>
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
    double budget_fraction = 0.2;   // k as fraction of bad triangles (Sec. 5 timing setting)
    unsigned k = 0;                 // optional explicit per-stage attempt budget
    bool strict_vertex_count = true; // roll back a round if split/collapse cannot balance
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
    std::size_t rolled_back_rounds = 0;
    std::size_t vertices_before = 0;
    std::size_t vertices_after = 0;
    std::size_t faces_before = 0;
    std::size_t faces_after = 0;
    std::size_t min_angle_violations = 0;
    std::size_t max_angle_violations = 0;
    Angles angles;
};

enum class RejectCause : std::size_t {
    stale = 0, boundary, feature, geometry, no_improvement, topology, count
};
struct RejectionStats {
    std::array<std::size_t, 7> rejected{};
    void record(RejectCause cause) { ++rejected[static_cast<std::size_t>(cause)]; }
};
inline const char* cause_name(std::size_t i) {
    const char* names[] = {"stale", "boundary", "feature", "geometry",
                           "no-improvement", "topology", "count"};
    return names[i];
}
inline void verify_delta(const Mesh& m,
                         std::size_t old_vertices, std::size_t old_faces,
                         int expected_vertices, int expected_faces,
                         const char* operation) {
    const auto expected_v = static_cast<std::ptrdiff_t>(old_vertices) + expected_vertices;
    const auto expected_f = static_cast<std::ptrdiff_t>(old_faces) + expected_faces;
    if (static_cast<std::ptrdiff_t>(m.number_of_vertices()) != expected_v ||
        static_cast<std::ptrdiff_t>(m.number_of_faces()) != expected_f) {
        throw std::runtime_error(std::string("Wang2019 ") + operation +
            " V=" + std::to_string(m.number_of_vertices()) +
            " expected=" + std::to_string(expected_v) +
            " F=" + std::to_string(m.number_of_faces()) +
            " expected=" + std::to_string(expected_f));
    }
}

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

// With a fixed edge midpoint, splitting an existing input triangle does not
// change its piecewise-linear geometry. Reprojecting the split vertex could.
inline bool is_feature_edge(const Mesh& m, Edge e, const Options& opt) {
    const Halfedge h = halfedge(e, m);
    if (is_border(h, m) || is_border(opposite(h, m), m)) return true;
    if (!opt.protect_features) return false;
    const auto v0 = vertices_of(m, face(h, m));
    const auto v1 = vertices_of(m, face(opposite(h, m), m));
    const Vector n0 = normal(m.point(v0[0]), m.point(v0[1]), m.point(v0[2]));
    const Vector n1 = normal(m.point(v1[0]), m.point(v1[1]), m.point(v1[2]));
    const double size = std::sqrt(norm2(n0) * norm2(n1));
    if (size <= 1e-30) return true;
    const double threshold = std::cos(opt.feature_angle * 3.14159265358979323846 / 180.0);
    return dot(n0, n1) / size < threshold;
}

inline double angle_energy(const Point& a, const Point& b, const Point& c) {
    // Squared angle deviation from the equilateral configuration (60 degrees).
    const double a0 = angle(b, a, c)-60.0;
    const double a1 = angle(a, b, c)-60.0;
    const double a2 = angle(a, c, b)-60.0;
    return a0*a0 + a1*a1 + a2*a2;
}

// Evaluate (without mutating) flipping this triangle-triangle edge.
// Used to select a local angle-optimal flip after the center-edge split.
// This does not yet enumerate the feature-specific Fig. 4(c-f) cases.
inline bool evaluate_angle_flip(const Mesh& m, Edge e, const Options& opt,
                                double& score) {
    if (m.is_removed(e) || is_feature_edge(m, e, opt)) return false;
    const Halfedge h = halfedge(e, m), ho = opposite(h, m);
    const Vertex a = source(h, m), b = target(h, m);
    const Vertex c = target(next(h, m), m), d = target(next(ho, m), m);
    if (c == d || a == d || b == c ||
        m.halfedge(c, d) != Mesh::null_halfedge() ||
        m.degree(a) < 4 || m.degree(b) < 4) return false;
    const Point pa=m.point(a), pb=m.point(b), pc=m.point(c), pd=m.point(d);
    const Vector n0=normal(pa,pb,pc), n1=normal(pb,pa,pd);
    if (!acceptable_triangle(pc,pd,pb,n0) ||
        !acceptable_triangle(pd,pc,pa,n1)) return false;
    const double old_max=(std::max)(triangle_angles(pa,pb,pc).max,
                                      triangle_angles(pb,pa,pd).max);
    const double new_max=(std::max)(triangle_angles(pc,pd,pb).max,
                                      triangle_angles(pd,pc,pa).max);
    if (new_max > (std::max)(old_max, opt.max_angle) + 1e-7) return false;
    score=angle_energy(pc,pd,pb)+angle_energy(pd,pc,pa);
    return true;
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

inline double large_angle_violation(const Point& a, const Point& b,
                                    const Point& c, double upper_bound) {
    const double error=(std::max)(0.0,
                                  triangle_angles(a,b,c).max-upper_bound);
    return error*error;
}

inline bool insert_at_large_angle(Mesh& m, Face f, const Tree&,
                                  const Options& opt, const std::set<Vertex>&,
                                  RejectionStats* reasons = nullptr,
                                  std::size_t* successful_angle_flips = nullptr) {
    auto reject = [&](RejectCause cause) {
        if (reasons) reasons->record(cause);
        return false;
    };
    if (m.is_removed(f) || angles_of(m,f).max<=opt.max_angle)
        return reject(RejectCause::stale);

    Halfedge chosen=halfedge(f,m);
    double longest=-1.0;
    for (const Halfedge h : CGAL::halfedges_around_face(halfedge(f,m),m)) {
        const double len2=norm2(m.point(source(h,m))-m.point(target(h,m)));
        if (len2>longest) { longest=len2; chosen=h; }
    }
    const Halfedge ho=opposite(chosen,m);
    if (is_border(chosen,m) || is_border(ho,m))
        return reject(RejectCause::boundary);
    if (is_feature_edge(m,edge(chosen,m),opt))
        return reject(RejectCause::feature);
    const Vertex va=source(chosen,m),vb=target(chosen,m),
                 vc=target(next(chosen,m),m),vd=target(next(ho,m),m);
    const Point a=m.point(va),b=m.point(vb),
                c=m.point(vc),d=m.point(vd);
    const Point p=CGAL::midpoint(a,b); // exact PL split; do not project
    const Vector n0=normal(a,b,c),n1=normal(b,a,d);
    if (!acceptable_triangle(a,p,c,n0) ||
        !acceptable_triangle(p,b,c,n0) ||
        !acceptable_triangle(b,p,d,n1) ||
        !acceptable_triangle(p,a,d,n1))
        return reject(RejectCause::geometry);

    const std::array<std::array<Point,3>,4> created{{
        {{a,p,c}},{{p,b,c}},{{b,p,d}},{{p,a,d}}
    }};
    // Each pair (u,v) is the oriented boundary edge of the corresponding
    // split triangle, i.e. that triangle has cyclic order (u,v,p).
    const std::array<std::pair<Vertex,Vertex>,4> outer{{
        {vc,va},{vb,vc},{vd,vb},{va,vd}
    }};

    int best_index=-1;
    double best_score=(std::numeric_limits<double>::infinity)();
    for (int i=0;i<4;++i) {
        const Vertex u=outer[i].first,v=outer[i].second;
        const Halfedge h=m.halfedge(u,v);
        if (h==Mesh::null_halfedge()) continue;
        const Edge e=edge(h,m);
        if (is_feature_edge(m,e,opt)) continue;
        // After the shared-edge split, c and d gain a new connection to p.
        const std::size_t du=m.degree(u)+(u==vc||u==vd?1:0);
        const std::size_t dv=m.degree(v)+(v==vc||v==vd?1:0);
        if (du<4 || dv<4) continue;
        Face outside=face(h,m);
        if (outside==face(chosen,m) || outside==face(ho,m))
            outside=face(opposite(h,m),m);
        if (outside==Mesh::null_face()) continue;
        const auto ov=vertices_of(m,outside);
        Vertex x=Mesh::null_vertex();
        for (const Vertex t : ov)
            if (t!=u && t!=v) x=t;
        if (x==Mesh::null_vertex() || x==va || x==vb) continue;
        const Point pu=m.point(u),pv=m.point(v),px=m.point(x);
        const Vector old_inner=normal(pu,pv,p);
        const Vector old_outer=normal(pv,pu,px);
        if (!acceptable_triangle(p,px,pv,old_inner) ||
            !acceptable_triangle(px,p,pu,old_outer)) continue;

        // Evaluate the 5-face patch (4 split faces plus one outside face)
        // AFTER the candidate flip, before mutating the input mesh.
        double after_violation=0.0, after_energy=0.0;
        for (int j=0;j<4;++j) if (j!=i) {
            const auto& t=created[j];
            after_violation+=large_angle_violation(t[0],t[1],t[2],opt.max_angle);
            after_energy+=angle_energy(t[0],t[1],t[2]);
        }
        after_violation+=large_angle_violation(p,px,pv,opt.max_angle);
        after_violation+=large_angle_violation(px,p,pu,opt.max_angle);
        after_energy+=angle_energy(p,px,pv)+angle_energy(px,p,pu);
        const double before_violation=
            large_angle_violation(a,b,c,opt.max_angle)+
            large_angle_violation(b,a,d,opt.max_angle)+
            large_angle_violation(m.point(ov[0]),m.point(ov[1]),
                                  m.point(ov[2]),opt.max_angle);
        if (after_violation + 1e-7 < before_violation &&
            after_energy < best_score) {
            best_score=after_energy;
            best_index=i;
        }
    }

    const double original_max=triangle_angles(a,b,c).max;
    const double split_target_max=(std::max)(
        triangle_angles(a,p,c).max,triangle_angles(p,b,c).max);
    const bool split_only_improves=split_target_max<original_max-1e-7;
    if (best_index<0 && !split_only_improves)
        return reject(RejectCause::no_improvement);

    const auto oldV=m.number_of_vertices(),oldF=m.number_of_faces();
    const Face f0=face(chosen,m),f1=face(ho,m);
    const Halfedge new_h=CGAL::Euler::split_edge(chosen,m);
    const Vertex vm=target(new_h,m);
    m.point(vm)=p;
    if (!split_quad(m,f0,vm) || !split_quad(m,f1,vm))
        throw std::runtime_error("Split did not generate two triangular pairs");
    verify_delta(m,oldV,oldF,1,2,"split");

    if (best_index>=0) {
        const auto uv=outer[best_index];
        const Halfedge h=m.halfedge(uv.first,uv.second);
        if (h==Mesh::null_halfedge())
            throw std::runtime_error("Chosen preflight flip edge vanished");
        CGAL::Euler::flip_edge(h,m);
        if (successful_angle_flips) ++*successful_angle_flips;
        verify_delta(m,oldV,oldF,1,2,"split+angle-flip");
    }
    // The simple split fallback is allowed only when the targeted obtuse
    // triangle improves. Full Fig. 4(c-f) feature patterns are still missing.
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
                                    const Options& opt, const std::set<Vertex>& fixed,
                                    RejectionStats* reasons = nullptr) {
    auto reject = [&](RejectCause cause) {
        if (reasons) reasons->record(cause);
        return false;
    };
    if (m.is_removed(f) || angles_of(m, f).min >= opt.min_angle)
        return reject(RejectCause::stale);
    Halfedge shortest = halfedge(f, m);
    double min_len2 = (std::numeric_limits<double>::max)();
    for (const Halfedge h : CGAL::halfedges_around_face(halfedge(f, m), m)) {
        const double d2 = norm2(m.point(source(h, m)) - m.point(target(h, m)));
        if (d2 < min_len2) { min_len2 = d2; shortest = h; }
    }
    const Edge e = edge(shortest, m);
    if (is_border(e, m)) return reject(RejectCause::boundary);
    if (!CGAL::Euler::does_satisfy_link_condition(e, m)) return reject(RejectCause::topology);
    const Vertex va = source(shortest, m), vb = target(shortest, m);
    if (is_protected(fixed, va) || is_protected(fixed, vb)) return reject(RejectCause::feature);
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
        if (!acceptable_triangle(na, nb, nc, normal(a, b, c))) return reject(RejectCause::geometry);
        const Angles newer = triangle_angles(na, nb, nc);
        min_after = (std::min)(min_after, newer.min);
        max_after = (std::max)(max_after, newer.max);
        any_remaining = true;
    }
    if (!any_remaining || min_after <= min_before + 1e-7 ||
        max_after > (std::max)(max_before, opt.max_angle + 5.0)) return reject(RejectCause::no_improvement);
    const auto old_vertices = m.number_of_vertices(), old_faces = m.number_of_faces();
    const Vertex kept = CGAL::Euler::collapse_edge(e, m);
    m.point(kept) = p;
    verify_delta(m, old_vertices, old_faces, -1, -2, "collapse");
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
    if (m.halfedge(c, d) != Mesh::null_halfedge()) return false;
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
            std::set<Face> adjacent;
            incident_faces(m, v, adjacent);
            if (adjacent.size() < 3) continue;

            // Eq. (1): neighboring triangle centroids, not the mean of
            // adjacent vertex coordinates. Triangle area is used as w_j here;
            // the paper does not define the precise weights in Sec. 4.2.3.
            double sum_x=0, sum_y=0, sum_z=0, weight_sum=0;
            Vector vertex_normal=CGAL::NULL_VECTOR;
            for (const Face f : adjacent) {
                const auto tri=vertices_of(m,f);
                const Point pa=m.point(tri[0]),pb=m.point(tri[1]),pc=m.point(tri[2]);
                const Vector normal_vector=normal(pa,pb,pc);
                const double twice_area=std::sqrt(norm2(normal_vector));
                if (twice_area <= 1e-16) continue;
                sum_x+=(pa.x()+pb.x()+pc.x())/3.0*twice_area;
                sum_y+=(pa.y()+pb.y()+pc.y())/3.0*twice_area;
                sum_z+=(pa.z()+pb.z()+pc.z())/3.0*twice_area;
                weight_sum+=twice_area;
                vertex_normal=vertex_normal+normal_vector;
            }
            const double nn=norm2(vertex_normal);
            if (weight_sum <= 0 || nn < 1e-28) continue;
            const Point old=m.point(v);
            const Vector direction=Point(sum_x/weight_sum,
                                         sum_y/weight_sum,
                                         sum_z/weight_sum)-old;
            const Vector tangent=direction-
                vertex_normal*(dot(direction,vertex_normal)/nn);
            const Point trial=old+tangent;
            proposed.emplace_back(v,opt.project?reference.closest_point(trial):trial);
        }
        for (const auto& move : proposed)
            if (vertex_move_valid(m,move.first,move.second))
                m.point(move.first)=move.second;
    }
}

inline Statistics analyze(const Mesh& m, const Options& opt) {
    Statistics s;
    s.vertices_after = m.number_of_vertices();
    s.faces_after = m.number_of_faces();
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
    if (!CGAL::is_triangle_mesh(mesh) || mesh.number_of_faces()==0 || !mesh.is_valid())
        throw std::invalid_argument("Wang2019 requires a valid, nonempty triangle mesh");
    if (!(opt.min_angle>0.0 && opt.min_angle<60.0 &&
          opt.max_angle>60.0 && opt.max_angle<180.0 &&
          opt.budget_fraction>0.0 && opt.budget_fraction<=1.0 &&
          opt.feature_angle>0.0 && opt.feature_angle<180.0))
        throw std::invalid_argument("Invalid angle or operation-budget bounds");

    const Mesh original(mesh);
    Tree reference(faces(original).first,faces(original).second,original);
    reference.accelerate_distance_queries();

    Statistics total;
    total.vertices_before=mesh.number_of_vertices();
    total.faces_before=mesh.number_of_faces();

    const auto feature_vertices = [&]() {
        std::set<Vertex> fixed;
        if (opt.protect_features)
            return protected_vertices(mesh,opt.feature_angle);
        for (const Vertex v : vertices(mesh))
            if (mesh.is_border(v)) fixed.insert(v);
        return fixed;
    };
    const auto valence_stage = [&]() {
        std::size_t flips=0;
        const std::set<Vertex> fixed=feature_vertices();
        std::vector<Edge> candidates;
        for (const Edge e : edges(mesh)) candidates.push_back(e);
        for (const Edge e : candidates)
            if (flip_if_better(mesh,e,fixed)) ++flips;
        tangential_smooth(mesh,reference,opt,feature_vertices());
        return flips;
    };

    for (unsigned round=0;round<opt.rounds;++round) {
        const Statistics initial=analyze(mesh,opt);
        if (!initial.min_angle_violations && !initial.max_angle_violations)
            break;

        const std::size_t requested_k=opt.k ?
            static_cast<std::size_t>(opt.k) :
            (std::max)(std::size_t(1),
                static_cast<std::size_t>(std::ceil(
                    opt.budget_fraction*(std::max)(
                        initial.min_angle_violations,initial.max_angle_violations))));
        // A stage may attempt at most the number of present large-angle faces.
        const std::size_t base_k=(std::min)(requested_k, initial.max_angle_violations);
        if (base_k==0) break;

        // The paper aims for k insertions and k deletions without changing N.
        // For this experimental backend, if the chosen batch cannot be
        // balanced, reduce k and retry the SAME round. This retains useful
        // smaller updates instead of always rolling back an entire iteration.
        const Mesh round_start=opt.strict_vertex_count ? Mesh(mesh) : Mesh();
        std::size_t trial_k=base_k;
        bool committed=false;
        std::size_t inserted=0, collapsed=0, flipped=0, used_k=0;
        RejectionStats final_split_rejections, final_collapse_rejections;
        double effective_min_angle=opt.min_angle;

        while (true) {
            if (opt.strict_vertex_count)
                mesh=round_start;
            const std::size_t vertices_before=mesh.number_of_vertices();
            const std::size_t faces_before=mesh.number_of_faces();
            RejectionStats split_rejections, collapse_rejections;
            std::size_t added=0, removed=0, angle_flips=0;

            std::vector<Face> large;
            for (const Face f : faces(mesh))
                if (angles_of(mesh,f).max>opt.max_angle) large.push_back(f);
            std::sort(large.begin(),large.end(),[&](Face x,Face y) {
                return angles_of(mesh,x).max>angles_of(mesh,y).max;
            });
            std::set<Vertex> fixed=feature_vertices();
            for (const Face f : large) {
                if (added>=trial_k) break;
                if (insert_at_large_angle(mesh,f,reference,opt,fixed,
                                          &split_rejections,&angle_flips))
                    ++added;
            }
            const std::size_t first_flips=valence_stage();

            // Under strict N, compensate only the insertions that actually
            // succeeded, not the requested attempts.
            const std::size_t target_removals=opt.strict_vertex_count ?
                                               added : trial_k;

            std::vector<Face> small;
            for (const Face f : faces(mesh)) small.push_back(f);
            std::sort(small.begin(),small.end(),[&](Face x,Face y) {
                return angles_of(mesh,x).min<angles_of(mesh,y).min;
            });
            // Sec. 4.2.4: temporary relaxation only if the candidate count
            // below the specified beta_min is insufficient for this batch.
            std::size_t below_user_bound=0;
            for (const Face f : small)
                if (angles_of(mesh,f).min<opt.min_angle) ++below_user_bound;
            double raised_beta=opt.min_angle;
            if (below_user_bound<target_removals && target_removals>0 &&
                !small.empty()) {
                const std::size_t last=(std::min)(target_removals,small.size())-1;
                raised_beta=std::nextafter(
                    (std::max)(opt.min_angle,angles_of(mesh,small[last]).min),
                    (std::numeric_limits<double>::infinity)());
            }
            fixed=feature_vertices();
            Options acceptance=opt;
            acceptance.min_angle=raised_beta;
            for (const Face f : small) {
                if (removed>=target_removals) break;
                if (collapse_at_small_angle(mesh,f,reference,acceptance,fixed,
                                            &collapse_rejections))
                    ++removed;
            }
            const std::size_t second_flips=valence_stage();
            if (!mesh.is_valid() || !CGAL::is_triangle_mesh(mesh))
                throw std::runtime_error("Wang2019: invalid mesh after iteration");

            const auto actual_dv=static_cast<std::ptrdiff_t>(mesh.number_of_vertices())-
                                 static_cast<std::ptrdiff_t>(vertices_before);
            const auto actual_df=static_cast<std::ptrdiff_t>(mesh.number_of_faces())-
                                 static_cast<std::ptrdiff_t>(faces_before);
            const auto expected_dv=static_cast<std::ptrdiff_t>(added)-
                                   static_cast<std::ptrdiff_t>(removed);
            if (actual_dv!=expected_dv || actual_df!=2*expected_dv)
                throw std::runtime_error("Wang2019: iteration V/F delta mismatch");

            if (!opt.strict_vertex_count || added==removed) {
                committed=true;
                inserted=added;
                collapsed=removed;
                flipped=first_flips+second_flips+angle_flips;
                used_k=trial_k;
                effective_min_angle=raised_beta;
                final_split_rejections=split_rejections;
                final_collapse_rejections=collapse_rejections;
                break;
            }

            ++total.rolled_back_rounds;
            if (opt.verbose)
                std::cout << "[Wang2019] round=" << round+1
                          << " retry: k=" << trial_k
                          << " inserted=" << added
                          << " collapsed=" << removed << '\n';
            if (trial_k==1) {
                mesh=round_start;
                final_split_rejections=split_rejections;
                final_collapse_rejections=collapse_rejections;
                break;
            }
            trial_k=(std::max)(std::size_t(1),trial_k/2);
        }
        if (!committed) {
            if (opt.verbose)
                std::cout << "[Wang2019] round=" << round+1
                          << " stopped: no balanced operation batch\n";
            break;
        }

        total.insertions+=inserted;
        total.collapses+=collapsed;
        total.flips+=flipped;
        const Statistics current=analyze(mesh,opt);
        if (opt.verbose) {
            std::cout << "[Wang2019] round=" << round+1
                      << " k=" << used_k
                      << " insert=" << inserted
                      << " collapse=" << collapsed
                      << " flip=" << flipped
                      << " V=" << current.vertices_after
                      << " min_angle=" << current.angles.min
                      << " max_angle=" << current.angles.max
                      << " below=" << current.min_angle_violations
                      << " above=" << current.max_angle_violations;
            if (effective_min_angle>opt.min_angle)
                std::cout << " temporary_min_angle=" << effective_min_angle;
            std::cout << '\n';
            for (std::size_t i=0;i<final_split_rejections.rejected.size();++i)
                if (final_split_rejections.rejected[i] ||
                    final_collapse_rejections.rejected[i])
                    std::cout << "[Wang2019] rejected_" << cause_name(i)
                              << " split=" << final_split_rejections.rejected[i]
                              << " collapse=" << final_collapse_rejections.rejected[i] << '\n';
        }
        if (!current.min_angle_violations && !current.max_angle_violations)
            break;
        if (!inserted && !collapsed && !flipped &&
            current.min_angle_violations>=initial.min_angle_violations &&
            current.max_angle_violations>=initial.max_angle_violations)
            break;
    }
    const Statistics last=analyze(mesh,opt);
    total.vertices_after=last.vertices_after;
    total.faces_after=last.faces_after;
    total.min_angle_violations=last.min_angle_violations;
    total.max_angle_violations=last.max_angle_violations;
    total.angles=last.angles;
    if (opt.strict_vertex_count && total.vertices_after!=total.vertices_before)
        throw std::runtime_error("Wang2019 strict-N invariant violated");
    return total;
}

} // namespace wang2019
} // namespace rar
