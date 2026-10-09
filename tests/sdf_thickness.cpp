#include "rar/SDFThicknessField.h"
#include "rar/OutputNaming.h"

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/iterator.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message); // active in Release
}
bool close(double a, double b, double tol = 1e-9)
{
    return std::abs(a-b) <= tol * (std::max)({1.0, std::abs(a), std::abs(b)});
}
template<class Fn> void rejects(Fn fn)
{
    bool rejected = false;
    try { fn(); } catch (const std::exception&) { rejected = true; }
    require(rejected, "invalid input was not rejected");
}

void add_box(rar::Mesh& mesh, double shift = 0.0)
{
    using rar::Point;
    std::vector<rar::Mesh::Vertex_index> v;
    for (const auto& p : {Point(0,0,0), Point(4,0,0), Point(4,4,0), Point(0,4,0),
                          Point(0,0,0.2), Point(4,0,0.2), Point(4,4,0.2), Point(0,4,0.2)})
        v.push_back(mesh.add_vertex(Point(p.x()+shift, p.y()+shift, p.z()+shift/10.0)));
    const int faces[][3] = {{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},
                            {1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    for (const auto& f : faces)
        require(mesh.add_face(v[f[0]],v[f[1]],v[f[2]]) != rar::Mesh::null_face(), "box face");
}

void formula_tests()
{
    using rar::sdf_thickness_target_length;
    const auto user_eps = sdf_thickness_target_length(2.0, 10.0, .01, .001, 10.0);
    require(user_eps.epsilon_limited && !user_eps.radius_limited && !user_eps.thickness_limited,
            "epsilon branch");
    require(close(user_eps.target_length, std::sqrt(.0297)), "epsilon branch value");
    const auto radius = sdf_thickness_target_length(2.0, 10.0, .5, .001, 10.0);
    require(radius.radius_limited && close(radius.epsilon_local, .0915), "radius cap");
    require(close(radius.target_length, std::sqrt(6*.0915/2 - 3*.0915*.0915)), "radius formula");
    const auto thin = sdf_thickness_target_length(2.0, .02, .5, .001, 10.0);
    require(thin.thickness_limited && close(thin.epsilon_local, .005), "thickness cap");
    require(thin.target_length < radius.target_length, "thickness must reduce target in safe formula range");
    const auto tie = sdf_thickness_target_length(2.0, .02, .005, .001, 10.0);
    require(tie.epsilon_limited && tie.thickness_limited, "limiter ties");
    require(sdf_thickness_target_length(2,.02,.5,.3,10).min_edge_clipped, "lower clipping");
    require(close(sdf_thickness_target_length(2,.02,.5,.3,10).target_length,.3), "lower bound");
    require(sdf_thickness_target_length(2,.02,.5,.001,.1).max_edge_clipped, "upper clipping");
    require(close(sdf_thickness_target_length(2,.02,.5,.001,.1).target_length,.1), "upper bound");
    const auto flat = sdf_thickness_target_length(0,.02,.5,.001,10);
    require(flat.target_length == 10 && flat.max_edge_clipped && flat.thickness_limited, "flat limit");
    for (double scale : {1e-100, 7.0, 1e100}) {
        const auto scaled = sdf_thickness_target_length(2/scale,.02*scale,.5*scale,.001*scale,10*scale);
        require(close(scaled.target_length/scale, thin.target_length), "formula scale covariance");
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double k : {-1.0, nan, inf}) rejects([&] { sdf_thickness_target_length(k,1,.1,.01,1); });
    for (double t : {-1.0, 0.0, nan, inf}) rejects([&] { sdf_thickness_target_length(1,t,.1,.01,1); });
    for (double eps : {-1.0, 0.0, nan, inf}) rejects([&] { sdf_thickness_target_length(1,1,eps,.01,1); });
    rejects([&] { sdf_thickness_target_length(1,1,.1,2,1); });
    rejects([&] { sdf_thickness_target_length(1,1,.1,0,1); });
    rejects([&] { sdf_thickness_target_length(1,1,.1,.01,inf); });
}

void field_tests()
{
    rar::Mesh mesh;
    add_box(mesh);
    rar::RemeshConfig cfg;
    cfg.field_type = rar::FieldType::CGALAdaptiveSDFThickness;
    cfg.epsilon = .5;
    cfg.min_edge_length = .0001;
    cfg.max_edge_length = 10;
    cfg.sdf_cone_angle_degrees = .01;
    const auto field = rar::make_sdf_thickness_field(mesh, cfg);
    require(mesh.number_of_vertices()==8 && mesh.number_of_faces()==12, "field changed topology");
    require(field.stats.thickness_limited_vertices > 0, "slab did not activate thickness limiter");
    // Near-normal rays on the broad faces of a closed slab measure its known
    // thickness. This detects normalization, wrong ray direction, or radius/diameter confusion.
    for (const auto f : mesh.faces()) {
        if (f.idx() < 4) require(close(field.face_thickness[f], .2, 1e-7), "raw slab SDF has wrong units");
    }
    for (const auto v : mesh.vertices()) {
        double sum = 0, area = 0;
        for (const auto h : CGAL::halfedges_around_target(v, mesh)) {
            const auto f = mesh.face(h);
            sum += field.face_area[f]*field.face_thickness[f];
            area += field.face_area[f];
        }
        require(close(field.thickness[v], sum/area), "face to vertex area transfer");
        require(close(field.vertex_area[v], area/3), "vertex area weights");
        require(field.samples[v].epsilon_local * field.curvature[v] <= .183 + 1e-14, "radius epsilon cap");
        require(field.samples[v].epsilon_local <= field.thickness[v]/4, "thickness epsilon cap");
        require(std::isfinite(field.target[v]) && field.target[v]>0, "invalid h0");
    }
    rar::Mesh scaled(mesh);
    for (const auto v : scaled.vertices()) {
        const auto p=scaled.point(v);
        scaled.point(v)=rar::Point(7*p.x(),7*p.y(),7*p.z());
    }
    auto scaled_cfg=cfg;
    scaled_cfg.epsilon*=7;
    scaled_cfg.min_edge_length*=7;
    scaled_cfg.max_edge_length*=7;
    const auto scaled_field=rar::make_sdf_thickness_field(scaled,scaled_cfg);
    double max_k_scale_residual=0.0, max_h_scale_residual=0.0;
    for (const auto v : mesh.vertices()) {
        require(close(scaled_field.thickness[v]/7,field.thickness[v]), "SDF physical scale");
        // The unchanged CGAL estimator adds K=1000*average_edge_length*area
        // before its eigensolve. That penalty is not scale covariant. Record
        // its effect separately instead of imposing an unsupported contract
        // on CGAL. Raw SDF and the scalar formula have strict scaling tests.
        max_k_scale_residual=(std::max)(max_k_scale_residual,
            std::abs(scaled_field.curvature[v]*7/field.curvature[v]-1));
        max_h_scale_residual=(std::max)(max_h_scale_residual,
            std::abs(scaled_field.target[v]/(7*field.target[v])-1));
        const double e=scaled_field.samples[v].epsilon_local;
        const double direct=std::sqrt(6*e/scaled_field.curvature[v]-3*e*e);
        require(close(scaled_field.target[v],direct), "integrated h0 disagrees with supplied curvature and SDF");
    }
    std::cout << std::setprecision(17) << "CGAL curvature scale residual: " << max_k_scale_residual
              << "; integrated h0 residual: " << max_h_scale_residual << '\n';
    const auto repeated=rar::make_sdf_thickness_field(mesh,cfg);
    for (const auto v : mesh.vertices())
        require(close(repeated.thickness[v],scaled_field.thickness[v]/7), "repeated field accumulation");
    auto reverse=mesh;
    CGAL::Polygon_mesh_processing::reverse_face_orientations(reverse);
    rejects([&] { rar::make_sdf_thickness_field(reverse,cfg); });
    rar::Mesh crossing;
    add_box(crossing);
    add_box(crossing,.5);
    rejects([&] { rar::validate_sdf_thickness_mesh(crossing); });
    rar::Mesh open;
    const auto a=open.add_vertex({0,0,0}), b=open.add_vertex({1,0,0}), c=open.add_vertex({0,1,0});
    open.add_face(a,b,c);
    rejects([&] { rar::make_sdf_thickness_field(open,cfg); });
    auto degenerate=mesh;
    degenerate.point(rar::Mesh::Vertex_index(1))=degenerate.point(rar::Mesh::Vertex_index(0));
    rejects([&] { rar::validate_sdf_thickness_mesh(degenerate); });
    cfg.sdf_rays=0;
    rejects([&] { rar::validate_sdf_thickness_config(cfg); });
    cfg.sdf_rays=25;
    for (double angle : {0.0,180.0,std::numeric_limits<double>::quiet_NaN()}) {
        cfg.sdf_cone_angle_degrees=angle;
        rejects([&] { rar::validate_sdf_thickness_config(cfg); });
    }
}
} // namespace

int main()
{
    try { formula_tests(); field_tests(); std::cout << "SDF thickness tests passed\n"; }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
