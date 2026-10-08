#include "rar/Types.h"
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/helpers.h>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    rar::Mesh mesh;
    namespace PMP = CGAL::Polygon_mesh_processing;
    if (!PMP::IO::read_polygon_mesh(argv[1],mesh) || !CGAL::is_triangle_mesh(mesh)) return 3;
    const bool closed = CGAL::is_closed(mesh);
    const bool self_intersects = PMP::does_self_intersect(mesh);
    const bool bounds_volume = closed && !self_intersects && PMP::does_bound_a_volume(mesh);
    std::cout << "{\"closed\":" << (closed ? "true":"false")
              << ",\"self_intersects\":" << (self_intersects ? "true":"false")
              << ",\"bounds_volume\":" << (bounds_volume ? "true":"false") << "}\n";
}
