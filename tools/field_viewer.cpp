#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct FieldMesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<std::size_t, 3>> faces;
    std::vector<double> curvature;
    std::vector<double> target_length;
};

FieldMesh read_field_ply(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open: " + path);
    }

    std::string line;
    if (!std::getline(in, line) || line != "ply") {
        throw std::runtime_error("Not a PLY file: " + path);
    }

    bool ascii = false;
    bool in_vertex = false;
    std::size_t vertex_count = 0;
    std::size_t face_count = 0;
    std::vector<std::string> vertex_properties;

    while (std::getline(in, line)) {
        if (line == "end_header") {
            break;
        }

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "format") {
            std::string format_name;
            iss >> format_name;
            ascii = (format_name == "ascii");
        } else if (token == "element") {
            std::string element_name;
            std::size_t count = 0;
            iss >> element_name >> count;
            in_vertex = (element_name == "vertex");

            if (element_name == "vertex") {
                vertex_count = count;
            } else if (element_name == "face") {
                face_count = count;
            }
        } else if (token == "property" && in_vertex) {
            std::string type_or_list;
            iss >> type_or_list;
            if (type_or_list == "list") {
                throw std::runtime_error(
                    "Vertex list properties are not supported.");
            }

            std::string property_name;
            iss >> property_name;
            vertex_properties.push_back(property_name);
        }
    }

    if (!ascii) {
        throw std::runtime_error(
            "Field viewer currently supports ASCII PLY only.");
    }

    auto property_index = [&](const std::string& name) -> std::size_t {
        const auto it = std::find(
            vertex_properties.begin(),
            vertex_properties.end(),
            name);
        if (it == vertex_properties.end()) {
            throw std::runtime_error(
                "Missing vertex property '" + name + "'.");
        }
        return static_cast<std::size_t>(
            std::distance(vertex_properties.begin(), it));
    };

    const std::size_t ix = property_index("x");
    const std::size_t iy = property_index("y");
    const std::size_t iz = property_index("z");
    const std::size_t ic = property_index("curvature");
    const std::size_t itarget = property_index("target_length");

    FieldMesh data;
    data.vertices.reserve(vertex_count);
    data.curvature.reserve(vertex_count);
    data.target_length.reserve(vertex_count);

    for (std::size_t i = 0; i < vertex_count; ++i) {
        if (!std::getline(in, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading vertices.");
        }

        std::istringstream iss(line);
        std::vector<double> values;
        double value = 0.0;
        while (iss >> value) {
            values.push_back(value);
        }

        if (values.size() < vertex_properties.size()) {
            throw std::runtime_error(
                "Vertex row has fewer values than the PLY header declares.");
        }

        data.vertices.push_back({
            values[ix],
            values[iy],
            values[iz]
        });
        data.curvature.push_back(values[ic]);
        data.target_length.push_back(values[itarget]);
    }

    data.faces.reserve(face_count);

    for (std::size_t i = 0; i < face_count; ++i) {
        if (!std::getline(in, line)) {
            throw std::runtime_error(
                "Unexpected end of file while reading faces.");
        }

        std::istringstream iss(line);
        std::size_t n = 0;
        iss >> n;

        if (n < 3) {
            continue;
        }

        std::vector<std::size_t> indices(n);
        for (std::size_t j = 0; j < n; ++j) {
            if (!(iss >> indices[j])) {
                throw std::runtime_error(
                    "Malformed face record in PLY.");
            }
        }

        // The exported field mesh is triangular, but fan triangulation keeps
        // the viewer tolerant of polygonal ASCII PLY inputs.
        for (std::size_t j = 1; j + 1 < n; ++j) {
            data.faces.push_back({
                indices[0],
                indices[j],
                indices[j + 1]
            });
        }
    }

    if (data.vertices.empty() || data.faces.empty()) {
        throw std::runtime_error(
            "PLY contains no usable surface mesh.");
    }

    return data;
}

std::vector<double> refinement_demand(
    const std::vector<double>& target_length)
{
    const auto [min_it, max_it] = std::minmax_element(
        target_length.begin(),
        target_length.end());

    const double lo = *min_it;
    const double hi = *max_it;
    const double span = hi - lo;

    std::vector<double> result(target_length.size(), 0.5);
    if (!(span > 0.0)) {
        return result;
    }

    for (std::size_t i = 0; i < target_length.size(); ++i) {
        // 1 = smallest target length = strongest refinement demand.
        result[i] = (hi - target_length[i]) / span;
    }
    return result;
}

void print_usage(const char* exe) {
    std::cerr
        << "Usage:\n  "
        << exe
        << " <field.ply>\n\n"
        << "The input should be the field.ply generated by rar_cgal.\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        const std::string input_path = argv[1];
        const FieldMesh data = read_field_ply(input_path);
        const std::vector<double> refinement =
            refinement_demand(data.target_length);

        polyscope::options::programName = "RAR Field Viewer";
        polyscope::options::verbosity = 0;
        polyscope::init();

        auto* mesh = polyscope::registerSurfaceMesh(
            "field mesh",
            data.vertices,
            data.faces);

        // Flat shading + visible edges makes remeshing artifacts easier to see.
        mesh->setSmoothShade(false);
        mesh->setEdgeWidth(1.0);

        auto* target_q = mesh->addVertexScalarQuantity(
            "target_length",
            data.target_length);
        target_q->setColorMap("turbo");
        target_q->setOnscreenColorbarEnabled(true);

        auto* curvature_q = mesh->addVertexScalarQuantity(
            "curvature",
            data.curvature);
        curvature_q->setColorMap("turbo");
        curvature_q->setOnscreenColorbarEnabled(true);

        auto* refinement_q = mesh->addVertexScalarQuantity(
            "refinement_demand",
            refinement);
        refinement_q->setColorMap("turbo");
        refinement_q->setMapRange({0.0, 1.0});
        refinement_q->setOnscreenColorbarEnabled(true);

        // Start in the most intuitive mode for sizing-field inspection:
        // red/high means stronger refinement demand.
        refinement_q->setEnabled(true);

        std::cout
            << "Loaded: " << input_path << "\n"
            << "Vertices: " << data.vertices.size() << "\n"
            << "Faces: " << data.faces.size() << "\n"
            << "\n"
            << "Viewer tips:\n"
            << "  - Switch target_length / curvature / refinement_demand "
               "in the mesh quantity panel.\n"
            << "  - Drag or type the scalar map min/max below the colorbar.\n"
            << "  - Change colormap from the quantity options.\n"
            << "  - Change edge width, smooth shading, and back-face policy "
               "from mesh options.\n"
            << "  - Click mesh elements to inspect values.\n";

        polyscope::show();
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
