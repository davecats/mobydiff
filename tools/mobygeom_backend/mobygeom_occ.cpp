#include <BRepBndLib.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <Bnd_Box.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Shape.hxx>
#include <TopAbs_State.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double SOLID = 1.0e30;
constexpr double CLASSIFY_TOL = 1.0e-14;
constexpr double INTERSECT_TOL = 1.0e-14;
constexpr int VAR_U = 1;
constexpr int VAR_V = 2;
constexpr int VAR_W = 3;

struct Options {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    double lx = 1.0;
    double ly = 1.0;
    double lz = 1.0;
    double re = 1.0;
};

[[noreturn]] void die(const std::string& msg) {
    throw std::runtime_error(msg);
}

double parse_double(const char* s, const std::string& name) {
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || *end != '\0') die("invalid real value for " + name + ": " + s);
    return v;
}

int parse_int(const char* s, const std::string& name) {
    char* end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0' || v <= 0) die("invalid positive integer for " + name + ": " + s);
    return static_cast<int>(v);
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

bool has_suffix(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    return lower(s.substr(s.size() - suffix.size())) == suffix;
}

TopoDS_Shape read_shape(const std::string& file) {
    if (has_suffix(file, ".step") || has_suffix(file, ".stp")) {
        STEPControl_Reader reader;
        const IFSelect_ReturnStatus stat = reader.ReadFile(file.c_str());
        if (stat != IFSelect_RetDone) die("could not read STEP file: " + file);
        const int n = reader.TransferRoots();
        if (n <= 0) die("STEP file did not transfer any roots: " + file);
        TopoDS_Shape shape = reader.OneShape();
        if (shape.IsNull()) die("STEP file produced a null shape: " + file);
        return shape;
    }

    if (has_suffix(file, ".iges") || has_suffix(file, ".igs")) {
        IGESControl_Reader reader;
        const IFSelect_ReturnStatus stat = reader.ReadFile(file.c_str());
        if (stat != IFSelect_RetDone) die("could not read IGES file: " + file);
        const int n = reader.TransferRoots();
        if (n <= 0) die("IGES file did not transfer any roots: " + file);
        TopoDS_Shape shape = reader.OneShape();
        if (shape.IsNull()) die("IGES file produced a null shape: " + file);
        return shape;
    }

    die("unsupported geometry format for file: " + file);
}

class CadBody {
public:
    explicit CadBody(const TopoDS_Shape& shape_in) : shape(shape_in) {
        classifier.Load(shape);
        intersector.Load(shape, INTERSECT_TOL);

        Bnd_Box box;
        BRepBndLib::Add(shape, box);
        double xmin, ymin, zmin, xmax, ymax, zmax;
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
        const double dx = xmax - xmin;
        const double dy = ymax - ymin;
        const double dz = zmax - zmin;
        ray_length = 10.0 * std::sqrt(dx*dx + dy*dy + dz*dz) + 1.0;
        if (!(ray_length > 0.0)) ray_length = 10.0;
    }

    bool is_inside(const std::array<double,3>& x) {
        const gp_Pnt p(x[0], x[1], x[2]);
        classifier.Perform(p, CLASSIFY_TOL);
        const TopAbs_State state = classifier.State();
        if (state == TopAbs_IN) return true;
        if (state == TopAbs_ON) return false;

        // IGES often imports a closed surface/shell rather than a TopoDS solid.
        // Fall back to odd/even ray casting through the CAD faces.
        return ray_inside(p);
    }

    bool segment_first_intersection_distance(const std::array<double,3>& xa,
                                             const std::array<double,3>& xb,
                                             double& distance) {
        const gp_Pnt p0(xa[0], xa[1], xa[2]);
        const gp_Pnt p1(xb[0], xb[1], xb[2]);
        gp_Vec v(p0, p1);
        const double len = v.Magnitude();
        if (!(len > 0.0)) return false;

        gp_Lin line(p0, gp_Dir(v));
        intersector.Perform(line, 0.0, len);

        std::vector<double> params;
        params.reserve(static_cast<std::size_t>(std::max(0, intersector.NbPnt())));
        for (int i = 1; i <= intersector.NbPnt(); ++i) {
            const double u = intersector.WParameter(i);
            if (u > INTERSECT_TOL && u < len + INTERSECT_TOL) params.push_back(std::min(u, len));
        }
        if (params.empty()) return false;
        std::sort(params.begin(), params.end());
        distance = params.front();
        return true;
    }

private:
    TopoDS_Shape shape;
    BRepClass3d_SolidClassifier classifier;
    IntCurvesFace_ShapeIntersector intersector;
    double ray_length = 10.0;

    bool ray_inside_one_direction(const gp_Pnt& p, const gp_Dir& dir) {
        gp_Lin ray(p, dir);
        intersector.Perform(ray, 0.0, ray_length);

        std::vector<double> params;
        params.reserve(static_cast<std::size_t>(std::max(0, intersector.NbPnt())));
        for (int i = 1; i <= intersector.NbPnt(); ++i) {
            const double u = intersector.WParameter(i);
            if (u > INTERSECT_TOL && u < ray_length) params.push_back(u);
        }
        if (params.empty()) return false;

        std::sort(params.begin(), params.end());
        int unique_count = 0;
        double previous = -1.0e300;
        const double unique_tol = 1.0e-9;
        for (double u : params) {
            if (std::abs(u - previous) > unique_tol) {
                ++unique_count;
                previous = u;
            }
        }
        return (unique_count % 2) == 1;
    }

    bool ray_inside(const gp_Pnt& p) {
        const std::array<gp_Dir, 5> dirs = {
            gp_Dir(0.9363291775690445, 0.2897841486884300, 0.1989788975644836),
            gp_Dir(0.2177215644661700, 0.9513782777304800, 0.2177215644661700),
            gp_Dir(0.1391731009600654, 0.3332056395026950, 0.9324867763917516),
            gp_Dir(-0.7814035730506970, 0.5372149566255140, 0.3177263539173010),
            gp_Dir(0.4375949744936837, -0.8123906677064160, 0.3850575775535940)
        };

        int inside_votes = 0;
        for (const gp_Dir& dir : dirs) {
            if (ray_inside_one_direction(p, dir)) ++inside_votes;
        }
        return inside_votes > static_cast<int>(dirs.size() / 2);
    }
};

bool is_face_staggered(int dir, int var) {
    return dir == var;
}

double coord(int idx, int dir, int var, int n, double length) {
    const double h = length / static_cast<double>(n);
    if (is_face_staggered(dir, var)) return (static_cast<double>(idx) - 1.0) * h;
    return (static_cast<double>(idx) - 0.5) * h;
}

std::array<double,3> point_for_index(int i, int j, int k, int var, const Options& opt) {
    return {coord(i, 1, var, opt.nx, opt.lx),
            coord(j, 2, var, opt.ny, opt.ly),
            coord(k, 3, var, opt.nz, opt.lz)};
}

void add_neighbor_coeff(double& coeff, CadBody& body,
                        const std::array<double,3>& xa,
                        const std::array<double,3>& xb) {
    if (!body.is_inside(xb)) return;

    const double dx = xb[0] - xa[0];
    const double dy = xb[1] - xa[1];
    const double dz = xb[2] - xa[2];
    const double d0 = std::sqrt(dx*dx + dy*dy + dz*dz);

    double d = 0.0;
    if (!body.segment_first_intersection_distance(xa, xb, d)) {
        die("failed to intersect a fluid-solid neighbor segment with the CAD body");
    }
    coeff += ((d0 - d) / d) / (d0 * d0);
}

void write_coeff_raw(const std::string& geometry, const std::string& output, const Options& opt) {
    TopoDS_Shape shape = read_shape(geometry);
    CadBody body(shape);

    const std::uint64_t ni = static_cast<std::uint64_t>(opt.nx + 2);
    const std::uint64_t nj = static_cast<std::uint64_t>(opt.ny + 2);
    const std::uint64_t nk = static_cast<std::uint64_t>(opt.nz + 2);
    const std::uint64_t count = 3ULL * ni * nj * nk;
    const double re_inv = 1.0 / opt.re;
    const double solid_coef = SOLID * re_inv;

    std::ofstream out(output, std::ios::binary);
    if (!out) die("could not open output raw file: " + output);

    const char magic[8] = {'M','O','B','Y','I','B','M','1'};
    const std::int32_t header[4] = {opt.nx, opt.ny, opt.nz, 3};
    out.write(magic, sizeof(magic));
    out.write(reinterpret_cast<const char*>(header), sizeof(header));
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (int var = VAR_U; var <= VAR_W; ++var) {
        for (int k = 0; k <= opt.nz + 1; ++k) {
            for (int j = 0; j <= opt.ny + 1; ++j) {
                for (int i = 0; i <= opt.nx + 1; ++i) {
                    const std::array<double,3> xa = point_for_index(i, j, k, var, opt);
                    double value = 0.0;
                    if (body.is_inside(xa)) {
                        value = solid_coef;
                    } else {
                        double coeff = 0.0;
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i-1, j,   k,   var, opt));
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i+1, j,   k,   var, opt));
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i,   j-1, k,   var, opt));
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i,   j+1, k,   var, opt));
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i,   j,   k-1, var, opt));
                        add_neighbor_coeff(coeff, body, xa, point_for_index(i,   j,   k+1, var, opt));
                        value = coeff * re_inv;
                    }
                    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
                }
            }
        }
    }

    if (!out) die("failed while writing output raw file: " + output);
}

int make_sphere_iges(int argc, char** argv) {
    if (argc != 7) die("usage: mobygeom_occ make-sphere-iges output.igs cx cy cz radius");
    const std::string output = argv[2];
    const double cx = parse_double(argv[3], "cx");
    const double cy = parse_double(argv[4], "cy");
    const double cz = parse_double(argv[5], "cz");
    const double radius = parse_double(argv[6], "radius");
    if (!(radius > 0.0)) die("sphere radius must be positive");

    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gp_Pnt(cx, cy, cz), radius).Shape();
    IGESControl_Writer writer("MM", 0);
    writer.AddShape(sphere);
    if (!writer.Write(output.c_str())) die("failed to write IGES sphere: " + output);
    return 0;
}


int make_sphere_step(int argc, char** argv) {
    if (argc != 7) die("usage: mobygeom_occ make-sphere-step output.step cx cy cz radius");
    const std::string output = argv[2];
    const double cx = parse_double(argv[3], "cx");
    const double cy = parse_double(argv[4], "cy");
    const double cz = parse_double(argv[5], "cz");
    const double radius = parse_double(argv[6], "radius");
    if (!(radius > 0.0)) die("sphere radius must be positive");

    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gp_Pnt(cx, cy, cz), radius).Shape();
    STEPControl_Writer writer;
    IFSelect_ReturnStatus stat = writer.Transfer(sphere, STEPControl_AsIs);
    if (stat != IFSelect_RetDone) die("failed to transfer STEP sphere: " + output);
    stat = writer.Write(output.c_str());
    if (stat != IFSelect_RetDone) die("failed to write STEP sphere: " + output);
    return 0;
}

int coeff(int argc, char** argv) {
    if (argc != 11) {
        die("usage: mobygeom_occ coeff geometry.step|iges output.raw nx ny nz lx ly lz re");
    }
    const std::string geometry = argv[2];
    const std::string output = argv[3];
    Options opt;
    opt.nx = parse_int(argv[4], "nx");
    opt.ny = parse_int(argv[5], "ny");
    opt.nz = parse_int(argv[6], "nz");
    opt.lx = parse_double(argv[7], "lx");
    opt.ly = parse_double(argv[8], "ly");
    opt.lz = parse_double(argv[9], "lz");
    opt.re = parse_double(argv[10], "re");
    if (!(opt.lx > 0.0 && opt.ly > 0.0 && opt.lz > 0.0)) die("domain lengths must be positive");
    if (!(opt.re > 0.0)) die("Reynolds number must be positive");

    write_coeff_raw(geometry, output, opt);
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) die("usage: mobygeom_occ <make-sphere-iges|coeff> ...");
        const std::string cmd = argv[1];
        if (cmd == "make-sphere-iges") return make_sphere_iges(argc, argv);
        if (cmd == "make-sphere-step") return make_sphere_step(argc, argv);
        if (cmd == "coeff") return coeff(argc, argv);
        die("unknown command: " + cmd);
    } catch (const Standard_Failure& e) {
        std::cerr << "OCCT error: " << e.GetMessageString() << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
