#include "VoronoiDensityShells.hpp"
#include "../BoundingBox.hpp"
#include "../Point.hpp"
#include "../KDTreeIndirect.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Random.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Named_function_parameters.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <functional>
#include <optional>
#include <CGAL/number_utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Voronoi {
namespace {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;

Vec3d vector(const Kernel::Point_3 &p) { return {CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())}; }
Kernel::Point_3 point(const Vec3d &p) { return {p.x(), p.y(), p.z()}; }

Eigen::Matrix3d fcc_basis(double angle)
{
    // FCC primitive vectors. The six-tetrahedron cell decomposition contains
    // two regular tetrahedra and four pieces of the neighbouring octahedra.
    Eigen::Matrix3d primitive, orientation;
    primitive << 1., -1., -1., 1., 0., 0., 0., -1., 1.;
    // Put a close-packed (111) plane parallel to the bed, then apply XY yaw.
    orientation.row(0) = Vec3d(1., -1., 0.).normalized();
    orientation.row(1) = Vec3d(1., 1., -2.).normalized();
    orientation.row(2) = Vec3d(1., 1., 1.).normalized();
    return Eigen::AngleAxisd(angle, Vec3d::UnitZ()).toRotationMatrix() * orientation * primitive / std::sqrt(2.);
}

Vec3d grid_position(const Eigen::Matrix3d &basis, const Vec3d &origin, double step, const Vec3d &index)
{
    return basis * (origin + step * index);
}

// Grid edges share mesh vertices, including diagonals, so adjoining cells
// stay connected. The coordinates come from the FCC primitive basis.
Mesh extract(const std::vector<float> &field, const Vec3d &origin, double step,
             const std::array<size_t, 3> &size, double level,
             const Eigen::Matrix3d &basis,
             const std::function<void(double)> &progress)
{
    Mesh mesh;
    std::map<std::pair<size_t, size_t>, Mesh::Vertex_index> crossings;
    constexpr int tetrahedra[6][4] = {{0,1,3,7}, {0,3,2,7}, {0,2,6,7}, {0,6,4,7}, {0,4,5,7}, {0,5,1,7}};
    const size_t plane = size[0] * size[1];
    for (size_t z = 0; z + 1 < size[2]; ++z) {
        for (size_t y = 0; y + 1 < size[1]; ++y)
            for (size_t x = 0; x + 1 < size[0]; ++x) {
                std::array<size_t, 8> ids;
                std::array<Vec3d, 8> positions;
                std::array<double, 8> values;
                for (size_t i = 0; i < 8; ++i) {
                    const size_t ix = x + (i & 1), iy = y + ((i >> 1) & 1), iz = z + ((i >> 2) & 1);
                    ids[i] = iz * plane + iy * size[0] + ix;
                    positions[i] = grid_position(basis, origin, step, Vec3d(double(ix), double(iy), double(iz)));
                    // Avoid crossings landing exactly at a shared grid vertex.
                    values[i] = double(field[ids[i]]) - level;
                    if (std::abs(values[i]) < 1e-7) values[i] = 1e-7;
                }
                if (std::all_of(values.begin(), values.end(), [](double v) { return v > 0.; }) ||
                    std::all_of(values.begin(), values.end(), [](double v) { return v < 0.; })) continue;
                const auto vertex = [&](int a, int b) {
                    const auto key = std::minmax(ids[a], ids[b]);
                    auto [it, inserted] = crossings.emplace(key, Mesh::null_vertex());
                    if (inserted)
                        it->second = mesh.add_vertex(point(positions[a] + (positions[b] - positions[a]) *
                                                           (values[a] / (values[a] - values[b]))));
                    return it->second;
                };
                for (const auto &tet : tetrahedra) {
                    std::vector<int> inside, outside;
                    for (int i : tet)
                        (values[i] < 0. ? inside : outside).push_back(i);
                    if (inside.empty() || outside.empty()) continue;
                    Vec3d direction = Vec3d::Zero();
                    for (int i : outside) direction += positions[i] / double(outside.size());
                    for (int i : inside) direction -= positions[i] / double(inside.size());
                    const auto triangle = [&](Mesh::Vertex_index a, Mesh::Vertex_index b, Mesh::Vertex_index c) {
                        const Vec3d pa = vector(mesh.point(a)), pb = vector(mesh.point(b)), pc = vector(mesh.point(c));
                        if ((pb - pa).cross(pc - pa).dot(direction) < 0.) std::swap(b, c);
                        if (mesh.add_face(a, b, c) == Mesh::null_face())
                            throw std::runtime_error("Non-manifold Voronoi density surface");
                    };
                    if (inside.size() == 2) {
                        const auto a = vertex(inside[0], outside[0]), b = vertex(inside[0], outside[1]);
                        const auto c = vertex(inside[1], outside[0]), d = vertex(inside[1], outside[1]);
                        triangle(a, b, d);
                        triangle(a, d, c);
                    } else {
                        if (inside.size() > outside.size()) std::swap(inside, outside);
                        triangle(vertex(inside[0], outside[0]), vertex(inside[0], outside[1]), vertex(inside[0], outside[2]));
                    }
                }
            }
        if (progress) progress(double(z + 1) / (size[2] - 1));
        // Intersections below the current slab cannot be needed again.
        for (auto it = crossings.begin(); it != crossings.end(); )
            if (it->first.second < (z + 1) * plane) it = crossings.erase(it); else ++it;
    }
    return mesh;
}

struct ShellFields {
    std::vector<float> raw, phase;
};

// Solve |gradient(phase)| = local density ratio. Equal phase intervals then
// have the same local spacing as the surface remesher. FCC contains four
// interleaved cubic grids; these offsets follow their orthogonal axes.
std::vector<float> spacing_field(const std::vector<float> &raw, const std::vector<float> &ratio,
                                 double step, const std::array<size_t, 3> &size,
                                 std::vector<float> phase,
                                 const std::function<void(double)> &progress)
{
    const float infinity = std::numeric_limits<float>::infinity();
    if (std::all_of(ratio.begin(), ratio.end(), [&](float r) { return r == ratio.front(); })) {
        for (size_t i = 0; i < raw.size(); ++i) phase[i] = raw[i] * ratio[i];
        progress(1.);
        return phase;
    }
    constexpr int directions[3][3] = {{0, 1, 1}, {2, 1, 1}, {0, -1, 1}};
    const double h = std::sqrt(2.) * step;
    const auto neighbour = [&](size_t i, int axis, int sign) -> size_t {
        const std::array<ptrdiff_t, 3> p {ptrdiff_t(i % size[0]), ptrdiff_t((i / size[0]) % size[1]),
                                        ptrdiff_t(i / (size[0] * size[1]))};
        size_t result = 0, stride = 1;
        for (int dimension = 0; dimension < 3; ++dimension) {
            const ptrdiff_t coordinate = p[dimension] + sign * directions[axis][dimension];
            if (coordinate < 0 || coordinate >= ptrdiff_t(size[dimension])) return size_t(-1);
            result += size_t(coordinate) * stride;
            stride *= size[dimension];
        }
        return result;
    };
    using Entry = std::pair<float, size_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pending;
    std::vector<bool> accepted(raw.size(), false);
    size_t completed = 0;
    for (size_t i = 0; i < raw.size(); ++i) {
        if (std::abs(raw[i]) <= h) pending.emplace(phase[i], i);
        if (i % 4096 == 0) progress(0.);
    }
    const auto update = [&](size_t i) {
        if (accepted[i] || std::abs(raw[i]) <= h) return;
        std::array<std::pair<double, double>, 3> values;
        for (int axis = 0; axis < 3; ++axis) {
            size_t nearest = size_t(-1);
            int direction = 0;
            for (int sign : {-1, 1}) {
                const size_t j = neighbour(i, axis, sign);
                if (j != size_t(-1) && accepted[j] && (raw[j] == 0. || (raw[i] > 0.) == (raw[j] > 0.)) &&
                    (nearest == size_t(-1) || phase[j] < phase[nearest])) {
                    nearest = j;
                    direction = sign;
                }
            }
            values[axis] = {nearest == size_t(-1) ? infinity : double(phase[nearest]), 1.};
            if (nearest == size_t(-1)) continue;
            const size_t second = neighbour(nearest, axis, direction);
            // Second-order upwinding reduces spacing bias through a decay
            // gradient. Fall back to first order beside a boundary or ridge.
            if (second != size_t(-1) && accepted[second] && (raw[second] == 0. || (raw[i] > 0.) == (raw[second] > 0.)) &&
                phase[second] <= phase[nearest])
                values[axis] = {(4. * phase[nearest] - phase[second]) / 3., 2.25};
        }
        std::sort(values.begin(), values.end());
        if (!std::isfinite(values[0].first)) return;
        const double cost = h * ratio[i];
        double weight = 0., sum = 0., squares = 0., candidate = infinity;
        for (const auto &[value, w] : values) {
            if (!std::isfinite(value) || candidate <= value) break;
            weight += w;
            sum += w * value;
            squares += w * value * value;
            candidate = (sum + std::sqrt(std::max(0., sum * sum - weight * (squares - cost * cost)))) / weight;
        }
        if (float(candidate) < phase[i]) {
            phase[i] = float(candidate);
            pending.emplace(phase[i], i);
        }
    };
    while (!pending.empty()) {
        const auto [value, i] = pending.top();
        pending.pop();
        if (accepted[i] || value != phase[i]) continue;
        accepted[i] = true;
        for (int axis = 0; axis < 3; ++axis)
            for (int sign : {-1, 1}) {
                const size_t j = neighbour(i, axis, sign);
                if (j != size_t(-1)) update(j);
            }
        if (++completed % 4096 == 0) progress(double(completed) / raw.size());
    }
    if (std::any_of(phase.begin(), phase.end(), [](float v) { return !std::isfinite(v); }))
        throw std::runtime_error("Incomplete Voronoi shell spacing field");
    for (size_t i = 0; i < raw.size(); ++i) if (raw[i] < 0.) phase[i] = -phase[i];
    progress(1.);
    return phase;
}

// Alternate extraction grids sample the continuous distance field directly;
// interpolating the shell grid again would blur away small inner surfaces.
ShellFields sample_grid(const WallDistanceSamples &sample, const Vec3d &origin, double step,
                        const std::array<size_t, 3> &size, const Eigen::Matrix3d &basis,
                        double decay,
                        const std::function<void(double)> &progress)
{
    const size_t plane = size[0] * size[1];
    const size_t count = plane * size[2];
    const double h = std::sqrt(2.) * step;
    ShellFields field {std::vector<float>(count), std::vector<float>(count, std::numeric_limits<float>::infinity())};
    std::vector<float> ratio(count);
    for (size_t z = 0; z < size[2]; ++z) {
        std::vector<Vec3d> points;
        points.reserve(plane);
        for (size_t y = 0; y < size[1]; ++y)
            for (size_t x = 0; x < size[0]; ++x)
                points.push_back(grid_position(basis, origin, step, Vec3d(double(x), double(y), double(z))));
        const auto values = sample(points);
        if (values.size() != points.size()) throw std::invalid_argument("Invalid Voronoi distance samples");
        for (size_t i = 0; i < values.size(); ++i) {
            const double depth = -values[i].distance;
            const size_t index = z * plane + i;
            field.raw[index] = float(depth);
            const auto &value = values[i];
            ratio[index] = float(depth > 0. ?
                WallDistancePointCloud::density_ratio(value.distance, value.smoothed, decay) :
                WallDistancePointCloud::density_ratio(0., value.smoothed - value.distance, decay));
            if (std::abs(field.raw[index]) <= h) {
                double integral = std::abs(depth) * ratio[index];
                if (depth > 0.) {
                    // Integrate the near-boundary band using the local corner
                    // correction; it varies slowly over one grid step.
                    integral = 0.;
                    for (int s = 0; s < 8; ++s) {
                        const double d = depth * (s + 0.5) / 8.;
                        integral += depth / 8. * WallDistancePointCloud::density_ratio(
                            -d, -d + value.smoothed - value.distance, decay);
                    }
                }
                field.phase[index] = float(integral);
            }
        }
        progress(0.5 * double(z + 1) / size[2]);
    }
    field.phase = spacing_field(field.raw, ratio, step, size, std::move(field.phase),
                                [&](double f) { progress(0.5 + 0.5 * f); });
    return field;
}

// A narrow component can disappear between contour levels. Fill only gaps
// larger than local spacing, using interior sites so exterior guards cannot
// mask a missing neck. The extraction grid sets the coverage resolution.
void cover_gaps(std::vector<Vec3d> &sites, const WallDistanceSamples &sample,
                const std::vector<float> &depth, const Vec3d &origin, double step,
                const std::array<size_t, 3> &size, double spacing, double decay,
                const Eigen::Matrix3d &basis,
                const std::function<void(double)> &progress)
{
    std::vector<size_t> interior;
    constexpr size_t batch = 4096;
    for (size_t first = 0; first < sites.size(); first += batch) {
        const size_t last = std::min(first + batch, sites.size());
        const auto fields = sample(std::vector<Vec3d>(sites.begin() + first, sites.begin() + last));
        for (size_t i = 0; i < fields.size(); ++i)
            if (fields[i].distance <= 1e-6) interior.push_back(first + i);
        if (progress) progress(0.);
    }
    const auto coordinate = [&](size_t i, size_t axis) { return sites[i][axis]; };
    KDTreeIndirect<3, double, decltype(coordinate)> tree(coordinate);
    const bool populated = !interior.empty();
    if (populated) tree.build(interior);
    const auto position = [&](size_t i) {
        return grid_position(basis, origin, step, Vec3d(double(i % size[0]), double((i / size[0]) % size[1]),
                                                        double(i / (size[0] * size[1]))));
    };
    struct Gap { size_t index; double radius; };
    std::vector<Gap> gaps;
    for (size_t first = 0; first < depth.size(); first += batch) {
        const size_t last = std::min(first + batch, depth.size());
        std::vector<Vec3d> points;
        std::vector<size_t> indices;
        for (size_t i = first; i < last; ++i)
            if (depth[i] > 0.) { points.push_back(position(i)); indices.push_back(i); }
        const auto fields = sample(points);
        for (size_t i = 0; i < points.size(); ++i) {
            const double radius = spacing / WallDistancePointCloud::density_ratio(fields[i].distance, fields[i].smoothed, decay);
            const size_t nearest = populated ? find_closest_point(tree, points[i]) : size_t(-1);
            if (nearest == size_t(-1) || (sites[nearest] - points[i]).squaredNorm() > radius * radius)
                gaps.push_back({indices[i], radius});
        }
        if (progress) progress(0.5 * double(last) / depth.size());
    }
    // Prefer the deepest part of each uncovered pocket, then suppress nearby
    // candidates rather than introducing a second dense grid of control points.
    std::sort(gaps.begin(), gaps.end(), [&](const Gap &a, const Gap &b) {
        if (depth[a.index] != depth[b.index]) return depth[a.index] > depth[b.index];
        return a.index < b.index;
    });
    const auto gap_coordinate = [&](size_t i, size_t axis) { return position(gaps[i].index)[axis]; };
    KDTreeIndirect<3, double, decltype(gap_coordinate)> gap_tree(gap_coordinate, gaps.size());
    std::vector<bool> covered(gaps.size(), false);
    CGAL::Random random(0x5348454c);
    for (size_t i = 0; i < gaps.size(); ++i) {
        if (!covered[i]) {
            Vec3d p = position(gaps[i].index);
            const double jitter = 0.15 * std::min(step, gaps[i].radius);
            Vec3d shifted = p;
            for (int axis = 0; axis < 3; ++axis) shifted[axis] += random.get_double(-jitter, jitter);
            if (sample({shifted}).at(0).distance < 0.) p = shifted;
            sites.push_back(p);
            for (size_t j : find_nearby_points(gap_tree, p, gaps[i].radius))
                if ((position(gaps[j].index) - p).squaredNorm() <= gaps[j].radius * gaps[j].radius)
                    covered[j] = true;
        }
        if (progress && (i % batch == 0 || i + 1 == gaps.size())) progress(0.5 + 0.5 * double(i + 1) / gaps.size());
    }
    if (progress) progress(1.);
}

// Adaptive sizing preserves corner boosts even when decay is disabled.
class SurfaceSpacing {
public:
    using FT = double;
    using Point_3 = Kernel::Point_3;
    SurfaceSpacing(const WallDistanceSamples &sample, double spacing, double decay, bool exterior, std::function<void()> cancel)
        : m_sample(sample), m_spacing(spacing), m_decay(decay), m_exterior(exterior), m_cancel(std::move(cancel)) {}
    FT at(Mesh::Vertex_index v, const Mesh &mesh) const { return target(mesh.point(v)); }
    std::optional<FT> is_too_long(Mesh::Vertex_index a, Mesh::Vertex_index b, const Mesh &mesh) const
    {
        const double length = (vector(mesh.point(a)) - vector(mesh.point(b))).norm();
        const double limit = 4. / 3. * target(CGAL::midpoint(mesh.point(a), mesh.point(b)));
        return length > limit ? std::optional<FT>((length / limit) * (length / limit)) : std::nullopt;
    }
    std::optional<FT> is_too_short(Mesh::Halfedge_index h, const Mesh &mesh) const
    {
        const auto a = mesh.source(h), b = mesh.target(h);
        const double length = (vector(mesh.point(a)) - vector(mesh.point(b))).norm();
        const double limit = 4. / 5. * target(CGAL::midpoint(mesh.point(a), mesh.point(b)));
        return length < limit ? std::optional<FT>((length / limit) * (length / limit)) : std::nullopt;
    }
    Point_3 split_placement(Mesh::Halfedge_index h, const Mesh &mesh) const
        { return CGAL::midpoint(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))); }
    void register_split_vertex(Mesh::Vertex_index, const Mesh &) const {}
private:
    double target(const Point_3 &p) const
    {
        if (++m_queries % 2048 == 0) m_cancel();
        auto field = m_sample({vector(p)}).at(0);
        if (m_exterior) {
            field.smoothed -= field.distance;
            field.distance = 0.;
        }
        return m_spacing / WallDistancePointCloud::density_ratio(field.distance, field.smoothed, m_decay);
    }
    const WallDistanceSamples &m_sample;
    double m_spacing, m_decay;
    bool m_exterior;
    std::function<void()> m_cancel;
    mutable size_t m_queries = 0;
};
} // namespace

DensityShellPointCloud::DensityShellPointCloud(WallDistanceSamples sample, const BoundingBoxf3 &bounds,
                                             double site_spacing, double decay, const std::function<void(double)> &progress,
                                             double angle)
    : WallDistancePointCloud(sample, decay)
{
    if (!std::isfinite(site_spacing) || site_spacing <= 0.) throw std::invalid_argument("Invalid Voronoi shell spacing");
    if (!std::isfinite(angle)) throw std::invalid_argument("Invalid Voronoi shell angle");
    if (!bounds.defined) return;
    // A triangular grid with spacing a and inter-shell separation a has site
    // intensity 2/(sqrt(3)*a^3). Match the random provider's nominal intensity.
    const double a = std::cbrt(2. / std::sqrt(3.)) * site_spacing;
    const Eigen::Matrix3d basis = fcc_basis(angle);
    const BoundingBoxf3 domain = BoundingBoxf3(bounds.min - Vec3d::Constant(2. * a), bounds.max + Vec3d::Constant(2. * a))
                                    .transformed(Transform3d(basis.inverse()));
    double step = std::clamp(a / 2., 1., 2.);
    std::array<size_t, 3> size;
    size_t count;
    for (;;) {
        count = 1;
        bool fits = true;
        for (int axis = 0; axis < 3; ++axis) {
            size[axis] = size_t(std::ceil((domain.max[axis] - domain.min[axis]) / step)) + 1;
            if (size[axis] > 2 * 1024 * 1024 / count) { fits = false; break; }
            count *= size[axis];
        }
        if (fits) break;
        step *= 2.;
    }
    const ShellFields fields = sample_grid(sample, domain.min, step, size, basis, decay,
        [&](double fraction) { if (progress) progress(0.1 * fraction); });
    const double maximum = *std::max_element(fields.phase.begin(), fields.phase.end());
    std::vector<double> levels {-a, 0.};
    for (double level = a; level < maximum; level += a) levels.push_back(level);
    ShellFields shifted;
    for (size_t shell = 0; shell < levels.size(); ++shell) {
        const auto report = [&](double fraction) {
            if (progress) progress(0.1 + 0.8 * (shell + fraction) / levels.size());
        };
        // The exterior shell supplies neighbours across the boundary.
        const Vec3d origin = domain.min + Vec3d::Constant(shell % 2 ? 0.5 * step : 0.);
        if (shell % 2 && shifted.raw.empty())
            shifted = sample_grid(sample, origin, step, size, basis, decay, [&](double f) { report(0.15 * f); });
        const auto &field = shell % 2 ? shifted : fields;
        Mesh mesh = extract(levels[shell] == 0. ? field.raw : field.phase, origin, step, size, levels[shell], basis,
                            [&](double f) { report(0.15 + 0.35 * f); });
        if (mesh.number_of_faces() == 0) { report(1.); continue; }
        for (int pass = 0; pass < 3; ++pass) {
            const SurfaceSpacing sizing(sample, a, decay, shell == 0, [&] { report(0.5 + pass / 6.); });
            CGAL::Polygon_mesh_processing::isotropic_remeshing(mesh.faces(), sizing, mesh,
                CGAL::parameters::number_of_iterations(1));
            report(0.5 + (pass + 1) / 6.);
        }
        for (auto vertex : mesh.vertices()) m_sites.push_back(vector(mesh.point(vertex)));
    }
    shifted = {};
    cover_gaps(m_sites, sample, fields.raw, domain.min, step, size, a, decay, basis,
               [&](double fraction) { if (progress) progress(0.9 + 0.1 * fraction); });
    for (const auto &p : m_sites) m_bounds.merge(p);
    std::sort(m_sites.begin(), m_sites.end(), [](const Vec3d &a, const Vec3d &b) {
        if (a.z() != b.z()) return a.z() < b.z();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.x() < b.x();
    });
    if (progress) progress(1.);
}

std::vector<Vec3d> DensityShellPointCloud::points_in(const BoundingBoxf3 &box, double) const
{
    std::vector<Vec3d> points;
    if (!box.defined) return points;
    const auto first = std::lower_bound(m_sites.begin(), m_sites.end(), box.min.z(), [](const Vec3d &p, double z) { return p.z() < z; });
    const auto last = std::lower_bound(first, m_sites.end(), box.max.z(), [](const Vec3d &p, double z) { return p.z() < z; });
    for (auto it = first; it != last; ++it)
        if ((it->array() >= box.min.array()).all() && (it->array() < box.max.array()).all()) points.push_back(*it);
    return points;
}
} // namespace Slic3r::Voronoi
