#include "VoronoiVoid.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "VoronoiPointCloud.hpp"
#include "VoronoiRelaxation.hpp"
#include "../libslic3r.h"
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Voronoi {
namespace {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Delaunay_triangulation_3<Kernel>;
Vec3d point(const Kernel::Point_3 &p) { return {p.x(), p.y(), p.z()}; }
Kernel::Point_3 point(const Vec3d &p) { return {p.x(), p.y(), p.z()}; }

// Clip a bounded dual face, preserving its cyclic vertex order.
std::vector<Vec3d> face_polygon(const Mesh &mesh, const Mesh::Edge &edge, const BoundingBoxf3 &bounds)
{
    std::vector<Vec3d> polygon;
    auto first = mesh.incident_cells(edge), cell = first;
    do {
        if (mesh.is_infinite(cell)) return {};
        const Vec3d p = point(mesh.dual(cell));
        if (!p.allFinite()) return {};
        polygon.push_back(p);
    } while (++cell != first);
    for (int axis = 0; axis < 3; ++axis)
        for (int side = 0; side < 2 && !polygon.empty(); ++side) {
            std::vector<Vec3d> clipped;
            const double limit = side ? bounds.max[axis] : bounds.min[axis];
            auto distance = [&](const Vec3d &p) { return (side ? -1. : 1.) * (p[axis] - limit); };
            Vec3d previous = polygon.back();
            double d0 = distance(previous);
            for (const auto &p : polygon) {
                const double d1 = distance(p);
                if ((d0 < 0.) != (d1 < 0.)) clipped.push_back(previous + (p - previous) * (d0 / (d0 - d1)));
                if (d1 >= 0.) clipped.push_back(p);
                previous = p;
                d0 = d1;
            }
            polygon = std::move(clipped);
        }
    polygon.erase(std::unique(polygon.begin(), polygon.end(), [](const Vec3d &a, const Vec3d &b) {
        return (a - b).norm() <= EPSILON;
    }), polygon.end());
    if (polygon.size() > 1 && (polygon.front() - polygon.back()).norm() <= EPSILON) polygon.pop_back();
    return polygon;
}

struct BadFace {
    Vec3d center, direction;
    double span, cost;
};

std::vector<BadFace> bad_faces(const Mesh &mesh, const BoundingBoxf3 &bounds,
                               const WallDistanceSamples &samples, const VoidRefinement &settings,
                               const std::function<void()> &cancel)
{
    std::vector<BadFace> result;
    if (mesh.dimension() != 3) return result;
    const double cosine = std::cos(settings.minimum_angle * std::acos(-1.) / 180.);
    size_t visited = 0;
    for (auto edge = mesh.finite_edges_begin(); edge != mesh.finite_edges_end(); ++edge) {
        if (++visited % 1024 == 0) cancel();
        const Vec3d a = point(edge->first->vertex(edge->second)->point());
        const Vec3d b = point(edge->first->vertex(edge->third)->point());
        const Vec3d delta = b - a;
        if (std::abs(delta.z()) <= cosine * delta.norm()) continue;
        const auto polygon = face_polygon(mesh, *edge, bounds);
        if (polygon.size() < 3) continue;
        Vec3d center = Vec3d::Zero();
        for (const auto &p : polygon) center += p;
        center /= polygon.size();
        Vec3d direction = Vec3d::Zero();
        double span = 0.;
        for (const auto &p : polygon)
            for (const auto &q : polygon) {
                const double length = (p - q).head<2>().norm();
                if (length > span) { span = length; direction = (p - q).normalized(); }
            }
        if (span <= settings.maximum_span) continue;
        // Integrate only triangles whose centre lies in the sliced model.
        std::vector<Vec3d> probes;
        for (size_t i = 0; i < polygon.size(); ++i)
            probes.push_back((center + polygon[i] + polygon[(i + 1) % polygon.size()]) / 3.);
        const auto values = samples(probes);
        if (values.size() != probes.size()) throw std::invalid_argument("Invalid Voronoi refinement samples");
        double area = 0.;
        for (size_t i = 0; i < polygon.size(); ++i)
            if (values[i].distance < 0.)
                area += 0.5 * (polygon[i] - center).cross(polygon[(i + 1) % polygon.size()] - center).norm();
        const double cost = area * (span / settings.maximum_span - 1.);
        if (cost > EPSILON) result.push_back({center, direction, span, cost});
    }
    // Equal-cost ties use geometry rather than CGAL traversal order.
    std::sort(result.begin(), result.end(), [](const BadFace &a, const BadFace &b) {
        if (a.cost != b.cost) return a.cost > b.cost;
        for (int axis = 0; axis < 3; ++axis)
            if (a.center[axis] != b.center[axis]) return a.center[axis] < b.center[axis];
        return false;
    });
    return result;
}

double cost(const std::vector<BadFace> &faces)
{
    double total = 0.;
    for (const auto &face : faces) total += face.cost;
    return total;
}
} // namespace

VoidPointCloud::VoidPointCloud(std::shared_ptr<const PointCloudProvider> source, const BoundingBoxf3 &bounds,
                               double site_spacing, WallDistanceSamples samples, double threshold,
                               const VoidRefinement &settings, const std::function<void(double)> &progress)
    : m_source(std::move(source))
{
    if (!m_source || !samples || !std::isfinite(threshold) || threshold <= 0. ||
        !std::isfinite(settings.minimum_angle) || settings.minimum_angle <= 0. || settings.minimum_angle >= 90. ||
        !std::isfinite(settings.maximum_span) || settings.maximum_span <= 0.)
        throw std::invalid_argument("Invalid Voronoi void settings");
    const auto report = [&](double fraction) { if (progress) progress(fraction); };
    report(0.);
    {
        const RelaxedPointCloud saved(m_source, bounds, site_spacing, 0,
            [&](double f) { report(0.1 * f); });
        m_bounds = *saved.extent();
        if (m_bounds.defined)
            m_sites = saved.points_in(BoundingBoxf3(m_bounds.min - Vec3d::Constant(EPSILON),
                                                   m_bounds.max + Vec3d::Constant(EPSILON)), site_spacing);
    }
    const auto values = samples(m_sites);
    if (values.size() != m_sites.size()) throw std::invalid_argument("Invalid Voronoi void samples");
    // Interpolation can place a constant field a few ulps below its cutoff.
    // Composited material densities are floats; keep points on the cutoff despite rounding.
    const double tolerance = 4. * std::numeric_limits<float>::epsilon() * std::max(1., threshold);
    size_t kept = 0;
    for (size_t i = 0; i < m_sites.size(); ++i)
        if (values[i].distance >= 0. || values[i].density + tolerance >= threshold) m_sites[kept++] = m_sites[i];
    m_result.removed = m_sites.size() - kept;
    m_sites.resize(kept);
    Mesh mesh;
    for (size_t i = 0; i < m_sites.size(); ++i) {
        mesh.insert(point(m_sites[i]));
        if (i % 4096 == 0) report(0.1);
    }
    const auto evaluate = [&] {
        return bad_faces(mesh, bounds, samples, settings, [&] { report(0.1 + 0.9 * double(m_result.inserted) / std::max(size_t(1), settings.maximum_insertions)); });
    };
    auto faces = bounds.defined ? evaluate() : std::vector<BadFace>{};
    double current = cost(faces);
    m_result.initial_cost = current;
    const double clearance = std::min(site_spacing / 4., settings.maximum_span / 10.);
    while (!faces.empty() && m_result.inserted < settings.maximum_insertions) {
        Vec3d best = Vec3d::Zero();
        double best_cost = current;
        // Consider the worst few faces; a locally blocked face must not prevent
        // an independent part of the constellation from improving.
        for (size_t i = 0; i < std::min(size_t(4), faces.size()); ++i) {
            const auto &face = faces[i];
            std::vector<Vec3d> candidates {face.center,
                face.center - 0.25 * face.span * face.direction,
                face.center + 0.25 * face.span * face.direction,
                face.center - Vec3d(0., 0., 0.25 * face.span)};
            const auto positions = samples(candidates);
            if (positions.size() != candidates.size()) throw std::invalid_argument("Invalid Voronoi refinement candidates");
            for (size_t j = 0; j < candidates.size(); ++j) {
                const Vec3d &candidate = candidates[j];
                if (positions[j].distance >= 0.) continue;
                const auto nearest = mesh.nearest_vertex(point(candidate));
                if ((point(nearest->point()) - candidate).norm() < clearance) continue;
                auto inserted = mesh.insert(point(candidate));
                const double trial = cost(evaluate());
                mesh.remove(inserted);
                if (trial < best_cost - 1e-6 * std::max(1., current)) {
                    best_cost = trial;
                    best = candidate;
                }
            }
        }
        if (best_cost == current) break;
        mesh.insert(point(best));
        ++m_result.inserted;
        faces = evaluate();
        current = cost(faces);
        report(0.1 + 0.9 * double(m_result.inserted) / settings.maximum_insertions);
    }
    m_result.remaining_faces = faces.size();
    m_result.final_cost = current;
    m_sites.clear();
    m_bounds = BoundingBoxf3();
    for (auto v = mesh.finite_vertices_begin(); v != mesh.finite_vertices_end(); ++v) {
        m_sites.push_back(point(v->point()));
        m_bounds.merge(m_sites.back());
    }
    std::sort(m_sites.begin(), m_sites.end(), [](const Vec3d &a, const Vec3d &b) {
        if (a.z() != b.z()) return a.z() < b.z();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.x() < b.x();
    });
    report(1.);
}

std::vector<Vec3d> VoidPointCloud::points_in(const BoundingBoxf3 &box, double) const
{
    std::vector<Vec3d> result;
    if (!box.defined) return result;
    auto first = std::lower_bound(m_sites.begin(), m_sites.end(), box.min.z(), [](const Vec3d &p, double z) { return p.z() < z; });
    for (; first != m_sites.end() && first->z() < box.max.z(); ++first)
        if ((first->array() >= box.min.array()).all() && (first->array() < box.max.array()).all()) result.push_back(*first);
    return result;
}
} // namespace Slic3r::Voronoi
