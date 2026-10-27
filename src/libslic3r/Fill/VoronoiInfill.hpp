#pragma once

#include "../BoundingBox.hpp"
#include "../ExPolygon.hpp"
#include "../Polyline.hpp"
#include "libslic3r/Point.hpp"
#include "VoronoiInfillCgal.hpp"
#include "VoronoiPointCloud.hpp"
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <cstddef>

namespace Slic3r {

// Printable sections of a 3D point field, with replaceable site sampling.
class VoronoiInfill
{
public:
    explicit VoronoiInfill(std::shared_ptr<const Voronoi::PointCloudProvider> provider = std::make_shared<Voronoi::PoissonPointCloud>())
        : m_provider(std::move(provider)), m_cloud(m_provider) {}

    struct GeometryParams {
        double extrusion_spacing; // Nominal single-line spacing used to calibrate point density.
        double density;
        double z;
        double width;
        double height;
        double reference_height = -1.;
        double row_spacing = -1.; // Actual flow spacing for floor rows.
        int wall_count = 1;
        double point_spacing = 0.; // Explicit cloud spacing overrides density calibration.
    };

    // Clip wall sections and floor bands to the fill region, before extrusion.
    Voronoi::Layer fill_geometry(const BoundingBox &object_bounds, const ExPolygon &region,
                                const GeometryParams &params);
    Polylines route_walls(Polylines walls, size_t layer_id, const Points &boundary_nodes = {}) const;

private:
    std::shared_ptr<const Voronoi::PointCloudProvider> m_provider;
    Voronoi::ActivePointCloud m_cloud;
    BoundingBox m_bounds;
    double      m_site_spacing { 0. };
    double      m_layer_z { std::numeric_limits<double>::quiet_NaN() };
    double      m_layer_width { 0. };
    double      m_layer_height { 0. };
    double      m_row_spacing { 0. };
    int m_wall_count { 0 };
    Voronoi::CellGeometry m_geometry;
    size_t m_geometry_revision = std::numeric_limits<size_t>::max();
    Voronoi::Section m_section;
    std::vector<Voronoi::Band> m_bands;
};

// One registry per ascending worker batch. Successive fillers reuse engines
// sequentially; concurrent workers own independent registries and networks.
struct VoronoiEngineRegistry {
    std::shared_ptr<VoronoiInfill> get(std::shared_ptr<const Voronoi::PointCloudProvider> provider, double spacing)
    {
        auto &engine = m_engines[{provider, spacing}];
        if (!engine) engine = provider ? std::make_shared<VoronoiInfill>(std::move(provider)) : std::make_shared<VoronoiInfill>();
        return engine;
    }
private:
    std::map<std::pair<std::shared_ptr<const Voronoi::PointCloudProvider>, double>, std::shared_ptr<VoronoiInfill>> m_engines;
};

} // namespace Slic3r
