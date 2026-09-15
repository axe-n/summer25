#ifndef BLACKRAY_RAYTRACER_HPP
#define BLACKRAY_RAYTRACER_HPP

#include "diffeqs.hpp"
#include "isco.hpp"

#include <vector>

namespace blackray {

enum class RayStatus {
    DiskHit,
    Horizon,
    Escaped,
    IntegrationFailure
};

struct CameraTetrad {
    // e[coordinate][local Lorentz index].
    double e[4][4]{};
};

struct Camera {
    double r_obs = 1.0e3;
    double theta_obs = 1.5707963267948966;
    double phi_obs = 0.0;
    double focal_length = 1.0;
};

struct RayTraceOptions {
    double r_out = 1.0e3;
    double escape_radius = 1.0e6;
    double maximum_affine_parameter = 1.0e7;
    int maximum_steps = 1000000;
    double disk_crossing_tolerance = 1.0e-10;
    IntegratorOptions integrator{};
};

struct RayHit {
    double x_obs = 0.0;
    double y_obs = 0.0;
    double r_disk = 0.0;
    double g_factor = 0.0;
    double cos_theta_e = 0.0;
    RayStatus status = RayStatus::IntegrationFailure;
};

CameraTetrad make_camera_tetrad(const Camera& camera, double a,
                                const DeformationParams& params);

// Maps a screen point to a past-directed null coordinate four-vector.
GeodesicState initial_ray(double x_obs, double y_obs, const Camera& camera,
                          const CameraTetrad& tetrad);

RayHit trace_ray(double x_obs, double y_obs, const Camera& camera,
                 double a, const DeformationParams& params,
                 double r_isco, const RayTraceOptions& options = {});

// Appends one result per screen coordinate to output, preserving input order.
void trace_rays(const std::vector<double>& x_obs,
                const std::vector<double>& y_obs,
                const Camera& camera, double a,
                const DeformationParams& params, double r_isco,
                std::vector<RayHit>& output,
                const RayTraceOptions& options = {});

} // namespace blackray

#endif // BLACKRAY_RAYTRACER_HPP