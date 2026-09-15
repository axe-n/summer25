#include "blackray_c_api.h"

#include "isco.hpp"
#include "raytracer.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>

extern "C" int blackray_trace_grid(const double* x_obs, const double* y_obs,
                                    int count, double observer_radius,
                                    double inclination, double spin,
                                    double alpha13, double alpha22,
                                    double alpha52, double epsilon3,
                                    double r_isco, double r_out,
                                    double* r_disk, double* g_factor,
                                    double* cos_theta_e, int* status) {
    try {
        if (x_obs == nullptr || y_obs == nullptr || r_disk == nullptr ||
            g_factor == nullptr || cos_theta_e == nullptr || status == nullptr ||
            count < 0) return 1;
        blackray::DeformationParams params{alpha13, alpha22, alpha52, epsilon3};
        if (!(r_isco > 0.0)) {
            r_isco = blackray::find_isco(spin, params).radius;
        }
        blackray::Camera camera;
        camera.r_obs = observer_radius;
        camera.theta_obs = inclination;
        camera.focal_length = observer_radius;
        blackray::RayTraceOptions options;
        options.r_out = r_out;
        options.escape_radius = std::max(2.0 * observer_radius, r_out * 2.0);
#pragma omp parallel for schedule(dynamic)
        for (int index = 0; index < count; ++index) {
            try {
                const blackray::RayHit hit = blackray::trace_ray(
                    x_obs[index], y_obs[index], camera, spin, params,
                    r_isco, options);
                r_disk[index] = hit.r_disk;
                g_factor[index] = hit.g_factor;
                cos_theta_e[index] = hit.cos_theta_e;
                status[index] = static_cast<int>(hit.status);
            } catch (const std::exception&) {
                r_disk[index] = 0.0;
                g_factor[index] = 0.0;
                cos_theta_e[index] = 0.0;
                status[index] = static_cast<int>(blackray::RayStatus::IntegrationFailure);
            }
        }
        return 0;
    } catch (const std::exception&) {
        return 2;
    }
}

extern "C" int blackray_max_null_norm(
    const double* x_obs, const double* y_obs, int count,
    double observer_radius, double inclination, double spin,
    double alpha13, double alpha22, double alpha52, double epsilon3,
    double step, int steps, double* maximum_norm) {
    try {
        if (x_obs == nullptr || y_obs == nullptr || maximum_norm == nullptr ||
            count < 0 || steps < 0 || !std::isfinite(step)) return 1;
        blackray::DeformationParams params{alpha13, alpha22, alpha52, epsilon3};
        blackray::Camera camera;
        camera.r_obs = observer_radius;
        camera.theta_obs = inclination;
        camera.focal_length = observer_radius;
        const blackray::CameraTetrad tetrad =
            blackray::make_camera_tetrad(camera, spin, params);
        blackray::IntegratorOptions options;
        options.null_tolerance = 1.0e-7;
        options.maximum_step = std::abs(step);
        double maximum = 0.0;
        for (int ray = 0; ray < count; ++ray) {
            blackray::GeodesicState state =
                blackray::initial_ray(x_obs[ray], y_obs[ray], camera, tetrad);
            maximum = std::max(maximum,
                               std::abs(blackray::null_norm(state, spin, params)));
            double current_step = std::abs(step);
            for (int iteration = 0; iteration < steps; ++iteration) {
                if (!blackray::rk4_step(state, current_step, spin, params, options)) {
                    return 3;
                }
                maximum = std::max(maximum,
                                   std::abs(blackray::null_norm(state, spin, params)));
            }
        }
        *maximum_norm = maximum;
        return 0;
    } catch (const std::exception&) {
        return 2;
    }
}

extern "C" int blackray_metric_diagnostics(
    double r, double theta, double spin, double alpha13, double alpha22,
    double alpha52, double epsilon3, double* determinant, double* g_phiphi) {
    try {
        if (determinant == nullptr || g_phiphi == nullptr) return 1;
        blackray::DeformationParams params{alpha13, alpha22, alpha52, epsilon3};
        double metric[4][4];
        blackray::get_full_metric_tensor(r, theta, spin, params, metric);
        *g_phiphi = metric[3][3];
        *determinant = metric[0][0] * metric[1][1] * metric[2][2] * metric[3][3] -
                       metric[0][3] * metric[0][3] * metric[1][1] * metric[2][2];
        return 0;
    } catch (const std::exception&) {
        return 2;
    }
}