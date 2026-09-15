#include "blackray_c_api.h"

#include "convolver.hpp"
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

extern "C" int blackray_find_isco(double spin, double alpha13, double alpha22,
                                  double alpha52, double epsilon3,
                                  double inner_radius, double outer_radius,
                                  double radial_scan_step, double root_tolerance,
                                  int maximum_iterations,
                                  double* r_isco, double* omega, double* u_t,
                                  double* energy, double* angular_momentum,
                                  double* radial_curvature,
                                  double* vertical_curvature,
                                  int* vertical_instability_warning) {
    try {
        if (r_isco == nullptr) return 1;
        blackray::DeformationParams params{alpha13, alpha22, alpha52, epsilon3};
        blackray::IscoOptions options;
        options.inner_radius = inner_radius;
        options.outer_radius = outer_radius;
        options.radial_scan_step = radial_scan_step;
        options.root_tolerance = root_tolerance;
        options.maximum_iterations = maximum_iterations;
        const blackray::IscoResult result =
            blackray::find_isco(spin, params, options);
        *r_isco = result.radius;
        if (omega) *omega = result.omega;
        if (u_t) *u_t = result.u_t;
        if (energy) *energy = result.energy;
        if (angular_momentum) *angular_momentum = result.angular_momentum;
        if (radial_curvature) *radial_curvature = result.radial_potential_curvature;
        if (vertical_curvature) *vertical_curvature = result.vertical_potential_curvature;
        if (vertical_instability_warning)
            *vertical_instability_warning = result.vertical_instability_warning ? 1 : 0;
        return 0;
    } catch (const std::exception&) {
        return 2;
    }
}

extern "C" int blackray_convolve_spectrum(
    const double* r_disk, const double* g_factor, const double* cos_theta_e,
    int ray_count,
    const double* spec_energies, const double* spec_cos_theta,
    const double* spec_flux, int n_energies, int n_angles,
    double emissivity_index, double fixed_cos_inclination, int angle_mode,
    int output_bins, double output_energy_min, double output_energy_max,
    double* energy_lo, double* energy_hi, double* flux) {
    try {
        if (r_disk == nullptr || g_factor == nullptr || cos_theta_e == nullptr ||
            ray_count < 0 || spec_energies == nullptr || spec_cos_theta == nullptr ||
            spec_flux == nullptr || n_energies < 1 || n_angles < 1 ||
            output_bins < 1 || energy_lo == nullptr || energy_hi == nullptr ||
            flux == nullptr) return 1;

        std::vector<blackray::RaySample> rays(ray_count);
        for (int i = 0; i < ray_count; ++i) {
            rays[i].r_disk = r_disk[i];
            rays[i].g_factor = g_factor[i];
            rays[i].cos_theta_e = cos_theta_e[i];
        }

        blackray::LocalSpectrumTable table;
        table.energy.assign(spec_energies, spec_energies + n_energies);
        table.cos_inclination.assign(spec_cos_theta, spec_cos_theta + n_angles);
        table.spectrum.assign(spec_flux,
            spec_flux + static_cast<std::size_t>(n_energies) * n_angles);

        blackray::ConvolutionOptions options;
        options.emissivity_index = emissivity_index;
        options.fixed_cos_inclination = fixed_cos_inclination;
        options.angle_mode = (angle_mode == 1)
            ? blackray::EmissionAngleMode::FixedInclination
            : blackray::EmissionAngleMode::RayTraced;
        options.output_bins = output_bins;
        options.output_energy_min = output_energy_min;
        options.output_energy_max = output_energy_max;

        const blackray::ConvolvedSpectrum result =
            blackray::convolve_spectrum(rays, table, options);

        if (static_cast<int>(result.flux.size()) != output_bins) return 2;
        for (int i = 0; i < output_bins; ++i) {
            energy_lo[i] = result.energy_lo[i];
            energy_hi[i] = result.energy_hi[i];
            flux[i] = result.flux[i];
        }
        return 0;
    } catch (const std::exception&) {
        return 2;
    }
}