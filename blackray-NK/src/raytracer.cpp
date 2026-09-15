#include "raytracer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace blackray {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void invert_metric(const double metric[4][4], double inverse[4][4]) {
    double augmented[4][8]{};
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            augmented[row][column] = metric[row][column];
            augmented[row][column + 4] = row == column ? 1.0 : 0.0;
        }
    }
    for (int column = 0; column < 4; ++column) {
        int pivot = column;
        for (int row = column + 1; row < 4; ++row) {
            if (std::abs(augmented[row][column]) >
                std::abs(augmented[pivot][column])) {
                pivot = row;
            }
        }
        if (!std::isfinite(augmented[pivot][column]) ||
            std::abs(augmented[pivot][column]) <=
                std::numeric_limits<double>::epsilon()) {
            throw std::domain_error("observer metric is not invertible");
        }
        for (int entry = 0; entry < 8; ++entry) {
            std::swap(augmented[column][entry], augmented[pivot][entry]);
        }
        const double scale = augmented[column][column];
        for (double& entry : augmented[column]) {
            entry /= scale;
        }
        for (int row = 0; row < 4; ++row) {
            if (row == column) {
                continue;
            }
            const double factor = augmented[row][column];
            for (int entry = 0; entry < 8; ++entry) {
                augmented[row][entry] -= factor * augmented[column][entry];
            }
        }
    }
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            inverse[row][column] = augmented[row][column + 4];
        }
    }
}

void validate_camera(const Camera& camera) {
    if (!(camera.r_obs > 0.0) || !(camera.focal_length > 0.0) ||
        !std::isfinite(camera.r_obs) || !std::isfinite(camera.theta_obs) ||
        !std::isfinite(camera.phi_obs) || !std::isfinite(camera.focal_length)) {
        throw std::invalid_argument("invalid camera parameters");
    }
}

void validate_options(const RayTraceOptions& options) {
    if (!(options.r_out > 0.0) || !(options.escape_radius > options.r_out) ||
        !(options.maximum_affine_parameter > 0.0) || options.maximum_steps < 1 ||
        !(options.disk_crossing_tolerance > 0.0)) {
        throw std::invalid_argument("invalid ray-tracing options");
    }
}

double covariant_component(const double metric[4][4], const GeodesicState& state,
                           int component) {
    double result = 0.0;
    for (int coordinate = 0; coordinate < 4; ++coordinate) {
        result += metric[component][coordinate] * state[coordinate + 4];
    }
    return result;
}

RayHit disk_observation(const GeodesicState& state, double x_obs, double y_obs,
                        double a, const DeformationParams& params) {
    double metric[4][4];
    get_full_metric_tensor(state[1], kPi / 2.0, a, params, metric);
    const double k_t = covariant_component(metric, state, 0);
    const double k_phi = covariant_component(metric, state, 3);
    const double omega = keplerian_omega(state[1], a, params);
    const double u_t = circular_u_t(state[1], omega, a, params);
    const double denominator = u_t * (1.0 - omega * (-k_phi / k_t));
    if (k_t == 0.0 || !std::isfinite(denominator) || denominator == 0.0) {
        throw std::domain_error("invalid disk redshift denominator");
    }
    const double emitter_energy = -(k_t + omega * k_phi);
    const double normal_projection =
        std::sqrt(metric[2][2]) * state[6];
    return {x_obs, y_obs, state[1], 1.0 / denominator,
            std::abs(normal_projection / emitter_energy), RayStatus::DiskHit};
}

} // namespace

CameraTetrad make_camera_tetrad(const Camera& camera, double a,
                                const DeformationParams& params) {
    validate_camera(camera);
    double metric[4][4], inverse[4][4];
    get_full_metric_tensor(camera.r_obs, camera.theta_obs, a, params, metric);
    invert_metric(metric, inverse);
    const double lapse_inverse = std::sqrt(-inverse[0][0]);
    if (!std::isfinite(lapse_inverse) || lapse_inverse <= 0.0 ||
        metric[1][1] <= 0.0 || metric[2][2] <= 0.0 || metric[3][3] <= 0.0) {
        throw std::domain_error("observer does not admit a spatial tetrad");
    }

    CameraTetrad tetrad;
    tetrad.e[0][0] = lapse_inverse;
    tetrad.e[3][0] = -inverse[0][3] / lapse_inverse;
    tetrad.e[1][1] = 1.0 / std::sqrt(metric[1][1]);
    tetrad.e[2][2] = 1.0 / std::sqrt(metric[2][2]);
    tetrad.e[3][3] = 1.0 / std::sqrt(metric[3][3]);
    return tetrad;
}

GeodesicState initial_ray(double x_obs, double y_obs, const Camera& camera,
                          const CameraTetrad& tetrad) {
    validate_camera(camera);
    if (!std::isfinite(x_obs) || !std::isfinite(y_obs)) {
        throw std::invalid_argument("screen coordinates must be finite");
    }
    const double norm = std::sqrt(x_obs * x_obs + y_obs * y_obs +
                                  camera.focal_length * camera.focal_length);
    const double nx = x_obs / norm;
    const double ny = y_obs / norm;
    const double nz = camera.focal_length / norm;
    GeodesicState state;
    state[0] = 0.0;
    state[1] = camera.r_obs;
    state[2] = camera.theta_obs;
    state[3] = camera.phi_obs;
    // Local k^(0) = -1 makes the ray past-directed; spatial direction points
    // inward through the negative radial tetrad leg.
    const double local[4] = {-1.0, -nz, ny, nx};
    for (int coordinate = 0; coordinate < 4; ++coordinate) {
        state[coordinate + 4] = 0.0;
        for (int leg = 0; leg < 4; ++leg) {
            state[coordinate + 4] += tetrad.e[coordinate][leg] * local[leg];
        }
    }
    return state;
}

RayHit trace_ray(double x_obs, double y_obs, const Camera& camera,
                 double a, const DeformationParams& params,
                 double r_isco, const RayTraceOptions& options) {
    validate_camera(camera);
    validate_options(options);
    if (!(r_isco > 0.0) || !(options.r_out >= r_isco)) {
        throw std::invalid_argument("invalid disk radial limits");
    }
    const CameraTetrad tetrad = make_camera_tetrad(camera, a, params);
    GeodesicState state = initial_ray(x_obs, y_obs, camera, tetrad);
    const double horizon = 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a));
    double step = options.integrator.maximum_step;
    double affine = 0.0;
    for (int iteration = 0; iteration < options.maximum_steps &&
         affine < options.maximum_affine_parameter; ++iteration) {
        const GeodesicState previous = state;
        const double previous_theta = previous[2];
        if (!rk4_step(state, step, a, params, options.integrator)) {
            return {x_obs, y_obs, 0.0, 0.0, 0.0, RayStatus::IntegrationFailure};
        }
        affine += std::abs(step);
        const double theta_product =
            (previous_theta - kPi / 2.0) * (state[2] - kPi / 2.0);
        if (theta_product <= 0.0 &&
            std::abs(state[2] - previous_theta) > options.disk_crossing_tolerance) {
            const double fraction = (kPi / 2.0 - previous_theta) /
                                    (state[2] - previous_theta);
            GeodesicState crossing = previous;
            for (int index = 0; index < kStateSize; ++index) {
                crossing[index] += fraction * (state[index] - previous[index]);
            }
            if (crossing[1] >= r_isco && crossing[1] <= options.r_out) {
                return disk_observation(crossing, x_obs, y_obs, a, params);
            }
        }
        if (state[1] <= horizon || !std::isfinite(state[1])) {
            return {x_obs, y_obs, 0.0, 0.0, 0.0, RayStatus::Horizon};
        }
        if (state[1] >= options.escape_radius ||
            !std::isfinite(state[2]) || !std::isfinite(state[3])) {
            return {x_obs, y_obs, 0.0, 0.0, 0.0, RayStatus::Escaped};
        }
    }
    return {x_obs, y_obs, 0.0, 0.0, 0.0, RayStatus::IntegrationFailure};
}

void trace_rays(const std::vector<double>& x_obs,
                const std::vector<double>& y_obs,
                const Camera& camera, double a,
                const DeformationParams& params, double r_isco,
                std::vector<RayHit>& output,
                const RayTraceOptions& options) {
    if (x_obs.size() != y_obs.size()) {
        throw std::invalid_argument("screen coordinate buffers must have equal size");
    }
    output.clear();
    output.reserve(x_obs.size());
    for (std::size_t index = 0; index < x_obs.size(); ++index) {
        output.push_back(trace_ray(x_obs[index], y_obs[index], camera, a,
                                   params, r_isco, options));
    }
}

} // namespace blackray