#include "isco.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace blackray {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct OrbitConstants {
    double energy;
    double angular_momentum;
};

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
            throw std::domain_error("metric tensor is not invertible");
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

void validate_options(const IscoOptions& options) {
    if (!(options.outer_radius > options.inner_radius) ||
        !(options.radial_scan_step > 0.0) || !(options.root_tolerance > 0.0) ||
        options.maximum_iterations < 1) {
        throw std::invalid_argument("invalid ISCO solver options");
    }
}

double effective_potential(double r, double theta, double a,
                           const DeformationParams& params,
                           const OrbitConstants& constants) {
    double metric[4][4], inverse[4][4];
    get_full_metric_tensor(r, theta, a, params, metric);
    invert_metric(metric, inverse);
    return 1.0 + inverse[0][0] * constants.energy * constants.energy -
           2.0 * inverse[0][3] * constants.energy *
               constants.angular_momentum +
           inverse[3][3] * constants.angular_momentum *
               constants.angular_momentum;
}

double radial_curvature(double r, double a, const DeformationParams& params,
                        const OrbitConstants& constants) {
    const double h = 1.0e-4 * std::max(1.0, std::abs(r));
    return (effective_potential(r - h, kPi / 2.0, a, params, constants) -
            2.0 * effective_potential(r, kPi / 2.0, a, params, constants) +
            effective_potential(r + h, kPi / 2.0, a, params, constants)) /
           (h * h);
}

double vertical_curvature(double r, double a, const DeformationParams& params,
                          const OrbitConstants& constants) {
    const double h = 1.0e-4;
    return (effective_potential(r, kPi / 2.0 - h, a, params, constants) -
            2.0 * effective_potential(r, kPi / 2.0, a, params, constants) +
            effective_potential(r, kPi / 2.0 + h, a, params, constants)) /
           (h * h);
}

OrbitConstants circular_constants(double r, double omega, double a,
                                  const DeformationParams& params) {
    double metric[4][4];
    get_full_metric_tensor(r, kPi / 2.0, a, params, metric);
    const double normalization = -metric[0][0] -
        2.0 * omega * metric[0][3] - omega * omega * metric[3][3];
    if (!(normalization > 0.0) || !std::isfinite(normalization)) {
        throw std::domain_error("circular orbit is not timelike");
    }
    const double u_t = 1.0 / std::sqrt(normalization);
    return {-u_t * (metric[0][0] + omega * metric[0][3]),
            u_t * (metric[0][3] + omega * metric[3][3])};
}

} // namespace

double keplerian_omega(double r, double a, const DeformationParams& params) {
    double dg_dr[3], dg_dtheta[3];
    get_metric_derivatives(r, kPi / 2.0, a, params, dg_dr, dg_dtheta);
    const double discriminant = dg_dr[G_TPHI] * dg_dr[G_TPHI] -
        dg_dr[G_TT] * dg_dr[G_PHIPHI];
    if (!(discriminant >= 0.0) || !std::isfinite(discriminant) ||
        dg_dr[G_PHIPHI] == 0.0) {
        throw std::domain_error("no real Keplerian angular velocity");
    }
    return (-dg_dr[G_TPHI] + std::sqrt(discriminant)) /
           dg_dr[G_PHIPHI];
}

double circular_u_t(double r, double omega, double a,
                    const DeformationParams& params) {
    double metric[4][4];
    get_full_metric_tensor(r, kPi / 2.0, a, params, metric);
    const double normalization = -metric[0][0] - 2.0 * omega * metric[0][3] -
                                 omega * omega * metric[3][3];
    if (!(normalization > 0.0) || !std::isfinite(normalization)) {
        throw std::domain_error("circular orbit is not timelike");
    }
    return 1.0 / std::sqrt(normalization);
}

IscoResult find_isco(double a, const DeformationParams& params,
                     const IscoOptions& options) {
    validate_options(options);
    const double horizon = 1.0 + std::sqrt(std::max(0.0, 1.0 - a * a));
    const double inner = std::max({options.inner_radius, horizon + 1.0e-3,
                                   1.0e-3});
    if (!(options.outer_radius > inner)) {
        throw std::invalid_argument("ISCO search interval is empty");
    }

    bool have_previous = false;
    double previous_radius = inner;
    double previous_curvature = 0.0;
    for (double radius = inner; radius <= options.outer_radius;
         radius += options.radial_scan_step) {
        try {
            const double omega = keplerian_omega(radius, a, params);
            const OrbitConstants constants = circular_constants(radius, omega, a, params);
            const double curvature = radial_curvature(radius, a, params, constants);
            if (have_previous && std::isfinite(curvature) &&
                previous_curvature * curvature <= 0.0) {
                double low = previous_radius;
                double high = radius;
                for (int iteration = 0; iteration < options.maximum_iterations;
                     ++iteration) {
                    const double midpoint = 0.5 * (low + high);
                    const double midpoint_omega = keplerian_omega(midpoint, a, params);
                    const OrbitConstants midpoint_constants =
                        circular_constants(midpoint, midpoint_omega, a, params);
                    const double midpoint_curvature =
                        radial_curvature(midpoint, a, params, midpoint_constants);
                    if (std::abs(high - low) <= options.root_tolerance) {
                        low = high = midpoint;
                        break;
                    }
                    if (previous_curvature * midpoint_curvature <= 0.0) {
                        high = midpoint;
                    } else {
                        low = midpoint;
                        previous_curvature = midpoint_curvature;
                    }
                }
                const double result_radius = 0.5 * (low + high);
                const double result_omega = keplerian_omega(result_radius, a, params);
                const OrbitConstants result_constants =
                    circular_constants(result_radius, result_omega, a, params);
                const double vertical = vertical_curvature(
                    result_radius, a, params, result_constants);
                return {result_radius, result_omega,
                        circular_u_t(result_radius, result_omega, a, params),
                        result_constants.energy, result_constants.angular_momentum,
                        radial_curvature(result_radius, a, params, result_constants),
                        vertical, vertical <= 0.0};
            }
            previous_radius = radius;
            previous_curvature = curvature;
            have_previous = true;
        } catch (const std::domain_error&) {
            have_previous = false;
        }
    }
    throw std::runtime_error("ISCO root was not found in the requested interval");
}

} // namespace blackray