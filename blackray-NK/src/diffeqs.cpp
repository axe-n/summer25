#include "diffeqs.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace blackray {
namespace {

constexpr int kCoordinateCount = 4;
constexpr int kPositionOffset = 0;
constexpr int kMomentumOffset = 4;

void validate_options(const IntegratorOptions& options) {
    if (!(options.absolute_tolerance > 0.0) ||
        !(options.relative_tolerance > 0.0) ||
        !(options.null_tolerance > 0.0) || !(options.minimum_step > 0.0) ||
        !(options.maximum_step >= options.minimum_step) ||
        options.maximum_rejections < 1) {
        throw std::invalid_argument("invalid geodesic integrator options");
    }
}

void metric_derivatives(double r, double theta, double a,
                        const DeformationParams& params,
                        double derivative[4][4][4]) {
    const double h_r = 1.0e-4 * std::max(1.0, std::abs(r));
    const double h_theta = 1.0e-5;
    for (int coordinate = 0; coordinate < kCoordinateCount; ++coordinate) {
        const double h = coordinate == 1 ? h_r : h_theta;
        if (coordinate < 1 || coordinate > 2) {
            for (int row = 0; row < kCoordinateCount; ++row) {
                for (int column = 0; column < kCoordinateCount; ++column) {
                    derivative[coordinate][row][column] = 0.0;
                }
            }
            continue;
        }
        double plus_one[4][4], plus_two[4][4];
        double minus_one[4][4], minus_two[4][4];
        auto evaluate = [&](double coordinate_value, double output[4][4]) {
            if (coordinate == 1) {
                get_full_metric_tensor(coordinate_value, theta, a, params, output);
            } else {
                get_full_metric_tensor(r, coordinate_value, a, params, output);
            }
        };
        evaluate(coordinate == 1 ? r + h : theta + h, plus_one);
        evaluate(coordinate == 1 ? r + 2.0 * h : theta + 2.0 * h, plus_two);
        evaluate(coordinate == 1 ? r - h : theta - h, minus_one);
        evaluate(coordinate == 1 ? r - 2.0 * h : theta - 2.0 * h, minus_two);
        for (int row = 0; row < kCoordinateCount; ++row) {
            for (int column = 0; column < kCoordinateCount; ++column) {
                derivative[coordinate][row][column] =
                    (minus_two[row][column] - 8.0 * minus_one[row][column] +
                     8.0 * plus_one[row][column] - plus_two[row][column]) /
                    (12.0 * h);
            }
        }
    }
}

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
        for (int entry = 0; entry < 8; ++entry) {
            augmented[column][entry] /= scale;
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

void derivative_at(const GeodesicState& state, double a,
                   const DeformationParams& params,
                   GeodesicState& derivative) {
    double gamma[4][4][4];
    christoffel_symbols(state[1], state[2], a, params, gamma);
    for (int coordinate = 0; coordinate < 4; ++coordinate) {
        derivative[coordinate] = state[coordinate + kMomentumOffset];
        double acceleration = 0.0;
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int beta = 0; beta < 4; ++beta) {
                acceleration += gamma[coordinate][alpha][beta] *
                                state[alpha + kMomentumOffset] *
                                state[beta + kMomentumOffset];
            }
        }
        derivative[coordinate + kMomentumOffset] = -acceleration;
    }
}

void rk4_trial(const GeodesicState& initial, double h, double a,
               const DeformationParams& params, GeodesicState& result) {
    GeodesicState k1, k2, k3, k4, trial;
    derivative_at(initial, a, params, k1);
    for (int index = 0; index < kStateSize; ++index) {
        trial[index] = initial[index] + 0.5 * h * k1[index];
    }
    derivative_at(trial, a, params, k2);
    for (int index = 0; index < kStateSize; ++index) {
        trial[index] = initial[index] + 0.5 * h * k2[index];
    }
    derivative_at(trial, a, params, k3);
    for (int index = 0; index < kStateSize; ++index) {
        trial[index] = initial[index] + h * k3[index];
    }
    derivative_at(trial, a, params, k4);
    for (int index = 0; index < kStateSize; ++index) {
        result[index] = initial[index] + h *
            (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index]) /
            6.0;
    }
}

double error_norm(const GeodesicState& full, const GeodesicState& half,
                  const GeodesicState& initial,
                  const IntegratorOptions& options) {
    double result = 0.0;
    for (int index = 0; index < kStateSize; ++index) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance *
            std::max(std::abs(initial[index]), std::abs(full[index]));
        result = std::max(result, std::abs(full[index] - half[index]) / scale);
    }
    return result;
}

} // namespace

void christoffel_symbols(double r, double theta, double a,
                         const DeformationParams& params,
                         double gamma[4][4][4]) {
    if (gamma == nullptr) {
        throw std::invalid_argument("Christoffel output pointer must not be null");
    }
    double metric[4][4], inverse[4][4], derivative[4][4][4];
    get_full_metric_tensor(r, theta, a, params, metric);
    invert_metric(metric, inverse);
    metric_derivatives(r, theta, a, params, derivative);
    for (int mu = 0; mu < 4; ++mu) {
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int beta = 0; beta < 4; ++beta) {
                gamma[mu][alpha][beta] = 0.0;
                for (int nu = 0; nu < 4; ++nu) {
                    gamma[mu][alpha][beta] += 0.5 * inverse[mu][nu] *
                        (derivative[beta][nu][alpha] +
                         derivative[alpha][nu][beta] -
                         derivative[nu][alpha][beta]);
                }
            }
        }
    }
}

void geodesic_derivatives(const GeodesicState& state, double a,
                          const DeformationParams& params,
                          GeodesicState& derivative) {
    derivative_at(state, a, params, derivative);
}

double null_norm(const GeodesicState& state, double a,
                 const DeformationParams& params) {
    double metric[4][4];
    get_full_metric_tensor(state[1], state[2], a, params, metric);
    double norm = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            norm += metric[mu][nu] * state[mu + kMomentumOffset] *
                    state[nu + kMomentumOffset];
        }
    }
    return norm;
}

bool rk4_step(GeodesicState& state, double& step, double a,
              const DeformationParams& params,
              const IntegratorOptions& options) {
    validate_options(options);
    if (!std::isfinite(step) || step == 0.0) {
        throw std::invalid_argument("integrator step must be finite and nonzero");
    }
    const double direction = step < 0.0 ? -1.0 : 1.0;
    double trial_step = direction * std::min(std::abs(step), options.maximum_step);
    const double initial_norm = null_norm(state, a, params);
    for (int rejection = 0; rejection < options.maximum_rejections; ++rejection) {
        GeodesicState full, midpoint, half;
        rk4_trial(state, trial_step, a, params, full);
        rk4_trial(state, trial_step * 0.5, a, params, midpoint);
        rk4_trial(midpoint, trial_step * 0.5, a, params, half);
        const double error = error_norm(full, half, state, options);
        const double norm_scale = std::max(1.0, std::abs(initial_norm));
        const bool null_ok = std::abs(null_norm(half, a, params) - initial_norm) <=
                             options.null_tolerance * norm_scale;
        if (std::isfinite(error) && error <= 1.0 && null_ok) {
            state = half;
            const double factor = error == 0.0 ? 2.0 :
                std::clamp(0.9 * std::pow(error, -0.2), 0.2, 2.0);
            step = direction * std::min(options.maximum_step,
                                        std::max(options.minimum_step,
                                                 std::abs(trial_step) * factor));
            return true;
        }
        if (std::abs(trial_step) <= options.minimum_step) {
            return false;
        }
        trial_step *= 0.5;
    }
    return false;
}

} // namespace blackray