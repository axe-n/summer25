#include "metric.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace blackray {
namespace {

using Real = long double;

struct Components {
    Real gtt;
    Real grr;
    Real gthth;
    Real gtphi;
    Real gphiphi;
};

void validate_inputs(Real r, Real theta, Real a,
                     const DeformationParams& params) {
    if (!std::isfinite(r) || !std::isfinite(theta) || !std::isfinite(a) ||
        r <= 0.0L || !std::isfinite(params.alpha13) ||
        !std::isfinite(params.alpha22) || !std::isfinite(params.alpha52) ||
        !std::isfinite(params.epsilon3)) {
        throw std::invalid_argument("metric inputs must be finite and r > 0");
    }
}

Components calculate_components(Real r, Real theta, Real a,
                                const DeformationParams& params) {
    validate_inputs(r, theta, a, params);

    const Real sin_theta = std::sin(theta);
    const Real cos_theta = std::cos(theta);
    const Real sin2 = sin_theta * sin_theta;
    const Real cos2 = cos_theta * cos_theta;
    const Real r2 = r * r;
    const Real a2 = a * a;
    const Real delta = r2 - 2.0L * r + a2;
    const Real a1 = 1.0L + static_cast<Real>(params.alpha13) / (r2 * r);
    const Real a2_func = 1.0L + static_cast<Real>(params.alpha22) / r2;
    const Real a5 = 1.0L + static_cast<Real>(params.alpha52) / r2;
    const Real f = static_cast<Real>(params.epsilon3) / r;
    const Real sigma = r2 + a2 * cos2 + f;
    const Real b = (r2 + a2) * a1 - a2 * a2_func * sin2;

    if (b == 0.0L) {
        throw std::domain_error("Johannsen metric is singular where B = 0");
    }

    const Real b2 = b * b;
    Components result{
        -sigma * (delta - a2 * a2_func * a2_func * sin2) / b2,
        sigma / (delta * a5),
        sigma,
        -a * sigma * sin2 * ((r2 + a2) * a1 * a2_func - delta) / b2,
        sigma * sin2 * ((r2 + a2) * (r2 + a2) * a1 * a1 -
                        a2 * delta * sin2) / b2};

    if (!std::isfinite(result.gtt) || !std::isfinite(result.grr) ||
        !std::isfinite(result.gthth) || !std::isfinite(result.gtphi) ||
        !std::isfinite(result.gphiphi)) {
        throw std::domain_error("Johannsen metric is singular at this point");
    }
    return result;
}

void store_compact(const Components& components, double g[3]) {
    g[G_TT] = static_cast<double>(components.gtt);
    g[G_TPHI] = static_cast<double>(components.gtphi);
    g[G_PHIPHI] = static_cast<double>(components.gphiphi);
}

template <typename Function>
Real five_point_derivative(Function&& function, Real x, Real scale) {
    const Real h = 1.0e-4L * std::max(1.0L, scale);
    return (function(x - 2.0L * h) - 8.0L * function(x - h) +
            8.0L * function(x + h) - function(x + 2.0L * h)) /
           (12.0L * h);
}

} // namespace

void get_metric_tensor(double r, double theta, double a,
                       const DeformationParams& params, double g[3]) {
    if (g == nullptr) {
        throw std::invalid_argument("metric output pointer must not be null");
    }
    store_compact(calculate_components(r, theta, a, params), g);
}

void get_inverse_metric(const double g[3], double g_inv[3]) {
    if (g == nullptr || g_inv == nullptr) {
        throw std::invalid_argument("metric pointers must not be null");
    }
    const double determinant = g[G_TT] * g[G_PHIPHI] -
                              g[G_TPHI] * g[G_TPHI];
    if (!std::isfinite(determinant) || determinant == 0.0) {
        throw std::domain_error("metric t-phi block is not invertible");
    }
    g_inv[G_TT] = g[G_PHIPHI] / determinant;
    g_inv[G_TPHI] = -g[G_TPHI] / determinant;
    g_inv[G_PHIPHI] = g[G_TT] / determinant;
}

void get_metric_derivatives(double r, double theta, double a,
                            const DeformationParams& params,
                            double dg_dr[3], double dg_dtheta[3]) {
    if (dg_dr == nullptr || dg_dtheta == nullptr) {
        throw std::invalid_argument("derivative output pointers must not be null");
    }
    validate_inputs(r, theta, a, params);

    const Real rr = r;
    const Real tt = theta;
    for (int component = 0; component < 3; ++component) {
        const auto radial = [=](Real value) {
            double compact[3];
            get_metric_tensor(static_cast<double>(value), theta, a, params,
                              compact);
            return static_cast<Real>(compact[component]);
        };
        const auto polar = [=](Real value) {
            double compact[3];
            get_metric_tensor(r, static_cast<double>(value), a, params,
                              compact);
            return static_cast<Real>(compact[component]);
        };
        dg_dr[component] = static_cast<double>(five_point_derivative(
            radial, rr, std::abs(rr)));
        dg_dtheta[component] = static_cast<double>(five_point_derivative(
            polar, tt, 1.0L));
    }
}

void get_full_metric_tensor(double r, double theta, double a,
                            const DeformationParams& params,
                            double g[4][4]) {
    if (g == nullptr) {
        throw std::invalid_argument("metric output pointer must not be null");
    }
    const Components components = calculate_components(r, theta, a, params);
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            g[row][column] = 0.0;
        }
    }
    g[0][0] = static_cast<double>(components.gtt);
    g[0][3] = static_cast<double>(components.gtphi);
    g[1][1] = static_cast<double>(components.grr);
    g[2][2] = static_cast<double>(components.gthth);
    g[3][0] = g[0][3];
    g[3][3] = static_cast<double>(components.gphiphi);
}

} // namespace blackray