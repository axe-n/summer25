#ifndef BLACKRAY_DIFFEQS_HPP
#define BLACKRAY_DIFFEQS_HPP

#include "metric.hpp"

namespace blackray {

constexpr int kStateSize = 8;

struct GeodesicState {
    double value[kStateSize]{};

    double& operator[](int index) { return value[index]; }
    const double& operator[](int index) const { return value[index]; }
};

struct IntegratorOptions {
    double absolute_tolerance = 1.0e-10;
    double relative_tolerance = 1.0e-8;
    double null_tolerance = 1.0e-9;
    double minimum_step = 1.0e-12;
    double maximum_step = 1.0;
    int maximum_rejections = 32;
};

// Computes dY/dlambda for Y = [t, r, theta, phi, kt, kr, ktheta, kphi].
void geodesic_derivatives(const GeodesicState& state, double a,
                          const DeformationParams& params,
                          GeodesicState& derivative);

// Computes Gamma^mu_(alpha beta) using the metric and its coordinate
// derivatives. The first two coordinate derivatives are zero by stationarity
// and axisymmetry; the radial and polar derivatives are calculated internally.
void christoffel_symbols(double r, double theta, double a,
                         const DeformationParams& params,
                         double gamma[4][4][4]);

// Takes one accepted adaptive RK4 step. On return, state is updated and step
// contains the next suggested step size. Returns false when the requested step
// cannot satisfy the tolerances before minimum_step.
bool rk4_step(GeodesicState& state, double& step, double a,
              const DeformationParams& params,
              const IntegratorOptions& options = {});

double null_norm(const GeodesicState& state, double a,
                 const DeformationParams& params);

} // namespace blackray

#endif // BLACKRAY_DIFFEQS_HPP