#ifndef BLACKRAY_ISCO_HPP
#define BLACKRAY_ISCO_HPP

#include "metric.hpp"

namespace blackray {

struct IscoOptions {
    double inner_radius = 0.0;
    double outer_radius = 50.0;
    double radial_scan_step = 0.02;
    double root_tolerance = 1.0e-10;
    int maximum_iterations = 100;
};

struct IscoResult {
    double radius = 0.0;
    double omega = 0.0;
    double u_t = 0.0;
    double energy = 0.0;
    double angular_momentum = 0.0;
    double radial_potential_curvature = 0.0;
    double vertical_potential_curvature = 0.0;
    bool vertical_instability_warning = false;
};

// Computes the prograde (+ square-root branch) equatorial Keplerian angular
// velocity from the radial metric derivatives.
double keplerian_omega(double r, double a, const DeformationParams& params);

// Computes u^t for a circular equatorial orbit at radius r.
double circular_u_t(double r, double omega, double a,
                    const DeformationParams& params);

// Finds the first radial marginal-stability root outside the horizon-like
// inner boundary and reports a possible vertical orbital instability.
IscoResult find_isco(double a, const DeformationParams& params,
                     const IscoOptions& options = {});

} // namespace blackray

#endif // BLACKRAY_ISCO_HPP