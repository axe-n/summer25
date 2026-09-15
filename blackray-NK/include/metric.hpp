#ifndef BLACKRAY_METRIC_HPP
#define BLACKRAY_METRIC_HPP

namespace blackray {

struct DeformationParams {
    double alpha13 = 0.0;
    double alpha22 = 0.0;
    double alpha52 = 0.0;
    double epsilon3 = 0.0;
};

// Compact ordering for the three independent t-phi block components.
enum MetricComponent : int {
    G_TT = 0,
    G_TPHI = 1,
    G_PHIPHI = 2
};

// Computes [g_tt, g_tphi, g_phiphi] in Boyer-Lindquist coordinates.
void get_metric_tensor(double r, double theta, double a,
                       const DeformationParams& params, double g[3]);

// Inverts the compact t-phi block returned by get_metric_tensor.
void get_inverse_metric(const double g[3], double g_inv[3]);

// Computes radial and polar derivatives in the same compact ordering.
void get_metric_derivatives(double r, double theta, double a,
                            const DeformationParams& params,
                            double dg_dr[3], double dg_dtheta[3]);

// Computes the complete 4D metric when radial and polar components are needed.
void get_full_metric_tensor(double r, double theta, double a,
                            const DeformationParams& params,
                            double g[4][4]);

} // namespace blackray

#endif // BLACKRAY_METRIC_HPP