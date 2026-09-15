#ifndef BLACKRAY_C_API_H
#define BLACKRAY_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

int blackray_trace_grid(const double* x_obs, const double* y_obs, int count,
                        double observer_radius, double inclination,
                        double spin, double alpha13, double alpha22,
                        double alpha52, double epsilon3, double r_isco,
                        double r_out, double* r_disk, double* g_factor,
                        double* cos_theta_e, int* status);

            int blackray_max_null_norm(const double* x_obs, const double* y_obs, int count,
                           double observer_radius, double inclination,
                           double spin, double alpha13, double alpha22,
                           double alpha52, double epsilon3, double step,
                           int steps, double* maximum_norm);

            int blackray_metric_diagnostics(double r, double theta, double spin,
                            double alpha13, double alpha22,
                            double alpha52, double epsilon3,
                            double* determinant, double* g_phiphi);

#ifdef __cplusplus
}
#endif

#endif /* BLACKRAY_C_API_H */