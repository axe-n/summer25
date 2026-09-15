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
                           double observer_radius, double inclination, double spin,
                           double alpha13, double alpha22, double alpha52,
                           double epsilon3, double step, int steps,
                           double* maximum_norm);

int blackray_metric_diagnostics(double r, double theta, double spin,
                                double alpha13, double alpha22,
                                double alpha52, double epsilon3,
                                double* determinant, double* g_phiphi);

/* Finds the ISCO for the given Johannsen deformation parameters.
   All output pointers are optional (may be NULL); only r_isco is required. */
int blackray_find_isco(double spin, double alpha13, double alpha22,
                       double alpha52, double epsilon3,
                       double inner_radius, double outer_radius,
                       double radial_scan_step, double root_tolerance,
                       int maximum_iterations,
                       double* r_isco, double* omega, double* u_t,
                       double* energy, double* angular_momentum,
                       double* radial_curvature, double* vertical_curvature,
                       int* vertical_instability_warning);

/* Convolves ray samples with a local spectrum table.
   angle_mode: 0 = RayTraced, 1 = FixedInclination
   Arrays are allocated by the caller; sizes must match. */
int blackray_convolve_spectrum(
    const double* r_disk, const double* g_factor, const double* cos_theta_e,
    int ray_count,
    const double* spec_energies, const double* spec_cos_theta,
    const double* spec_flux, int n_energies, int n_angles,
    double emissivity_index, double fixed_cos_inclination, int angle_mode,
    int output_bins, double output_energy_min, double output_energy_max,
    double* energy_lo, double* energy_hi, double* flux);

#ifdef __cplusplus
}
#endif

#endif /* BLACKRAY_C_API_H */
