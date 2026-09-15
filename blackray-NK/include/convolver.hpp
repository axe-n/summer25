#ifndef BLACKRAY_CONVOLVER_HPP
#define BLACKRAY_CONVOLVER_HPP

#include <vector>

namespace blackray {

enum class EmissionAngleMode {
    RayTraced,
    FixedInclination
};

struct RaySample {
    double r_disk;
    double g_factor;
    double cos_theta_e;
};

struct LocalSpectrumTable {
    std::vector<double> energy;
    std::vector<double> cos_inclination;
    // Row-major: spectrum[angle_index * energy.size() + energy_index].
    std::vector<double> spectrum;
};

struct ConvolutionOptions {
    double fixed_cos_inclination = 0.5;
    double emissivity_index = -3.0;
    EmissionAngleMode angle_mode = EmissionAngleMode::RayTraced;
    int output_bins = 2000;
    double output_energy_min = 0.00035;
    double output_energy_max = 2000.0;
};

struct ConvolvedSpectrum {
    std::vector<double> energy_lo;
    std::vector<double> energy_hi;
    std::vector<double> flux;
};

ConvolvedSpectrum convolve_spectrum(const std::vector<RaySample>& rays,
                                    const LocalSpectrumTable& table,
                                    const ConvolutionOptions& options = {});

} // namespace blackray

#endif // BLACKRAY_CONVOLVER_HPP