#include "convolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace blackray {
namespace {

double interpolate(const std::vector<double>& x, const std::vector<double>& y,
                   double value) {
    if (x.size() != y.size() || x.empty()) {
        throw std::invalid_argument("spectrum axes and values must be nonempty");
    }
    if (value <= x.front()) return y.front();
    if (value >= x.back()) return y.back();
    const auto upper = std::upper_bound(x.begin(), x.end(), value);
    const std::size_t hi = static_cast<std::size_t>(upper - x.begin());
    const std::size_t lo = hi - 1;
    const double fraction = (value - x[lo]) / (x[hi] - x[lo]);
    return y[lo] + fraction * (y[hi] - y[lo]);
}

double local_flux(const LocalSpectrumTable& table, double energy,
                  double cosine) {
    if (table.cos_inclination.empty() || table.energy.empty() ||
        table.spectrum.size() != table.energy.size() *
                                  table.cos_inclination.size()) {
        throw std::invalid_argument("invalid local spectrum table");
    }
    const double clamped = std::clamp(std::abs(cosine),
                                      table.cos_inclination.front(),
                                      table.cos_inclination.back());
    std::vector<double> angle_values(table.energy.size());
    for (std::size_t angle = 0; angle < table.cos_inclination.size(); ++angle) {
        std::vector<double> row(table.energy.size());
        std::copy_n(table.spectrum.begin() + angle * table.energy.size(),
                    table.energy.size(), row.begin());
        angle_values[angle] = interpolate(table.energy, row, energy);
    }
    return interpolate(table.cos_inclination, angle_values, clamped);
}

} // namespace

ConvolvedSpectrum convolve_spectrum(const std::vector<RaySample>& rays,
                                    const LocalSpectrumTable& table,
                                    const ConvolutionOptions& options) {
    if (options.output_bins < 1 || !(options.output_energy_min > 0.0) ||
        !(options.output_energy_max > options.output_energy_min)) {
        throw std::invalid_argument("invalid convolution energy grid");
    }
    ConvolvedSpectrum result;
    result.energy_lo.resize(options.output_bins);
    result.energy_hi.resize(options.output_bins);
    result.flux.assign(options.output_bins, 0.0);
    const double log_min = std::log(options.output_energy_min);
    const double log_width = (std::log(options.output_energy_max) - log_min) /
                             options.output_bins;
    for (int bin = 0; bin < options.output_bins; ++bin) {
        result.energy_lo[bin] = std::exp(log_min + bin * log_width);
        result.energy_hi[bin] = std::exp(log_min + (bin + 1) * log_width);
    }
    for (const RaySample& ray : rays) {
        if (!(ray.r_disk > 0.0) || !(ray.g_factor > 0.0) ||
            !std::isfinite(ray.r_disk) || !std::isfinite(ray.g_factor)) {
            continue;
        }
        const double cosine = options.angle_mode == EmissionAngleMode::RayTraced
            ? ray.cos_theta_e : options.fixed_cos_inclination;
        const double radial_weight = std::pow(ray.r_disk, options.emissivity_index);
        for (int bin = 0; bin < options.output_bins; ++bin) {
            const double observed_energy = 0.5 *
                (result.energy_lo[bin] + result.energy_hi[bin]);
            const double local_energy = observed_energy / ray.g_factor;
            const double local = local_flux(table, local_energy, cosine);
            result.flux[bin] += radial_weight * std::pow(ray.g_factor, 3.0) * local;
        }
    }
    return result;
}

} // namespace blackray