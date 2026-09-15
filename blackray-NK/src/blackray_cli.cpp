#include "isco.hpp"
#include "raytracer.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Config {
    int rstep = 128;
    int pstep = 128;
    double screen_radius = 20.0;
    double observer_radius = 1.0e3;
    double inclination = 1.0;
    double spin = 0.5;
    double r_out = 50.0;
    std::string output = "ray_hits.csv";
    std::string config_file;
};

double parse_double(const char* value, const char* name) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0') {
        throw std::invalid_argument(std::string("invalid value for ") + name);
    }
    return parsed;
}

int parse_int(const char* value, const char* name) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < 1) {
        throw std::invalid_argument(std::string("invalid value for ") + name);
    }
    return static_cast<int>(parsed);
}

void print_usage(const char* executable) {
    std::cout << "Usage: " << executable << " [options]\n"
              << "  --rstep N             radial screen samples (default 128)\n"
              << "  --pstep N             azimuthal screen samples (default 128)\n"
              << "  --screen-radius R     screen half-width (default 20)\n"
              << "  --observer-radius R   observer radius (default 1000)\n"
              << "  --inclination I       observer theta (default 1)\n"
              << "  --spin A              dimensionless spin (default 0.5)\n"
              << "  --r-out R             outer disk radius (default 50)\n"
              << "  --output FILE         CSV output path (default ray_hits.csv)\n"
              << "  --config FILE         key=value configuration file\n"
              << "  --help                show this help\n";
}

void apply_option(Config& config, const std::string& option,
                  const std::string& value) {
    if (option == "rstep" || option == "--rstep") config.rstep = parse_int(value.c_str(), "rstep");
    else if (option == "pstep" || option == "--pstep") config.pstep = parse_int(value.c_str(), "pstep");
    else if (option == "screen-radius" || option == "--screen-radius") config.screen_radius = parse_double(value.c_str(), "screen-radius");
    else if (option == "observer-radius" || option == "--observer-radius") config.observer_radius = parse_double(value.c_str(), "observer-radius");
    else if (option == "inclination" || option == "--inclination") config.inclination = parse_double(value.c_str(), "inclination");
    else if (option == "spin" || option == "--spin") config.spin = parse_double(value.c_str(), "spin");
    else if (option == "r-out" || option == "--r-out") config.r_out = parse_double(value.c_str(), "r-out");
    else if (option == "output" || option == "--output") config.output = value;
    else throw std::invalid_argument("unknown configuration key: " + option);
}

void read_config(Config& config, const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open config file: " + path);
    std::string line;
    while (std::getline(input, line)) {
        const std::size_t comment = line.find('#');
        line = line.substr(0, comment);
        const std::size_t separator = line.find('=');
        if (separator == std::string::npos) {
            if (!line.empty()) throw std::invalid_argument("config entries must use key=value");
            continue;
        }
        std::string key = line.substr(0, separator);
        std::string value = line.substr(separator + 1);
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t\r") + 1);
        if (!key.empty() && !value.empty()) apply_option(config, key, value);
    }
}

Config parse_arguments(int argc, char** argv) {
    Config config;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--help") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        }
        if (index + 1 >= argc) throw std::invalid_argument("missing value for " + option);
        const std::string value = argv[++index];
        if (option == "--config") {
            config.config_file = value;
            read_config(config, value);
        } else {
            apply_option(config, option, value);
        }
    }
    if (!(config.screen_radius > 0.0) || !(config.observer_radius > 0.0) ||
        !(config.r_out > 0.0) || std::abs(config.spin) >= 1.0 ||
        config.output.empty()) {
        throw std::invalid_argument("invalid CLI configuration");
    }
    return config;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_arguments(argc, argv);
        blackray::DeformationParams params;
        blackray::IscoOptions isco_options;
        isco_options.outer_radius = config.r_out;
        const blackray::IscoResult isco =
            blackray::find_isco(config.spin, params, isco_options);

        blackray::Camera camera;
        camera.r_obs = config.observer_radius;
        camera.theta_obs = config.inclination;
        camera.focal_length = config.observer_radius;

        blackray::RayTraceOptions trace_options;
        trace_options.r_out = config.r_out;
        trace_options.escape_radius = config.observer_radius * 2.0;
        trace_options.integrator.maximum_step = 1.0;

        const std::size_t sample_count =
            static_cast<std::size_t>(config.rstep) * config.pstep;
        std::vector<blackray::RayHit> hits(sample_count);
        const double radial_step = 2.0 * config.screen_radius / config.rstep;
        const double azimuthal_step = 2.0 * config.screen_radius / config.pstep;

#pragma omp parallel for schedule(dynamic)
        for (int linear = 0; linear < static_cast<int>(sample_count); ++linear) {
            const int radial_index = linear / config.pstep;
            const int azimuthal_index = linear % config.pstep;
            const double x = -config.screen_radius +
                (radial_index + 0.5) * radial_step;
            const double y = -config.screen_radius +
                (azimuthal_index + 0.5) * azimuthal_step;
            hits[linear] = blackray::trace_ray(x, y, camera, config.spin,
                                               params, isco.radius, trace_options);
        }

        std::ofstream output(config.output);
        if (!output) throw std::runtime_error("cannot open output file");
        output << "x_obs,y_obs,r_disk,g_factor,cos_theta_e,status\n";
        output << std::setprecision(17);
        for (const auto& hit : hits) {
            output << hit.x_obs << ',' << hit.y_obs << ',' << hit.r_disk << ','
                   << hit.g_factor << ',' << hit.cos_theta_e << ','
                   << static_cast<int>(hit.status) << '\n';
        }
        std::cerr << "ISCO radius: " << isco.radius << "\n"
                  << "Wrote " << hits.size() << " rays to " << config.output << "\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "blackray_cli: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}