#ifndef UTOPIA_MODELS_PCPVERTEX_UTILS_HH
#define UTOPIA_MODELS_PCPVERTEX_UTILS_HH

namespace Utopia {
namespace Models {
namespace PCPVertex {

std::lognormal_distribution<double> get_lognormal_distribution
(double mean, double stddev)
{
    if (mean < 0.) {
        throw std::invalid_argument(fmt::format("Cannot construct "
            "lognormal distribution with mean {} < 0!", mean));
    }
    if (stddev <= 0.) {
        throw std::invalid_argument(fmt::format("Cannot construct "
            "lognormal distribution with stddev {} <= 0!", stddev));
    }

    double m = log(  std::pow(mean, 2)
                   / sqrt(std::pow(stddev, 2) + std::pow(mean, 2)));
    double s = sqrt(log(std::pow(stddev, 2) / std::pow(mean, 2) + 1));

   return std::lognormal_distribution<double>(m, s);
}

} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif