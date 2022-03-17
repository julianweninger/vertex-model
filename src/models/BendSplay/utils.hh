#ifndef UTOPIA_MODELS_BENDSPLAY_UTILS_HH
#define UTOPIA_MODELS_BENDSPLAY_UTILS_HH

#include <cmath>
#include <random>
#include <numeric>
#include <vector>
#include <limits>

#include <armadillo>
#include <spdlog/spdlog.h>




namespace Utopia::Models::BendSplay {


// -- Angle-related tools -----------------------------------------------------

/// Returns a uniformly random angle value in [-π, +π)
template<class RNG, class T=double>
T random_angle (const std::shared_ptr<RNG>& rng) {
    return std::uniform_real_distribution<T>(-M_PI, +M_PI)(*rng);
}

/// Constrains an angle value to interval [-π, +π)
template<class T>
T constrain_angle (T angle) {
    angle = std::fmod(angle + M_PI, 2 * M_PI);
    while (angle < 0.) {
        angle += 2 * M_PI;
    }
    return angle - M_PI;
}

// /// In-place constrains all angles in a container to interval [-π, +π)
// void constrain_angles (arma::mat& angles) {
//     angles.transform([](double val) { constrain_angle(val); });
// }


// -- Circular Statistics -----------------------------------------------------

/// Compute sum of sine and cosine values from angles in a container
auto _circular_sin_cos_sum (const arma::mat& angles) {
    const auto sin_sum = arma::accu(arma::sin(angles));
    const auto cos_sum = arma::accu(arma::cos(angles));

    return std::make_pair(sin_sum, cos_sum);
}


/// Computes the circular mean from a sample of (constrained) angles
/** Uses circular statistics to compute the mean.
  * Assumes angles to be in radians. While it does not matter in which
  * interval they are, the resulting mean value will be in [-π, +π).
  *
  * Returns NaN if the given container is empty.
  *
  * See scipy implementation for reference:
  *   https://github.com/scipy/scipy/blob/v1.7.1/scipy/stats/morestats.py#L3474
  */
auto circular_mean (const arma::mat& angles) {
    const auto&& [sin_sum, cos_sum] = _circular_sin_cos_sum(angles);
    return constrain_angle(std::atan2(sin_sum, cos_sum));
}


/// Computes the circular mean and std from a sample of (constrained) angles
/** Uses circular statistics to compute the mean and standard deviation.
  * Assumes angles to be in radians. While it does not matter in which
  * interval they are, the resulting mean value will be in [-π, +π).
  *
  * Returns NaN if the given container is empty.
  *
  * See scipy implementation for reference:
  *   https://github.com/scipy/scipy/blob/v1.7.1/scipy/stats/morestats.py#L3595
  */
auto circular_mean_and_std (const arma::mat& angles) {
    const auto&& [sin_sum, cos_sum] = _circular_sin_cos_sum(angles);
    const auto mean = constrain_angle(std::atan2(sin_sum, cos_sum));

    const auto r = std::min(
        1., std::sqrt(sin_sum*sin_sum + cos_sum*cos_sum) / angles.n_elem
    );
    const auto std = std::sqrt(-2. * std::log(r));

    return std::make_pair(mean, std);
}


} // namespace Utopia::Models::Bendsplay

#endif // UTOPIA_MODELS_SIMPLEFLOCKING_UTILS_HH