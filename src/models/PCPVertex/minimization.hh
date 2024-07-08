#ifndef UTOPIA_MODELS_PCPVERTEX_MINIMIZATION_HH
#define UTOPIA_MODELS_PCPVERTEX_MINIMIZATION_HH


namespace Utopia {
namespace Models {
namespace PCPVertex {

/// The collection of parameters used for energy minimization
struct MinimizationParams {
    /// The tolerance in energy change 
    double tolerance;

    /// The default timestep
    double dt;

    /// The maximum number of steps per minimization
    std::size_t min_steps;

    /// The maximum number of steps per minimization
    std::size_t max_steps;

    /// Iterate for a fixed number of steps
    /** If num_steps = 0, energy minimized for tolerance
     */
    std::size_t num_steps;

    /// Linetension fluctuation parameter
    /** Parameter fluctuations are implemented as Ornstein-Uhlenbeck process
     * 
     *  \f$  \frac{d\Lambda_{mn}}{dt} = - \frac{1}{\tau_\Lambda}
     *      (\Lambda_{mn}(t) - \Lambda_0)
     *      + \Delta \Lambda \sqrt{2 / \tau_\Lambda} \Theta_{mn}(t)
     *  \f$
     * 
     *  with the first value \f$ \tau \f$ and the second value
     *  \f$ \Delta \Lambda \f$.
     */
    std::pair<double, double> linetension_fluctuations;

    /// The number of jiggling the vertices
    /** \details The first jiggle is applied before the first minimization,
     *           then the energy is minimized up to `jiggle_tolerance`.
     *           This is repeated up to the final minimization, which is 
     *           done up `tolerance`.
     * 
     *  \note Must be > 0
     */
    std::size_t num_repeat;

    /// The reduced tolerance for the preliminary minimizations
    /** \details Only the final minimization is done with `tolerance`.
     *  \note   Optional parameter, default is `tolerance`.
     *  \note   The value of `jiggle_tolerance` cannot be smaller than
     *          `tolerance`.
     */
    double jiggle_tolerance;

    /// The intensity of the jiggling
    /** \details The intensity scales relative to the typical lengthscale of a
     *           cell \f$ l = \sqrt{A_{domain} / num_cells} \f$
     */
    double jiggle_intensity;

    template <typename Config>
    MinimizationParams(const Config& cfg)
    :
        tolerance(get_as<double>("tolerance", cfg)),
        dt(get_as<double>("dt", cfg)),
        min_steps(get_as<std::size_t>("min_steps", cfg, 0)),
        max_steps(get_as<std::size_t>("max_steps", cfg)),
        num_steps(get_as<std::size_t>("num_steps", cfg, 0)),
        linetension_fluctuations(
            std::make_pair(
                get_as<double>("linetension_fluctuations_tau", cfg, 1.),
                get_as<double>("linetension_fluctuations", cfg, 0.)
            )
        ),
        num_repeat(get_as<std::size_t>("num_repeat", cfg, 1)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg, tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg, 0.))
    {
        if (num_repeat == 0) {
            throw Utopia::KeyError("num_repeat", cfg, fmt::format(
                "Value must be larger than 0, but was {}", num_repeat));
        }
        if (num_repeat > 1 and jiggle_tolerance < tolerance) {
            throw Utopia::KeyError("jiggle_tolerance", cfg, fmt::format(
                "Value must be larger or equal to 'tolerance', but was {} < {}",
                jiggle_tolerance, tolerance));
        }

        if (min_steps > max_steps) {
            throw Utopia::KeyError("min_steps", cfg, fmt::format(
                "Value must be smaller or equal `max_steps`, but was {} > {}.",
                min_steps, max_steps));
        }
    }

    /// Initialize from config and inherit not defined values from default
    template <typename Config>
    MinimizationParams(const Config& cfg, const MinimizationParams& defaults)
    :
        tolerance(get_as<double>("tolerance", cfg, defaults.tolerance)),
        dt(get_as<double>("dt", cfg, defaults.dt)),
        min_steps(get_as<std::size_t>("min_steps", cfg, defaults.min_steps)),
        max_steps(get_as<std::size_t>("max_steps", cfg, defaults.max_steps)),
        num_steps(get_as<std::size_t>("num_steps", cfg, defaults.num_steps)),
        linetension_fluctuations(
            std::make_pair(
                get_as<double>("linetension_fluctuations_tau", cfg,
                    std::get<0>(defaults.linetension_fluctuations)),
                get_as<double>("linetension_fluctuations", cfg,
                    std::get<1>(defaults.linetension_fluctuations))
            )
        ),
        num_repeat(get_as<std::size_t>("num_repeat", cfg, defaults.num_repeat)),
        jiggle_tolerance(get_as<double>("jiggle_tolerance", cfg,
                                        defaults.jiggle_tolerance)),
        jiggle_intensity(get_as<double>("jiggle_intensity", cfg,
                                        defaults.jiggle_intensity))
    {        
        if (num_repeat == 0) {
            throw Utopia::KeyError("num_repeat", cfg, fmt::format(
                "Value must be larger than 0, but was {}", num_repeat));
        }
        if (jiggle_tolerance < tolerance) {
            throw Utopia::KeyError("jiggle_tolerance", cfg, fmt::format(
                "Value must be larger or equal to 'tolerance', but was {} < {}",
                jiggle_tolerance, tolerance));
        }

        if (min_steps > max_steps) {
            throw Utopia::KeyError("min_steps", cfg, fmt::format(
                "Value must be smaller or equal `max_steps`, but was {} > {}.",
                min_steps, max_steps));
        }
    }
};


} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif