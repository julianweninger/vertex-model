#ifndef UTOPIA_MODELS_PCPVERTEX_SPACE_HH
#define UTOPIA_MODELS_PCPVERTEX_SPACE_HH

#include <cmath>
#include <limits>

#include <armadillo>

#include <utopia/core/space.hh>
#include <utopia/data_io/cfg_utils.hh>
#include <utopia/core/types.hh>

namespace Utopia::Models::PCPVertex::Space {

template<std::size_t num_dims>
struct CustomSpace : public Utopia::Space<num_dims> {
public:
    /// The dimensionality of the space
    using Space = Utopia::Space<num_dims>;

    /// The type of a vector in space
    using SpaceVec = typename Space::SpaceVec;

    /// whether the space is bound by the extent
    const bool bound;

private:
    /// The scaling of the domain
    /** A remapping vector to a space extend in absolute coordinates.
     * 
     *  \note other than the extent this can be subject to Topological changes
     */
    SpaceVec _domain_scale;

    /// Lees-Edwards boundary condition, i.e. skew
    SpaceVec _skew;

public:
    CustomSpace(const DataIO::Config& cfg)
    :
        Space(cfg),
        bound(this->periodic),
        _domain_scale(arma::fill::ones),
        _skew(arma::fill::zeros)
    {
        if (bound != this->periodic) {
            throw std::invalid_argument("The space can only be unbound in "
                "non-periodic boundary conditions! If the space is non-periodic"
                "the Vertex model relies on unbound space configuration.");
        }
    }

    /// Constructor without any arguments, i.e. constructing a default space
    /** \details The default space is non-periodic and has default extent of 1.
      *         into each dimension.
      */
    CustomSpace()
    :
        Space(),
        bound(this->periodic),
        _domain_scale(arma::fill::ones),
        _skew(arma::fill::zeros)
    {
        if (bound != this->periodic) {
            throw std::invalid_argument("The space can only be unbound in "
                "non-periodic boundary conditions! If the space is non-periodic"
                "the Vertex model relies on unbound space configuration.");
        }
    }

    /// Whether this space contains the given coordinate (without mapping it)
    /** \note   No distinction is made between periodic and non-periodic space.
      *
      * \tparam  include_high_value_boundary  Whether to check the closed or
      *         the half-open interval. The latter case is useful when working
      *         with periodic grids, allowing to map values on the high-value
      *         boundary back to the low-value boundary.
      */
    template<bool include_high_value_boundary=true>
    bool contains(const SpaceVec& pos) const {
        if (not bound) {
            return true;
        }
        else {
            return Space::contains(pos / _domain_scale);
        }
    }

    /// Map a position (potentially outside space's extent) back into space
    /** \details This is intended for use with periodic grids. It will also work
      *         with non-periodic grids, but the input value should not have
      *         been permitted in the first place.
      *
      * \note   The high-value boundary is is mapped back to the low-value
      *         boundary, such that all points are well-defined.
      */
    SpaceVec map_into_space(const SpaceVec& pos) const {
        if (not bound) {
            return pos;
        }
        else if (arma::norm(_skew) > 1.e-12){
            SpaceVec domain = this->get_domain_size();

            SpaceVec add_skew = arma::floor(
                arma::shift(pos / domain, -1) ) % _skew;

            return _domain_scale % Space::map_into_space(
                (pos - add_skew) / _domain_scale
            );
        }
        else {
            return Space::map_into_space(pos / _domain_scale) % _domain_scale;
        }
    }
    
    /// The displacement between 2 coordinates
    /** \details Calculates vector pointing from pos_0 to pos_1.
     *           In periodic boundary it calculates the shorter displacement.
     *  
     *  \warning The displacement of two coordinates in periodic boundary can be
     *           maximum of half the domain size, i.e. moving away from a
     *           coordinate in a certain direction will decrease the
     *           displacement once reached half the domain size.
     */
    SpaceVec displacement(const SpaceVec& pos_0, const SpaceVec& pos_1) const {
        if (not this->periodic) {
            return _domain_scale % Space::displacement(
                pos_0 / _domain_scale,
                pos_1 / _domain_scale
            );
        }
        else if (arma::norm(_skew) > 1.e-12) {
            // the virtual position (might be outside the domain)
            SpaceVec displ = pos_1 - pos_0;

            SpaceVec domain = get_domain_size();

            SpaceVec add_skew = (  arma::shift(arma::round(displ / domain), -1)
                                 % (-1. * _skew));

            return _domain_scale % Space::displacement(
                pos_0 / _domain_scale,
                (pos_1 + add_skew) / _domain_scale
            );
        }
        else {
            return _domain_scale % Space::displacement(
                pos_0 / _domain_scale,
                pos_1 / _domain_scale
            );
        }
    }

    /// The distance of 2 coordinates in space
    /** \details Calculates the distance of 2 coordinates using the norm
     *           implemented within Armadillo.
     *           In periodic boundary it calculates the shorter distance.
     * 
     *  \warning The distance of two coordinates in periodic boundary can be
     *           maximum of half the domain size wrt every dimension,
     *           i.e. moving away from a coordinate in a certain direction will
     *           decrease the distance once reached half the domain size.
     * 
     *  \param   p   The norm used to compute the distance, see arma::norm(X, p).
     *               Can be either an integer >= 1 or one of "-inf", "inf", "fro"
     */
    template<class NormType=std::size_t>
    auto distance(const SpaceVec& pos_0, const SpaceVec& pos_1,
                  const NormType p=2) const
    {
        return arma::norm(displacement(pos_0, pos_1), p);
    }

    
    /// Intersection of two (finite) lines
    /** \param pos_0    The origin of first line
     *  \param vec_0    The direction of first line
     *  \param pos_1    The origin of second line
     *  \param vec_1    The direction of second line
     *  \param finite_0 Whether the first line is finite, 
     *                  i.e. from pos_0 to pos_0 + vec_0, or infinit
     *  \param finite_1 Whether the second line is finite, 
     *                  i.e. from pos_1 to pos_1 + vec_1, or infinit
     * 
     *  \return intersection, success   The coordinate of the intersection
     *                                  and whether the lines have an
     *                                  intersection
     */
    template <typename std::enable_if_t<num_dims == 2, int> = 0>
    std::pair<SpaceVec, bool> intersection(SpaceVec pos_0, SpaceVec vec_0,
            SpaceVec pos_1, SpaceVec vec_1,
            bool finite_0=true, bool finite_1=true) const
    {
        // work with a copy of pos_1 relative to pos_0
        pos_1 = pos_0 + displacement(pos_0, pos_1);

        // define line as y0 = a + b*x
        double a, b;
        if (fabs(vec_0.at(0)) < 1e-14) { vec_0 += SpaceVec({2.e-14, 0}); }
        b = vec_0.at(1) / vec_0.at(0);
        a = pos_0.at(1) - b * pos_0.at(0);

        // define line as y1 = c + d*x
        double c, d;
        if (fabs(vec_1.at(0)) < 1e-14) { vec_1 += SpaceVec({2.e-14, 0}); }
        d = vec_1.at(1) / vec_1.at(0);
        c = pos_1.at(1) - d * pos_1.at(0);

        double x = (c-a) / (b-d);
        SpaceVec intersection = {x, a + b*x};

        // parallel lines
        if (fabs(d - b) < 1e-14) {
            intersection = this->map_into_space(intersection);
            return std::make_pair(intersection, false);
        }
                
        // consider finite length
        if (finite_0) {
            double x0, x1;
            x0 = pos_0.at(0);
            x1 = x0 + vec_0.at(0);

            if (x0 <= x1 and (intersection.at(0) < x0 - 1e-14 or 
                              intersection.at(0) > x1 + 1e-14)) {
                intersection = this->map_into_space(intersection);
                return std::make_pair(intersection, false);
            }
            if (x0 > x1 and (intersection.at(0)  > x0 + 1e-14 or 
                             intersection.at(0)  < x1 - 1e-14)) {
                intersection = this->map_into_space(intersection);
                return std::make_pair(intersection, false);
            }
        }
        if (finite_1) {
            double x0, x1;
            x0 = pos_1.at(0);
            x1 = x0 + vec_1.at(0);

            if (x0 <= x1 and (intersection.at(0) < x0 - 1e-14 or 
                              intersection.at(0) > x1 + 1e-14)) {
                intersection = this->map_into_space(intersection);
                return std::make_pair(intersection, false);
            }
            if (x0 > x1 and (intersection.at(0)  > x0 + 1e-14 or 
                             intersection.at(0)  < x1 - 1e-14)) {
                intersection = this->map_into_space(intersection);
                return std::make_pair(intersection, false);
            }
        }

        intersection = this->map_into_space(intersection);
        return std::make_pair(intersection, true);
    }

    /// Getter for the size of the domain
    SpaceVec get_domain_size() const {
        return this->extent % _domain_scale;
    }

    /// The d-dimensional volume of the domain
    double get_domain_volume() const {
        return arma::prod(get_domain_size());
    }

    /// Getter for the scale of the domain
    /** \details positions in space are remapped to stretches of the domain size
     */
    void get_domain_scale() const {
       return _domain_scale;
    }

    /// Setter for the scale of the domain
    /** \details positions in space are remapped to stretches of the domain size
     */
    void set_domain_size(const SpaceVec& domain_size) {
        _domain_scale = domain_size / this->extent;
    }

    /// Getter for skew boundary condition
    const SpaceVec& get_skew() const {
        return _skew;
    }

    /// Setter for skew boundary condition
    void set_skew(const SpaceVec& skew) {
        if (not this->periodic and arma::norm(skew) > 1.e-12) {
            throw std::runtime_error(fmt::format(
                "Cannot set non-zero skew ({}, {}) in non-periodic boundary "
                "conditions!", skew[0], skew[1]
            ));
        }
        if (fabs(skew[0]) > 1.e-12 and fabs(skew[1]) > 1.e-12)
        {
            throw std::runtime_error(fmt::format(
                "Cannot set skew in x ({}) and y ({}) using Lee-Edwards periodic "
                "boundary condition", skew[0], skew[1]
            ));
        }

        _skew = skew;
    }
}; // struct CustomSpace

} // namespace Utopia::Models::PCPVertex::Space
#endif
