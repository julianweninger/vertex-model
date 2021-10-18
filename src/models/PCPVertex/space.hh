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
    /** \details Checks whether the given coordinate is within this space's
      *         extent by computing the relative position and checking whether
      *         it is within [0, 1] or [0, 1) for all elements.
      *
      * \note   No distinction is made between periodic and non-periodic space.
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
            return Space::contains(this->map_to_relative_space(pos));
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
            SpaceVec _pos = this->map_to_absolute_space(
                Space::map_into_space(this->map_to_relative_space(pos)));
            SpaceVec __pos = _pos;
            if (_pos[1] < pos[1] - 1.e-10) {
                // crossed upper boundary
                __pos[0] -= _skew[0];
            }
            else if (_pos[1] > pos[1] + 1.e-10) {
                // crossed lower boundary
                __pos[0] += _skew[0];
            }
            if (_pos[0] < pos[0] - 1.e-10) {
                // crossed left boundary
                __pos[1] -= _skew[1];
            }
            else if (_pos[0] > pos[0] + 1.e-10) {
                // crossed right boundary
                __pos[1] += _skew[1];
            }

            return this->map_to_absolute_space(
                Space::map_into_space(this->map_to_relative_space(__pos)));
        }
        else {
            return this->map_to_absolute_space(Space::map_into_space(
                    this->map_to_relative_space(pos)));
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
        SpaceVec displ = this->map_to_absolute_space(Space::displacement(
                    this->map_to_relative_space(pos_0),
                    this->map_to_relative_space(pos_1)));

        if (arma::norm(_skew) > 1.e-12) {
            // the virtual position (might be outside the domain)
            SpaceVec _pos_1 = pos_0 + displ;

            if (_pos_1[0] > get_domain_size()[0]) {
                // right boundary
                displ += _skew % SpaceVec({0., 1.});
            }
            else if (_pos_1[0] < 0) {
                // left boundary
                displ -= _skew % SpaceVec({0., 1.});
            }

            if (_pos_1[1] > get_domain_size()[1]) {
                // upper boundary
                displ += _skew % SpaceVec({1., 0.});
            }
            else if (_pos_1[1] < 0) {
                // lower boundary
                displ -= _skew % SpaceVec({1., 0.});
            }
        }

        return displ;
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
    
    /// Map the position in to absolute coordinates
    SpaceVec map_to_absolute_space(const SpaceVec& pos) const {
        return pos % _domain_scale;
    }
    
    
    /// Map the position in to relative coordinates
    SpaceVec map_to_relative_space(const SpaceVec& abs_pos) const {
        return abs_pos / _domain_scale;
    }  

    /// Getter for the size of the domain
    SpaceVec get_domain_size() const {
        return this->extent % _domain_scale;
    }

    /// Setter for the scale of the domain
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
    SpaceVec get_skew() const {
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
        SpaceVec domain = get_domain_size();

        if (skew[0] > domain[0] / 3. or skew[1] > domain[1] / 3.) {
            throw std::runtime_error(fmt::format(
                "Cannot set skew to ({}, {}) as it exceeds 1/3 of the domain "
                "size ({}, {})!", skew[0], skew[1], domain[0], domain[1]
            ));
        }

        _skew = skew;
    }
}; // struct CustomSpace

} // namespace Utopia::Models::PCPVertex::Space
#endif
