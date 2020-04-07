#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_SPACE_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_SPACE_HH

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

private:
    /// The scaling of the domain
    /** A remapping vector to a space extend in absolute coordinates.
     * 
     *  \note other than the extent this can be subject to Topological changes
     */
    SpaceVec _domain_scale;

public:
    CustomSpace(const DataIO::Config& cfg)
    :
        Space(cfg),
        _domain_scale(arma::fill::ones)
    { }

    /// Constructor without any arguments, i.e. constructing a default space
    /** \details The default space is non-periodic and has default extent of 1.
      *         into each dimension.
      */
    CustomSpace()
    :
        Space(),
        _domain_scale(arma::fill::ones)
    { }

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
        return Space::contains(this->map_to_relative_space(pos));
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
        return this->map_to_absolute_space(Space::map_into_space(
                this->map_to_relative_space(pos)));
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
        return this->map_to_absolute_space(Space::displacement(
                    this->map_to_relative_space(pos_0),
                    this->map_to_relative_space(pos_1)));
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
    std::pair<SpaceVec, bool> intersection(SpaceVec pos_0, SpaceVec vec_0,
            SpaceVec pos_1, SpaceVec vec_1,
            bool finite_0=true, bool finite_1=true)
    {
        static_assert(this->dim == 2, "Intersection of lines not implemented " 
                      "in other than 2 dimensional space! ");
                      // Space dimension was "+ std::to_string(this->dim)+ ".");
        
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
        intersection = this->map_into_space(intersection);
        
        // parallel lines
        if (fabs(d - b) < 1e-14) {
            return std::make_pair(intersection, false);
        }
                
        // consider finite length
        if (finite_0) {
            double x0, x1;
            x0 = pos_0.at(0);
            x1 = x0 + vec_0.at(0);

            if (x0 <= x1 and (intersection.at(0) < x0 - 1e-14 or 
                              intersection.at(0) > x1 + 1e-14)) {
                return std::make_pair(intersection, false);
            }
            if (x0 > x1 and (intersection.at(0)  > x0 + 1e-14 or 
                             intersection.at(0)  < x1 - 1e-14)) {
                return std::make_pair(intersection, false);
            }
        }
        if (finite_1) {
            double x0, x1;
            x0 = pos_1.at(0);
            x1 = x0 + vec_1.at(0);

            if (x0 <= x1 and (intersection.at(0) < x0 - 1e-14 or 
                              intersection.at(0) > x1 + 1e-14)) {
                return std::make_pair(intersection, false);
            }
            if (x0 > x1 and (intersection.at(0)  > x0 + 1e-14 or 
                             intersection.at(0)  < x1 - 1e-14)) {
                return std::make_pair(intersection, false);
            }
        }

        return std::make_pair(intersection, true);
    }
    
    /// Map the position in to absolute coordinates
    SpaceVec map_to_absolute_space(const SpaceVec& pos) const {
        return pos % _domain_scale;
    }
    
    
    /// Map the position in to absolute coordinates
    SpaceVec map_to_relative_space(const SpaceVec& abs_pos) const {
        return abs_pos / _domain_scale;
    }  

    /// Getter for the size of the domain
    auto get_domain_size() const {
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
}; // struct CustomSpace

} // namespace Utopia::Models::PCPVertex::Space
#endif
