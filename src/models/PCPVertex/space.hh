#ifndef UTOPIA_MODELS_SAVANNAHETEROGENEOUS_SPACE_HH
#define UTOPIA_MODELS_SAVANNAHETEROGENEOUS_SPACE_HH

#include <cmath>
#include <limits>

#include <armadillo>

#include <utopia/data_io/cfg_utils.hh>
#include <utopia/core/types.hh>

namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace Space {

template<std::size_t num_dims>
struct CustomSpace : Utopia::Space<num_dims> {
    using Space = Utopia::Space<num_dims>;

    /// The type of a vector in space
    using SpaceVec = typename Space::SpaceVec;

    CustomSpace(const DataIO::Config& cfg)
    :
        Space(cfg)
    { }

    /// Constructor without any arguments, i.e. constructing a default space
    /** \details The default space is non-periodic and has default extent of 1.
      *         into each dimension.
      */
    CustomSpace()
    :
        Space()
    { }

    /// The displacement between two coordinates
    SpaceVec displacement(const SpaceVec& pos_0, const SpaceVec& pos_1) const {
        SpaceVec dx = pos_1 - pos_0;

        // The general case
        // Use a temporary vector to apply the transformation to
        auto mdx = dx;

        // Loop over displacement and extent and store the result of the
        // transformation in the temporary vector
        std::transform(dx.begin(), dx.end(),
                       this->extent.begin(), mdx.begin(),
            [](const double& d, const double& e){
                // Given the position in one dimension and the
                // corresponding extent, calculate the shorter displacement
                // (in units of the extend) to be back inside the space.
                return (d - std::round(d / e) * e);
            }
        );

        return mdx;
    }

    /// The distance of 2 coordinates in space
    /** \param p    the norm
     */
    double distance(const SpaceVec& pos_0, const SpaceVec& pos_1, int p=2) const
    {
        return norm(displacement(pos_0, pos_1), p);
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

}; // struct CustomSpace

template<std::size_t num_dims>
struct AbsoluteCustomSpace : CustomSpace<num_dims> {
    using Space = CustomSpace<num_dims>;

    /// The type of a vector in space
    using SpaceVec = typename Space::SpaceVec;

    /// The extend of the domain
    /** A remapping vector to a space extend in absolute coordinates.
     * 
     *  \note other than the extent this can be subject to Topological changes
     */
    SpaceVec domain_size;

    AbsoluteCustomSpace(const DataIO::Config& cfg)
    :
        Space(cfg),
        domain_size(this->extent)
    { }

    /// Constructor without any arguments, i.e. constructing a default space
    /** \details The default space is non-periodic and has default extent of 1.
      *         into each dimension.
      */
    AbsoluteCustomSpace()
    :
        Space(),
        domain_size(this->extent)
    { }

    /// Map the position in to the domain size
    SpaceVec convert_absolute(const SpaceVec& pos) const {
        return pos / this->extent % domain_size;
    }

    /// The absolute distance of 2 coordinates in space
    /** \param p    the norm
     */
    double distance(const SpaceVec& pos_0, const SpaceVec& pos_1,
                    int p=2) const
    {
        return norm(convert_absolute(this->displacement(pos_0, pos_1)), p);
    }
}; // struct AbsoluteCustomSpace

} // namespace Space
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia
#endif
