#ifndef UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH
#define UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH

#include "work_function.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {

/// @brief Boundary quadratic potential in shape of a rectangle
/// @tparam Model  
/** Use the WorkFunctionVertexTerm interface for a quadratic potential 
 *  penalising vertices on the outside of a rectangular box by their distance
 *  to the boundary.
 * 
 *  Parameters:
 *      - `elastic_constant` (double): The elastic constant
 *      - `width` (double): The width (along x-axis) of the rectangle
 *      - `height` (double): The height (along y-axis) of the rectangle
 *      - `origin` (SpaceVec<2>, optional): The center of the rectangle.
 *              If not provided the rectangle is centred to the barycenter of 
 *              the tissue's boundary.
 *      - `deform_plastic` (bool, default: false): If true, the width and height
 *              of a rectangle enclosing all vertices is taken and vertex 
 *              positions rescaled to the new width and height in a plastic 
 *              deformation. Note, this parameter is not inherited in 
 *              `update_parameters`.
 */
template <typename Model>
class BoundaryStripePotential : public WorkFunctionVertexTerm<Model>
{
public:
    using Base = WorkFunctionVertexTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

private:
    double _elastic_constant;

    double _width;

    double _height;

    SpaceVec _origin;

    SpaceVec get_lower_left () const {
        return _origin - 0.5 * SpaceVec({_width, _height});
    }

    SpaceVec get_extent () const {
        return SpaceVec({_width, _height});
    }

    /// The vector from a position to the closest boundary, if outside
    SpaceVec get_displacement (const SpaceVec& pos) const {
        const SpaceVec ll_corner = get_lower_left();
        const SpaceVec ur_corner = ll_corner + get_extent();

        SpaceVec displ({0., 0.});
        for (std::size_t axis = 0; axis <= 1; axis++) {
            if (pos[axis] < ll_corner[axis]) {
                displ[axis] = ll_corner[axis] - pos[axis];
            }
            else if (pos[axis] > ur_corner[axis]) {
                displ[axis] = ur_corner[axis] - pos[axis];
            }
        }

        return displ;
    }

    void deform_plastic () {
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::lowest();

        for (const auto &v : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(v);
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }

        double L = x_max - x_min;
        double H = y_max - y_min;

        SpaceVec extent = get_extent();
        if (_width > std::numeric_limits<double>::max() - 1) {
            extent[0] = L;
        }
        if (_height > std::numeric_limits<double>::max() - 1) {
            extent[1] = H;
        }
        SpaceVec scaling = extent / SpaceVec({L, H});

        for (const auto& v : this->_am.vertices()) {
            SpaceVec displ = this->_am.position_of(v) - _origin;
            this->_am.move_by(v, displ % (scaling - SpaceVec({1., 1.})));
        }
    }

public:
    BoundaryStripePotential (const DataIO::Config& cfg, const Model& model)
    :
        Base(cfg, model),
        _elastic_constant(get_as<double>("elastic_constant", cfg)),
        _width(get_as<double>("width", cfg)),
        _height(get_as<double>("height", cfg))
    {
        if (this->_am.get_space()->periodic) {
            throw std::runtime_error("Cannot setup stripe boundary potential "
                "in periodic space!");
        }

        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }
        else {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }

        if (_height < 1.e-8) {
            _height = std::numeric_limits<double>::max();
        }
        if (_width < 1.e-8) {
            _width = std::numeric_limits<double>::max();
        }

        if (get_as<bool>("deform_plastic", cfg, false)) {
            deform_plastic();
        }
    }

    SpaceVec compute_force(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec pos = this->_am.position_of(vertex);
        return _elastic_constant * get_displacement(pos);
    }

    double compute_energy(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec pos = this->_am.position_of(vertex);
        double distance = arma::norm(get_displacement(pos));
        return 0.5 * _elastic_constant * std::pow(distance, 2);
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _elastic_constant = get_as<double>("elastic_constant", cfg,
                                           _elastic_constant);
        _width = get_as<double>("width", cfg, _width);
        _height = get_as<double>("height", cfg, _height);
        
        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }

        if (get_as<bool>("deform_plastic", cfg, false)) {
            deform_plastic();
        }
    }
};


/// @brief Boundary quadratic potential in shape of a ring
/// @tparam Model  
/** Use the WorkFunctionVertexTerm interface for a quadratic potential 
 *  penalising vertices on the outside of a circular stripe of given thickness,
 *  length, and curvature.
 * 
 *  Parameters:
 *      - `elastic_constant` (double): The elastic constant
 *      - `curvature` (double): The curvature at the stripe's midline.
 *      - `width` (double): The width (along azimuthal-axis) of the stripe
 *      - `height` (double): The height (along radial-axis) of the stripe
 *      - `origin` (SpaceVec<2>, optional): The center of the stripe.
 *              If not provided the stripe is centred to the barycenter of 
 *              the tissue's boundary.
 *      - `deform_plastic` (bool, default: false): If true, the width and height
 *              of a stripe enclosing all vertices is taken and vertex 
 *              positions rescaled to the new width, height, and curvature
 *              in a plastic deformation. Note, this parameter is not inherited 
 *              in `update_parameters`.
 */
template <typename Model>
class BoundaryCurvedStripePotential : public WorkFunctionVertexTerm<Model>
{
public:
    using Base = WorkFunctionVertexTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

private:
    double _elastic_constant;

    double _curvature;

    double _width;

    double _height;

    SpaceVec _origin;

    /// The vector from a position to the closest boundary, if outside
    SpaceVec get_displacement (SpaceVec pos) const {
        if (fabs(_curvature) < 1.e-8) {
            const SpaceVec ll_corner = _origin - SpaceVec({_width, _height})/2.;
            const SpaceVec ur_corner = ll_corner + SpaceVec({_width, _height});

            SpaceVec displ({0., 0.});
            for (std::size_t axis = 0; axis <= 1; axis++) {
                if (pos[axis] < ll_corner[axis]) {
                    displ[axis] = ll_corner[axis] - pos[axis];
                }
                else if (pos[axis] > ur_corner[axis]) {
                    displ[axis] = ur_corner[axis] - pos[axis];
                }
            }

            return displ;
        }
        else {
            double R = 1. / _curvature;
            SpaceVec center = _origin - SpaceVec({0., R});

            SpaceVec coords = pos - center;

            double theta = atan2(coords[0], coords[1]);
            double r = arma::norm(coords);

            double boundary_r, boundary_theta;
            if (r - R < -_height / 2.) {
                boundary_r = R - _height / 2.;
            }
            else if (r - R > _height / 2.) {
                boundary_r = R + _height / 2.;
            }
            else {
                boundary_r = r;
            }

            double theta_max = _width / 2. / R;
            if (theta > theta_max) {
                boundary_theta = theta_max;
            }
            else if (theta < -theta_max) {
                boundary_theta = -theta_max;
            }
            else {
                boundary_theta = theta;
            }

            SpaceVec boundary = boundary_r * SpaceVec({sin(boundary_theta),
                                                       cos(boundary_theta)});

            return boundary - coords;
        }
    }

    void deform_plastic (double previous_curvature = 0.) {
        
        // remove the curvature and map to a rectangle
        if (fabs(previous_curvature) > 1.e-8) {
            double R = 1. / previous_curvature;

            SpaceVec center = _origin - SpaceVec({0., R});
            SpaceVec new_center = _origin - SpaceVec({0., 1.e8});

            for (const auto& vertex : this->_am.vertices()) {
                SpaceVec pos = this->_am.position_of(vertex) + center;

                double theta = atan2(pos[0], pos[1]);
                // the new radius with a curvature of 1.e-8
                double r = arma::norm(pos) - R + 1.e8;

                this->_am.move_to(
                    vertex,
                    r * SpaceVec({sin(theta), cos(theta)}) - new_center
                );
            }
        }


        // fit a rectangle around the vertices and rescale to width and height
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();
        double y_min = std::numeric_limits<double>::max();
        double y_max = std::numeric_limits<double>::lowest();

        for (const auto &v : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(v);
            x_min = std::min(x_min, pos[0]);
            x_max = std::max(x_max, pos[0]);
            y_min = std::min(y_min, pos[1]);
            y_max = std::max(y_max, pos[1]);
        }

        double L = x_max - x_min;
        double H = y_max - y_min;

        SpaceVec extent = SpaceVec({_width, _height});
        if (_width > std::numeric_limits<double>::max() - 1) {
            extent[0] = L;
        }
        if (_height > std::numeric_limits<double>::max() - 1) {
            extent[1] = H;
        }
        SpaceVec scaling = extent / SpaceVec({L, H});

        for (const auto& v : this->_am.vertices()) {
            SpaceVec displ = this->_am.position_of(v) - _origin;
            this->_am.move_by(v, displ % (scaling - SpaceVec({1., 1.})));
        }

        if (fabs(_curvature) < 1.e-8) {
            return;
        }

        // map to radial coordinates using curvature
        double R = 1. / _curvature;

        if (_width > 2 * M_PI * R) {
            throw std::runtime_error(fmt::format("Cannot map "
                "tissue to circle of radius {} (curvature {}) "
                "as it is longer than the circumference of the "
                "circle: {} > {}!",
                R, _curvature, _width, 2 * M_PI * R));
        }

        SpaceVec center = _origin - SpaceVec({0., R});

        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex) - _origin;
            // NOTE has no curvature here

            // the radius 
            double delta_r = pos[1];
            double r = delta_r + R;

            // the angle with arc-length delta_x
            double theta = pos[0] / R;

            SpaceVec new_pos = SpaceVec({r * sin(theta), r * cos(theta) - R});

            this->_am.move_to(
                vertex, new_pos + _origin);
        }
    }

public:
    BoundaryCurvedStripePotential (const DataIO::Config& cfg,
                                   const Model& model)
    :
        Base(cfg, model),
        _elastic_constant(get_as<double>("elastic_constant", cfg)),
        _curvature(get_as<double>("curvature", cfg)),
        _width(get_as<double>("width", cfg)),
        _height(get_as<double>("height", cfg))
    {
        if (this->_am.get_space()->periodic) {
            throw std::runtime_error("Cannot setup stripe boundary potential "
                "in periodic space!");
        }

        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }
        else {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }

        if (_height < 1.e-8) {
            _height = std::numeric_limits<double>::max();
        }
        if (_width < 1.e-8) {
            _width = std::numeric_limits<double>::max();
        }

        if (get_as<bool>("deform_plastic", cfg, false)) {
            model.get_logger()->warn(
                "Deforming the tissue plastically to a stripe with radius {}",
                1./_curvature
            );
            
            deform_plastic();
        }
    }

    SpaceVec compute_force(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec pos = this->_am.position_of(vertex);
        return _elastic_constant * get_displacement(pos);
    }

    double compute_energy(const std::shared_ptr<Vertex>& vertex) const final {
        SpaceVec pos = this->_am.position_of(vertex);
        double distance = arma::norm(get_displacement(pos));
        return 0.5 * _elastic_constant * std::pow(distance, 2);
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _elastic_constant = get_as<double>("elastic_constant", cfg,
                                           _elastic_constant);
        _width = get_as<double>("width", cfg, _width);
        _height = get_as<double>("height", cfg, _height);
        
        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }

        if (get_as<bool>("deform_plastic", cfg, false)) {
            deform_plastic(_curvature);
        }
        _curvature = get_as<double>("curvature", cfg, _curvature);
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
