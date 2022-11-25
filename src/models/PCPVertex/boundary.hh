#ifndef UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH
#define UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH

#include "work_function.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {


/// @brief Boundary quadratic potential in shape of a rectangle or ring
/// @tparam Model  
/** Use the WorkFunctionVertexTerm interface for a quadratic potential 
 *  penalising vertices on the outside of a rectangle or circular stripe of
 *  given thickness, length, and curvature.
 * 
 *  Parameters:
 *      - `elastic_constant` (double): The elastic constant
 *      - `curvature` (double, default: 0): The curvature at the stripe's
 *              midline.
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
 *      - `stretch` (bool): Whether to pull boundary vertices that are inside
 *              the domain towards the closest boundary.
 */
template <typename Model>
class BoundaryStripePotential : public WorkFunctionTerm<Model>
{
public:
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Base::Edge;

    using Cell = typename Base::Cell;

private:
    double _elastic_constant;

    double _curvature;

    double _width;

    double _height;

    SpaceVec _origin;

    bool _stretch;

    /// The vector from a position to the closest boundary, if outside
    SpaceVec get_displacement (const std::shared_ptr<Vertex>& vertex) const {
        SpaceVec pos = this->_am.position_of(vertex);
        if (fabs(_curvature) < 1.e-8) {
            const SpaceVec ll_corner = _origin - SpaceVec({_width, _height})/2.;
            const SpaceVec ur_corner = ll_corner + SpaceVec({_width, _height});

            SpaceVec displ({0., 0.});
            bool apply_stretch = this->_am.is_boundary(vertex);
            for (std::size_t axis = 0; axis <= 1; axis++) {
                if (pos[axis] < ll_corner[axis]) {
                    displ[axis] = ll_corner[axis] - pos[axis];
                    apply_stretch = false;
                }
                else if (pos[axis] > ur_corner[axis]) {
                    displ[axis] = ur_corner[axis] - pos[axis];
                    apply_stretch = false;
                }
            }

            // use only axis along which distance is shorter
            if (_stretch and apply_stretch) {
                for (std::size_t axis = 0; axis <= 1; axis++) {
                    if (pos[axis] < _origin[axis]) {
                        displ[axis] = ll_corner[axis] - pos[axis];
                    }
                    else {
                        displ[axis] = ur_corner[axis] - pos[axis];
                    }
                }

                if (fabs(displ[1]) > fabs(displ[0])) {
                    displ[1] = 0.;
                }
                else {
                    displ[0] = 0.;
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
            bool apply_stretch = this->_am.is_boundary(vertex);
            if (r < R - _height / 2.) {
                boundary_r = R - _height / 2.;
                apply_stretch = false;
            }
            else if (r > R + _height / 2.) {
                boundary_r = R + _height / 2.;
                apply_stretch = false;
            }
            else {
                boundary_r = r;
            }

            double theta_max = _width / 2. / R;
            if (theta > theta_max) {
                boundary_theta = theta_max;
                apply_stretch = false;
            }
            else if (theta < -theta_max) {
                boundary_theta = -theta_max;
                apply_stretch = false;
            }
            else {
                boundary_theta = theta;
            }
            
            if (_stretch and apply_stretch) {
                if (r < R) {
                    boundary_r = R - _height / 2.;
                }
                else {
                    boundary_r = R + _height / 2.;
                }

                if (theta < 0) {
                    boundary_theta = -theta_max;
                }
                else {
                    boundary_theta = theta_max;
                }

                // get closer boundary
                double dr = fabs(r - boundary_r);
                double ds = r * fabs(theta - boundary_theta);

                if (dr > ds) {
                    boundary_r = r;
                }
                else {
                    boundary_theta = theta;
                }
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
                SpaceVec pos = this->_am.position_of(vertex) - _origin;

                double delta_r = arma::norm(pos + SpaceVec({0., R})) - R;
                double theta = atan2(pos[0], pos[1] + R);

                this->_am.move_to(
                    vertex,
                    SpaceVec({theta * R, delta_r}) + _origin
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

        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex) - _origin;
            // NOTE has no curvature here

            // the radius 
            double delta_r = pos[1];
            double r = delta_r + R;

            // the angle with arc-length delta_x
            double theta = pos[0] / R;

            this->_am.move_to(
                vertex, 
                SpaceVec({r * sin(theta), r * cos(theta) - R}) + _origin
            );
        }
    }

public:
    BoundaryStripePotential (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elastic_constant(get_as<double>("elastic_constant", cfg)),
        _curvature(get_as<double>("curvature", cfg, 0.)),
        _width(get_as<double>("width", cfg)),
        _height(get_as<double>("height", cfg)),
        _origin(SpaceVec({0., 0.})),
        _stretch(get_as<bool>("stretch", cfg))
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
            model.get_logger()->info(
                "Deforming the tissue plastically to a stripe "
                "with curvature {}",
                _curvature
            );
            
            deform_plastic();
        }
    }

    void compute_and_set_forces () final {
        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.add_force(compute_force(vertex));
        }
    }

    SpaceVec compute_force(const std::shared_ptr<Vertex>& vertex) const final {
        return _elastic_constant * get_displacement(vertex);
    }

    double compute_energy(const std::shared_ptr<Vertex>& vertex) const final {
        double distance = arma::norm(get_displacement(vertex));
        return 0.5 * _elastic_constant * std::pow(distance, 2);
    }

    double compute_energy (
        const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double energy = 0.;
        for (const auto& vertex : vertices) {
            energy += compute_energy(vertex);
        }
        return energy;
    }

    void update_parameters (const DataIO::Config& cfg) final {
        _elastic_constant = get_as<double>("elastic_constant", cfg,
                                           _elastic_constant);
        _width = get_as<double>("width", cfg, _width);
        _height = get_as<double>("height", cfg, _height);
        
        if (cfg["origin"]) {
            _origin = get_as_SpaceVec<2>("origin", cfg);
        }

        double tmp_curvature = _curvature;
        _curvature = get_as<double>("curvature", cfg, _curvature);
        if (get_as<bool>("deform_plastic", cfg, false)) {
            deform_plastic(tmp_curvature);
        }
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
