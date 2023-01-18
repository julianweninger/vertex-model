#ifndef UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH
#define UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH

#include "work_function.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {



/// @brief Fixes the boundary vertices in space
template <typename Model>
class BoundaryFixed : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using Vertex = typename Base::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

public:
    BoundaryFixed (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model)
    { }

    void compute_and_set_forces () final {
        for (const auto& vertex : this->_am.vertices()) {
            if (this->_am.is_boundary(vertex)) {
                vertex->state.fix_in_space = true;
            }
        }
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        return 0.;
    }

    void update_parameters ([[maybe_unused]] const DataIO::Config& cfg) final
    { }
};


/// @brief Fixes a segment of the boundary vertices in space
/** 
 *  Parameters:
 *      - `lower_limit` (double): The lower limit of the fixed boundary. 
 *              0 being the point to the right of the barycenter of the tisuse 
 *              on the boundary. Can be in [0, 1], referencing a relative 
 *              distance from the point 0, measuring the length of the boundary.
 *      - `upper_limit` (double): The upper limit of the fixed boundary.
 *              Analogous to `lower_limit`.
 */
template <typename Model>
class BoundaryFixedPartial : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

    /// The lower limit of the fixed boundary.
    /** 0 being the point to the right of the barycenter of the tisuse on the 
     *  boundary. Can be in [0, 1], referencing a relative distance from the
     *  point 0, measuring the length of the boundary
    */
    double _lower_limit;
    
    /// The upper limit of the fixed boundary. See also lower limit.
    double _upper_limit;


public:
    BoundaryFixedPartial (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _lower_limit(get_as<double>("lower_limit", cfg)),
        _upper_limit(get_as<double>("upper_limit", cfg))
    { }

    void compute_and_set_forces () final {
        if (_lower_limit < 0 or _lower_limit > 1) {
            throw std::invalid_argument(fmt::format(
                "In BoundaryFixedPartial, 'lower_limit' has to be in [0, 1], "
                "but was {}", _lower_limit));
        }
        if (_upper_limit < 0 or _upper_limit > 1) {
            throw std::invalid_argument(fmt::format(
                "In BoundaryFixedPartial, 'lower_limit' has to be in [0, 1], "
                "but was {}", _upper_limit));
        }

        const auto boundary_edges = this->_am.get_boundary_edges();
        double boundary_length = 0.;
        std::vector<double> angles{};
        angles.reserve(boundary_edges.size());
        
        SpaceVec origin = this->_am.barycenter_of(boundary_edges);

        std::size_t i = 0;
        for (const auto& [edge, flip] : boundary_edges) {
            boundary_length += this->_am.length_of(edge);
            SpaceVec pos;
            if (not flip) {
                pos = this->_am.position_of(edge->custom_links().a);
            }
            else {
                pos = this->_am.position_of(edge->custom_links().b);
            }

            angles.push_back(fabs(atan2((pos - origin)[1], (pos - origin)[0])));
        }
        auto min = std::distance(
            angles.begin(),
            std::min_element(angles.begin(), angles.end())
        );

        double lower_limit = _lower_limit * boundary_length;
        double upper_limit = _upper_limit * boundary_length;
        
        double distance = 0.;
        for (i = 0; i < boundary_edges.size(); i++) {
            auto [edge, flip] = boundary_edges[(i + min)%boundary_edges.size()];
            
            if (distance >= lower_limit and distance <= upper_limit) {
                if (not flip) {
                    edge->custom_links().a->state.fix_in_space = true;
                }
                else {
                    edge->custom_links().b->state.fix_in_space = true;
                }
                distance += this->_am.length_of(edge);
            }
            else {
                break;
            }
        }
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        return 0.;
    }

    void update_parameters ([[maybe_unused]] const DataIO::Config& cfg) final
    {
        _lower_limit = get_as<double>("lower_limit", cfg, _lower_limit);
        _upper_limit = get_as<double>("upper_limit", cfg, _upper_limit);
    }
};

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
 * 
 *  TODO: Currently forces are applied in tangential / normal direction. 
 *        However, the derivative of these directors with displacement along
 *        the arc are not considered in a first order approximation.
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

    enum Quadrant {
        none,
        inner_right,
        inner_upper,
        inner_left,
        inner_lower,
        outer_right,
        outer_upper_right,
        outer_upper,
        outer_upper_left,
        outer_left,
        outer_lower_left,
        outer_lower,
        outer_lower_right,
    };

    bool is_inner_quadrant(Quadrant quadrant) const {
        if (quadrant == inner_right) {
            return true;
        }
        if (quadrant == inner_upper) {
            return true;
        }
        if (quadrant == inner_left) {
            return true;
        }
        if (quadrant == inner_lower) {
            return true;
        }
        return false;
    }


private:
    double _elastic_constant;

    double _curvature;

    double _width;

    double _height;

    SpaceVec _origin;

    bool _stretch;

    std::pair<SpaceVec, SpaceVec> corners () const {
        if (fabs(_curvature) < 1.e-8) {
            return std::make_pair(
                _origin - SpaceVec({_width, _height})/2.,
                _origin + SpaceVec({_width, _height})/2.
            );
        }
        throw std::runtime_error("");
    }
    std::pair<SpaceVec, SpaceVec> other_corners () const {
        if (fabs(_curvature) < 1.e-8) {
            const auto [ll_corner, ur_corner] = corners();
            return std::make_pair(
                ll_corner + SpaceVec({0, _height}),
                ur_corner - SpaceVec({0, _height})
            );
        }
        throw std::runtime_error("");
    }

    Quadrant get_quadrant (const std::shared_ptr<Vertex>& vertex) const {
        SpaceVec pos = this->_am.position_of(vertex);

        if (fabs(_curvature) < 1.e-8) {
            const auto [ll_corner, ur_corner] = corners();

            bool inside=true;
            Quadrant quadrant = none;

            if (pos[0] < ll_corner[0]) {
                inside=false;
                quadrant = outer_left;
            }
            else if (pos[0] > ur_corner[0]) {
                inside=false;
                quadrant = outer_right;
            }

            if (pos[1] < ll_corner[1]) {
                inside=false;
                if (quadrant==outer_left) {
                    quadrant = outer_lower_left;
                }
                else if (quadrant==outer_right) {
                    quadrant = outer_lower_right;
                }
                else {
                    quadrant = outer_lower;
                }
            }
            else if (pos[1] > ur_corner[1]) {
                inside=false;
                if (quadrant==outer_left) {
                    quadrant = outer_upper_left;
                }
                else if (quadrant==outer_right) {
                    quadrant = outer_upper_right;
                }
                else {
                    quadrant = outer_upper;
                }
            }

            if (inside) {
                double dx;
                if (pos[0] < _origin[0]) {
                    quadrant = inner_left;
                    dx = ll_corner[0] - pos[0];
                }
                else {
                    quadrant = inner_right;
                    dx = ur_corner[0] - pos[0];
                }
                
                if (pos[1] < _origin[1]) {
                    if (fabs(ll_corner[1] - pos[1]) < fabs(dx)) {
                        quadrant = inner_lower;
                    }
                    // else: keep q=left/right quadrant
                }
                else {
                    if (fabs(ur_corner[1] - pos[1]) < fabs(dx)) {
                        quadrant = inner_upper;
                    }
                    // else: keep q=left/right quadrant
                }
            }

            return quadrant;
        }
        
        // with curvature
        double R = 1. / _curvature;
        SpaceVec origin = _origin - SpaceVec({0., R});

        SpaceVec coords = pos - origin;

        double theta = atan2(coords[0], coords[1]);
        double r = arma::norm(coords);

        bool inside=true;
        Quadrant quadrant = none;

        if (r < R - _height / 2.) {
            inside = false;
            quadrant = outer_lower;
        }
        else if (r > R + _height / 2.) {
            inside = false;
            quadrant = outer_upper;
        }

        double theta_max = _width / 2. / R;
        if (theta > theta_max) {
            inside = false;
            if (quadrant==outer_lower) {
                quadrant = outer_lower_right;
            }
            else if (quadrant==outer_upper) {
                quadrant = outer_upper_right;
            }
            else {
                quadrant = outer_right;
            }
        }
        else if (theta < -theta_max) {
            inside = false;
            if (quadrant==outer_lower) {
                quadrant = outer_lower_left;
            }
            else if (quadrant==outer_upper) {
                quadrant = outer_upper_left;
            }
            else {
                quadrant = outer_left;
            }
        }
        
        if (inside) {
            double dr;
            if (r < R) {
                quadrant = inner_lower;
                dr = (R - _height / 2.) - r;
            }
            else {
                quadrant = inner_upper;
                dr = (R + _height / 2.) - r;
            }

            if (theta < 0) {
                double ds = r * (theta + theta_max);
                if (fabs(ds) < fabs(dr)) {
                    quadrant = inner_left;
                }
            }
            else {
                double ds = r * (theta - theta_max);
                if (fabs(ds) < fabs(dr)) {
                    quadrant = inner_right;
                }
            }
        }

        return quadrant;
    }

    SpaceVec get_normal (const std::shared_ptr<Vertex>& vertex) const {
        SpaceVec pos = this->_am.position_of(vertex);
        Quadrant quadrant = get_quadrant(vertex);
        
        if (fabs(_curvature) < 1.e-8) {
            const auto [ll_corner, ur_corner] = corners();
            const auto [ul_corner, lr_corner] = other_corners();

            if (quadrant == inner_right) {
                return SpaceVec({-1., 0.});
            }
            else if (quadrant == inner_upper) {
                return SpaceVec({0., -1.});
            }
            else if (quadrant == inner_left) {
                return SpaceVec({1., 0.});
            }
            else if (quadrant == inner_lower) {
                return SpaceVec({0., 1.});
            }
            else if (quadrant == outer_right) {
                return SpaceVec({1., 0.});
            }
            else if (quadrant == outer_upper_right) {
                SpaceVec displ = ur_corner - pos;
                return displ / arma::norm(displ);
            }
            else if (quadrant == outer_upper) {
                return SpaceVec({0, 1.});
            }
            else if (quadrant == outer_upper_left) {
                SpaceVec displ = ul_corner - pos;
                return displ / arma::norm(displ);
            }
            else if (quadrant == outer_left) {
                return SpaceVec({-1., 0.});
            }
            else if (quadrant == outer_lower_left) {
                SpaceVec displ = ll_corner - pos;
                return displ / arma::norm(displ);
            }
            else if (quadrant == outer_lower) {
                return SpaceVec({0., 1.});
            }
            else if (quadrant == outer_lower_right) {
                SpaceVec displ = lr_corner - pos;
                return displ / arma::norm(displ);
            }
            else {
                throw std::runtime_error("Unknown quadrant!");
            }
        }
        
        // with curvature
        double R = 1. / _curvature;
        SpaceVec origin = _origin - SpaceVec({0., R});

        SpaceVec rpos = pos - origin;
        double theta = std::atan2(rpos[0], rpos[1]);

        // vectors of the arc
        SpaceVec normal( {sin(- theta), -cos(- theta)});
        SpaceVec tangent({cos(- theta),  sin(- theta)});

        // opening angle of the arc
        double theta_max = _width / 2. / R;

        if (quadrant == inner_right) {
            return tangent;
        }
        else if (quadrant == inner_upper) {
            return -1. * normal;
        }
        else if (quadrant == inner_left) {
            return -1. * tangent;
        }
        else if (quadrant == inner_lower) {
            return normal;
        }
        else if (quadrant == outer_right) {
            return -1. * tangent;
        }
        else if (quadrant == outer_upper_right) {
            double uR = R + _height / 2.;
            SpaceVec ur_corner = uR * SpaceVec({sin(theta_max),
                                                cos(theta_max)});
            SpaceVec displ = ur_corner - rpos;
            return displ / arma::norm(displ);
        }
        else if (quadrant == outer_upper) {
            return normal;
        }
        else if (quadrant == outer_upper_left) {
            double uR = R + _height / 2.;
            SpaceVec ul_corner = uR * SpaceVec({sin(-theta_max),
                                                cos(-theta_max)});
            SpaceVec displ = ul_corner - rpos;
            return displ / arma::norm(displ);
        }
        else if (quadrant == outer_left) {
            return tangent;
        }
        else if (quadrant == outer_lower_left) {
            double lR = R - _height / 2.;
            SpaceVec ll_corner = lR * SpaceVec({sin(-theta_max),
                                                cos(-theta_max)});
            SpaceVec displ = ll_corner - rpos;
            return displ / arma::norm(displ);
        }
        else if (quadrant == outer_lower) {
            return -1. * normal;
        }
        else if (quadrant == outer_lower_right) {
            double lR = R - _height / 2.;
            SpaceVec lr_corner = lR * SpaceVec({sin(theta_max),
                                                cos(theta_max)});
            SpaceVec displ = lr_corner - rpos;
            return displ / arma::norm(displ);
        }
        else {
            throw std::runtime_error("Unknown quadrant!");
        }
    }


    /// The vector from a position to the closest boundary, if outside
    SpaceVec get_displacement (const std::shared_ptr<Vertex>& vertex) const {
        SpaceVec pos = this->_am.position_of(vertex);
        Quadrant quadrant = get_quadrant(vertex);
        if (fabs(_curvature) < 1.e-8) {
            const auto [ll_corner, ur_corner] = corners();
            const auto [ul_corner, lr_corner] = other_corners();

            const SpaceVec right = 0.5 * (lr_corner + ur_corner);
            const SpaceVec upper = 0.5 * (ul_corner + ur_corner);
            const SpaceVec left = 0.5 * (ll_corner + ul_corner);
            const SpaceVec lower = 0.5 * (ll_corner + lr_corner);

            SpaceVec normal = get_normal(vertex);

            if (quadrant == inner_right or quadrant == outer_right) {
                return arma::dot((pos - right), normal) * normal;
            }
            else if (quadrant == inner_upper or quadrant == outer_upper) {
                return arma::dot((pos - upper), normal) * normal;
            }
            else if (quadrant == inner_left or quadrant == outer_left) {
                return arma::dot((pos - left), normal) * normal;
            }
            else if (quadrant == inner_lower or quadrant == outer_lower) {
                return arma::dot((pos - lower), normal) * normal;
            }
            else if (quadrant == outer_upper_right) {
                return arma::dot((pos - ur_corner), normal) * normal;
            }
            else if (quadrant == outer_upper_left) {
                return arma::dot((pos - ul_corner), normal) * normal;
            }
            else if (quadrant == outer_lower_left) {
                return arma::dot((pos - ll_corner), normal) * normal;
            }
            else if (quadrant == outer_lower_right) {
                return arma::dot((pos - lr_corner), normal) * normal;
            }
            else {
                throw std::runtime_error("Unknown quadrant!");
            }
            return SpaceVec({0., 0.});
        }
        
        // with curvature
        double R = 1. / _curvature;
        SpaceVec origin = _origin - SpaceVec({0., R});

        SpaceVec rpos = pos - origin;
        double r = arma::norm(rpos);
        double theta = std::atan2(rpos[0], rpos[1]);
        double theta_max = _width / 2. / R;

        SpaceVec normal = get_normal(vertex);

        if (quadrant == inner_right) {
            return - r * fabs(theta_max - theta) * normal;
        }
        else if (quadrant == inner_upper) {
            return - fabs((R + _height / 2.) - r) * normal;
        }
        else if (quadrant == inner_left) {
            return - r * fabs(-theta_max - theta) * normal;
        }
        else if (quadrant == inner_lower) {
            return - fabs((R - _height / 2. - r)) * normal;
        }
        else if (quadrant == outer_right) {
            return - r * fabs(theta - theta_max) * normal;
        }
        else if (quadrant == outer_upper_right) {
            double uR = R + _height / 2.;
            SpaceVec ur_corner = uR * SpaceVec({sin(theta_max),
                                                cos(theta_max)});
            return rpos - ur_corner;
        }
        else if (quadrant == outer_upper) {
            return - fabs(r - (R + _height / 2.)) * normal;
        }
        else if (quadrant == outer_upper_left) {
            double uR = R + _height / 2.;
            SpaceVec ul_corner = uR * SpaceVec({sin(-theta_max),
                                                cos(-theta_max)});
            return rpos - ul_corner;
        }
        else if (quadrant == outer_left) {
            return - r * fabs(theta + theta_max) * normal;
        }
        else if (quadrant == outer_lower_left) {
            double lR = R - _height / 2.;
            SpaceVec ll_corner = lR * SpaceVec({sin(-theta_max),
                                                cos(-theta_max)});
            return rpos - ll_corner;
        }
        else if (quadrant == outer_lower) {
            return - fabs(r - (R - _height / 2.)) * normal;
        }
        else if (quadrant == outer_lower_right) {
            double lR = R - _height / 2.;
            SpaceVec lr_corner = lR * SpaceVec({sin(theta_max),
                                                cos(theta_max)});
            return rpos - lr_corner;
        }
        else {
            throw std::runtime_error("Unknown quadrant!");
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
        else if (not model.get_space()->periodic) {
            _origin = this->_am.barycenter_of(this->_am.get_boundary_edges());
        }
        else {
            _origin = model.get_space()->get_domain_size() / 2.;
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
        if (not this->_am.is_boundary(vertex)) {
            return SpaceVec({0., 0.});
        }

        const auto quadrant = get_quadrant(vertex);
        if (not _stretch and is_inner_quadrant(quadrant)) {
            return SpaceVec({0., 0.});
        }
        
        if (fabs(_curvature) < 1.e-10) {
            return -1. * _elastic_constant * get_displacement(vertex);
        }
        
        // with curvature
        if (this->_am.get_space()->periodic) {
            throw std::runtime_error("No curvature in periodic BC implemented");
        }

        SpaceVec displ = get_displacement(vertex);
        return -1. * _elastic_constant * get_displacement(vertex);

        // TODO --->

        // derivative of axis to displacement of vertex

        // SpaceVec origin = _origin - SpaceVec({0., 1. / _curvature});
        // SpaceVec rpos = this->_am.position_of(vertex) - origin;
        
        // SpaceVec normal = displ / arma::norm(displ);

        // // SpaceVec tmp;
        // // if (quadrant == inner_right) {
        // //     double x = rpos[0];
        // //     double y = rpos[1];
        // //     double x2 = x*x;
        // //     double y2 = y*y;
        // //     double _tmp = pow(x2 / y2 + 1, 1.5);
        // //     auto dtx_dx = -x  / (y2 * _tmp)
        // //     auto dty_dx =  x2 / (y2 * y * _tmp);

        // //     _tmp = (sqrt(x/y2 + 1) (rpos[0]^2 + y2));
        // //     auto dtx_dy = -rpos[1] / _tmp;
        // //     auto dty_dy =  rpos[0] / _tmp;
        // // }
        
        // // // use tangent ..

        // double sum_x = rpos[0];
        // double sum_y = rpos[1];
        // double arg = sum_x / sum_y;

        // double tmp_0 = sum_y * (pow(arg, 2) + 1);
        // double tmp_1 = pow(sum_y, 2) * (pow(arg, 2) + 1);

        // double dpx_dx =  normal[1] / tmp_0;
        // double dpy_dx = -normal[0] / tmp_0;

        // double dpx_dy = -(sum_x * normal[1]) / tmp_1;
        // double dpy_dy =  (sum_x * normal[0]) / tmp_1;

        // SpaceVec tmp({
        //     displ[0] * dpx_dx + displ[1] * dpy_dx,
        //     displ[0] * dpx_dy + displ[1] * dpy_dy,
        // });
        // double _dE_dx = _elastic_constant;

        // // d/dx(1/sqrt(x^2/y^2 + 1)) = -x/(y^2 (x^2/y^2 + 1)^(3/2))
        // // d/dy(1/sqrt(x^2/y^2 + 1)) = x^2/(y^3 (x^2/y^2 + 1)^(3/2))   

        // // d/dx(-x/(y sqrt(x^2/y^2 + 1))) = -y/(sqrt(x^2/y^2 + 1) (x^2 + y^2))
        // // d/dy(-x/(y sqrt(x^2/y^2 + 1))) = x/(sqrt(x^2/y^2 + 1) (x^2 + y^2))

        // return - _dE_dx * (displ - tmp);
    }

    double compute_energy(const std::shared_ptr<Vertex>& vertex) const final {
        if (not this->_am.is_boundary(vertex)) {
            return 0.;
        }
        if (not _stretch and is_inner_quadrant(get_quadrant(vertex))) {
            return 0.;
        }

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
