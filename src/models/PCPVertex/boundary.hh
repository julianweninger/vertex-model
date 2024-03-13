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
 *      - `from_left` (double, default: 0): The domain size from the left most
 *          vertex in which all boundary vertices are fixed.
 *      - `from_right` (double, default: 0): Same as `from_left`, measured from
 *          the right most vertex.
 *      - `from_top` (double, default: 0): Same as `from_left`, measured from
 *          the top most vertex.
 *      - `from_bottom` (double, default: 0): Same as `from_left`, measured from
 *          the bottom most vertex.
 *      - `fix` (bool): If false, the selection of vertices is updated every 
 *          iteration.
 */
template <typename Model>
class BoundaryFixedPartial : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::AgentManager::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

private:
    /// The distance from the right most vertex which are fixed
    double _right;
    
    double _top;

    double _left;

    double _bottom;

    bool _fix;

    const std::string _param_fixed;

    void update_fixation () {
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::min();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::min();
        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex);

            min_x = std::min(min_x, pos[0]);
            max_x = std::max(max_x, pos[0]);
            min_y = std::min(min_y, pos[1]);
            max_y = std::max(max_y, pos[1]);
        }

        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex);
            if (   pos[0] + 1.e-8 < min_x + _left
                or pos[0] - 1.e-8 > max_x - _right
                or pos[1] + 1.e-8 < min_y + _bottom
                or pos[1] - 1.e-8 > max_y - _top
            ) {
                if (this->_am.is_boundary(vertex)) {
                    vertex->state.update_parameter(_param_fixed, true);
                }
                else {
                    vertex->state.update_parameter(_param_fixed, false);
                }
            }
            else {
                vertex->state.update_parameter(_param_fixed, false);
            }
        }
    }

public:
    BoundaryFixedPartial (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),

        _right(get_as<double>("from_right", cfg, 0.)),
        _top(get_as<double>("from_top", cfg, 0.)),
        _left(get_as<double>("from_left", cfg, 0.)),
        _bottom(get_as<double>("from_bottom", cfg, 0.)),
        _fix(get_as<bool>("fix", cfg)),
        _param_fixed(this->_name + "_fixed")
    {
        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.register_parameter(_param_fixed, false);
        }

        update_fixation();
    }

    ~BoundaryFixedPartial() {
        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.unregister_parameter(_param_fixed);
        }
    }

    void compute_and_set_forces () final {
        if (not _fix) {
            update_fixation();
        }

        for (const auto& vertex : this->_am.vertices()) {
            if (not vertex->state.has_parameter(_param_fixed)) {
                vertex->state.register_parameter(_param_fixed, false);
            }
            else if (vertex->state.get_parameter(_param_fixed)) {
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
    {
        bool change = false;

        if (cfg["from_right"]) {
            _right = get_as<double>("from_right", cfg);
            change = true;
        }
        if (cfg["from_top"]) {
            _top = get_as<double>("from_top", cfg);
            change = true;
        }
        if (cfg["from_left"]) {
            _left = get_as<double>("from_left", cfg);
            change = true;
        }
        if (cfg["from_bottom"]) {
            _bottom = get_as<double>("from_bottom", cfg);
            change = true;
        }
        if (cfg["fix"]) {
            _fix = get_as<bool>("fix", cfg);
            change = true;
        }

        if (change) {
            update_fixation();
        }
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
 *      - `skip_quadrants` (list[int]): Which quadrants to exclude from bc.
 *              inner quadrants 1-4 from right anti-clockwise. outer quadrants
 *                 5-12 from right, including corners.
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
        inner_right=1,
        inner_upper=2,
        inner_left=3,
        inner_lower=4,
        outer_right=5,
        outer_upper_right=6,
        outer_upper=7,
        outer_upper_left=8,
        outer_left=9,
        outer_lower_left=10,
        outer_lower=11,
        outer_lower_right=12,
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


protected:
    double _elastic_constant;

    double _curvature;

    double _width;

    double _height;

    SpaceVec _origin;

    bool _stretch;

    std::set<Quadrant> _skip_quadrants;

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
        if (not this->_am.is_boundary(vertex)) {
            return Quadrant::none;
        }

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
        _stretch(get_as<bool>("stretch", cfg)),
        _skip_quadrants{}
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

        auto skip_quadrants = get_as<std::vector<std::size_t>>(
            "skip_quadrants", cfg, {}
        );
        for (auto q : skip_quadrants) {
            _skip_quadrants.insert(Quadrant(q));
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

        if (_skip_quadrants.find(quadrant) != _skip_quadrants.end()) {
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
        const auto quadrant = get_quadrant(vertex);
        if (not _stretch and is_inner_quadrant(quadrant)) {
            return 0.;
        }

        if (_skip_quadrants.find(quadrant) != _skip_quadrants.end()) {
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

    void update_parameters (const DataIO::Config& cfg) override {
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

/// @brief  BoundaryStripePotential with pure shear variation of domain size
/// @tparam Model 
/** Width and height of the BoundaryStripePotential are changed with a pure 
 *  shear velocity, conserving total domain area.
 * 
 *  Arguments:
 *      - `pure_shear_velocity` (double): The pure shear velocity vxx.
 *      - `width` (double): Width of the domain. Cannot be updated.
 *      - `height` (double): Height of the domain. Cannot be updated.
 *      - `curvature` (double, fixed value: 0): Restricted to zero curvature
 *              compared to BoundaryStripePotential.
*/
template <typename Model>
class BoundaryStripePotentialPureShear : public BoundaryStripePotential<Model>
{
    using Base = BoundaryStripePotential<Model>;

    using SpaceVec = typename Base::SpaceVec;

protected:
    /// @brief The pure shear velocity v_xx
    /** In units of model time */
    double _pure_shear_velocity;

    /// The cumulative pure shear
    double _cum_pure_shear;

    /// @brief  The domain width at t = 0
    const double _original_width;

    /// @brief  The domain width at t = 0.
    const double _original_height;

public:
    BoundaryStripePotentialPureShear (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _pure_shear_velocity(get_as<double>("pure_shear_velocity", cfg)),
        _cum_pure_shear(0.),
        _original_width(this->_width),
        _original_height(this->_height)
    {
        if (fabs(this->_curvature) > 1.e-6) {
            throw std::runtime_error("Cannot operate "
                "BoundaryStripePotentialPureShear WF with non zero curvature!"
            );
        }
    }

    void update (double dt) override {
        _cum_pure_shear += _pure_shear_velocity * dt;
        this->_width = _original_width * exp(_cum_pure_shear);
        this->_height = _original_height * exp(-_cum_pure_shear);
    }

    void update_parameters (const DataIO::Config& cfg) override {
        if (cfg["width"] or cfg["height"]) {
            throw std::runtime_error("Cannot update width or height in "
                "BoundaryStripePotentialPureShear WF! Please use the "
                "pure_shear_velocity to do so."
            );
        }

        Base::update_parameters(cfg);
        if (fabs(this->_curvature) > 1.e-6) {
            throw std::runtime_error("Cannot operate "
                "BoundaryStripePotentialPureShear WF with non zero curvature!"
            );
        }

        _pure_shear_velocity = get_as<double>(
            "pure_shear_velocity", cfg, _pure_shear_velocity
        );
    }
};


/// @brief  Classifies the boundary edges to 4 quadrants
/// @tparam Model 
/** Classifier the boundary edges to types BoundaryClassifier::Quadrant.
 * 
 *  Initially, the boundaries are identified by the center point of each 
 *  boundary edge and the distance to the extreme extents of the tissue in
 *  x, y calculated. In order, top, bottom, right, and left, the edges
 *  are classified if they are within a certain distance of each boundary.
 * 
 *  In every update, boundary edges with type Quadrant::none are updated to 
 *  the quadrant of the edge that precedes it in anti-clockwise orientation.
 *  Note, that initially unassigned boundary edges will eventually inherit a
 *  classification by this update.
 * 
 *  Note, for best results, initialise this work function in a state of a
 *  rectangular tissue.
 * 
 *  Arguments:
 *      - `from_top`: Distance from maximum coordinate in y
 *      - `from_bottom`: Distance from minimum coordinate in y
 *      - `from_right`: Distance from maximum coordinate in x
 *      - `from_left`: Distance from minimum coordinate in x
*/
template <typename Model>
class BoundaryClassifier : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

public:
    enum Quadrant {
        none=0,
        top=1,
        right=2,
        bottom=3,
        left=4
    };

protected:
    /// The name of the parameter where the classification is stored
    const std::string _classifier;

public:
    BoundaryClassifier (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _classifier(name + "_boundary_class")
    {

        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::min();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::min();
        for (const auto& edge : this->_am.edges()) {
            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            SpaceVec pos = 0.5 * (  this->_am.position_of(a) 
                                  + this->_am.position_of(b));

            min_x = std::min(min_x, pos[0]);
            max_x = std::max(max_x, pos[0]);
            min_y = std::min(min_y, pos[1]);
            max_y = std::max(max_y, pos[1]);
        }

        double top = get_as<double>("from_top", cfg);
        double bottom = get_as<double>("from_bottom", cfg);
        double right = get_as<double>("from_right", cfg);
        double left = get_as<double>("from_left", cfg);

        for (const auto& edge : this->_am.edges()) {
            const auto& a = edge->custom_links().a;
            const auto& b = edge->custom_links().b;

            SpaceVec pos = 0.5 * (  this->_am.position_of(a) 
                                  + this->_am.position_of(b));

            if (not this->_am.is_1_cell_boundary_edge(edge)) {
                edge->state.register_parameter(_classifier, Quadrant::none);
            }
            else if (pos[1] - 1.e-8 > max_y - top) {
                edge->state.register_parameter(_classifier, Quadrant::top);
            }
            else if (pos[1] + 1.e-8 < min_y + bottom) {
                edge->state.register_parameter(_classifier, Quadrant::bottom);
            }
            else if (pos[0] - 1.e-8 > max_x - right) {
                edge->state.register_parameter(_classifier, Quadrant::right);
            }
            else if (pos[0] + 1.e-8 < min_x + left) {
                edge->state.register_parameter(_classifier, Quadrant::left);
            }
            else {
                edge->state.register_parameter(_classifier, Quadrant::none);
            }
        }

        for (std::size_t i = 0; i < 5; i++) {
            this->update(0.);
        }

        // for (const auto& edge : this->_am.edges()) {
        //     auto [ca, cb] = this->_am.adjoints_of(edge);

        //     if (ca and cb) {
        //         // bulk
        //         edge->state.register_parameter(_classifier, 0);
        //         continue;
        //     }

        //     if (not ca) { std::swap(ca, cb); }

        //     SpaceVec displ({0., 1.});
        //     for (const auto& [e, flip] : ca->custom_links().edges) {
        //         if (e == edge) {
        //             displ = this->_am.displacement(edge, flip);
        //             break;
        //         }
        //     }

        //     double angle = atan2(displ[1], displ[0]) / M_PI * 180;
        //     double crit_angle = get_as<double>("separation_angle", cfg);
        //     if (fabs(angle) > 90 + crit_angle) {
        //         // left oriented junction (top)
        //         edge->state.register_parameter(_classifier, Quadrant::top);
        //     }
        //     else if (angle > 90 - crit_angle) {
        //         // top oriented junction (right)
        //         edge->state.register_parameter(_classifier, Quadrant::right);
        //     }
        //     else if (angle > -90 + crit_angle) {
        //         // right oriented junction (bottom)
        //         edge->state.register_parameter(_classifier, Quadrant::bottom);
        //     }
        //     else {
        //         // down oriented junction (left)
        //         edge->state.register_parameter(_classifier, Quadrant::left);
        //     }
        // }
    }

    ~BoundaryClassifier() {
        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.unregister_parameter(_classifier);
        }
        for (const auto& edge : this->_am.edges()) {
            edge->state.unregister_parameter(_classifier);
        }
    }

    void compute_and_set_forces () final { }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        return 0.;
    }

    void update ([[maybe_unused]] double dt) final {
        for (const auto& edge : this->_am.edges()) {
            if (not edge->state.has_parameter(_classifier)) {
                edge->state.register_parameter(_classifier, 0);
            }
        }

        const auto edges = this->_am.get_boundary_edges();
        
        // every new junction inherits its quadrant from the previous junction
        // where junctions are order anti-clockwise along the boundary
        const auto& [edge, flip] = edges.back();
        Quadrant q = (Quadrant)(edge->state.get_parameter(_classifier) + .5);        
        for (const auto& [edge, flip] : edges) {
            Quadrant _q = (Quadrant)(edge->state.get_parameter(_classifier)+.5);

            if (_q == 0) {
                _q = q;
            }

            edge->state.update_parameter(_classifier, _q);
            
            q = _q;
        }
    }

    void update_parameters ([[maybe_unused]] const DataIO::Config& cfg) final
    {
        return;
    }

    std::vector<std::string> write_task_edge_properties_names () const final {
        return std::vector<std::string>({
            "boundary_class"
        });
    }

    std::vector<std::vector<double>> write_edge_properties () const final {
        std::vector<double> contractilities({});
        for (const auto& edge : this->_am.edges()) {
            contractilities.push_back(edge->state.get_parameter(_classifier));
        }
        return std::vector<std::vector<double>>({contractilities});
    }

    Quadrant get_class (const std::shared_ptr<Edge>& edge) const {
        return (Quadrant)(edge->state.get_parameter(_classifier)+.5);
    }
};


/// @brief Bending elasticity on boundary
/** \f$ W = 1/2 B \kappa^2 \f$ 
 * 
 *  Parameters:
 *      - `from_left` (double, default: 0): The domain size from the left most
 *          vertex in which all boundary vertices are fixed.
 *      - `from_right` (double, default: 0): Same as `from_left`, measured from
 *          the right most vertex.
 *      - `from_top` (double, default: 0): Same as `from_left`, measured from
 *          the top most vertex.
 *      - `from_bottom` (double, default: 0): Same as `from_left`, measured from
 *          the bottom most vertex.
*/
template <typename Model>
class BoundaryBendElastic : public WorkFunctionTerm<Model>
{
    using Base = WorkFunctionTerm<Model>;

    using SpaceVec = typename Base::SpaceVec;

    using Vertex = typename Base::Vertex;

    using Edge = typename Model::Edge;

    using Cell = typename Model::Cell;

private:
    /// The elastic modulus
    double _elasticity;

    /// The distance from the right most vertex which are considered
    double _right;
    
    double _top;

    double _left;

    double _bottom;

    const std::string _param_bend;

    enum BC {
        periodic,
        open,
        fixed
    } _boundary_condition;

    /// The angles at which to fix the junction at the boundary
    std::pair<double, double> _fixed_bc;

    /// the anchoring points in fixed BC that comply with specified angles
    std::shared_ptr<std::pair<SpaceVec, SpaceVec>> _fixed_bc__origin;

    void update_contribution () {
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::min();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::min();
        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex);

            min_x = std::min(min_x, pos[0]);
            max_x = std::max(max_x, pos[0]);
            min_y = std::min(min_y, pos[1]);
            max_y = std::max(max_y, pos[1]);
        }

        for (const auto& vertex : this->_am.vertices()) {
            SpaceVec pos = this->_am.position_of(vertex);
            if (   pos[0] + 1.e-8 < min_x + _left
                or pos[0] - 1.e-8 > max_x - _right
                or pos[1] + 1.e-8 < min_y + _bottom
                or pos[1] - 1.e-8 > max_y - _top
            ) {
                if (this->_am.is_boundary(vertex)) {
                    vertex->state.update_parameter(_param_bend, true);
                }
                else {
                    vertex->state.update_parameter(_param_bend, false);
                }
            }
            else {
                vertex->state.update_parameter(_param_bend, false);
            }
        }
    }

    auto get_boundary () const {
        auto boundary = this->_am.get_boundary_edges();

        if (_boundary_condition == BC::periodic) {
            return boundary;
        }

        std::size_t index = 0;
        for (std::size_t i = 0; i < boundary.size(); i++) {
            const auto& [edge, flip] = boundary[i];
            std::shared_ptr<Vertex> vertex;
            if (not flip) {
                vertex = edge->custom_links().a;
            }
            else {
                vertex = edge->custom_links().b;
            }

            if (not vertex->state.has_parameter(_param_bend)) {
                continue;
            }
            if (not vertex->state.get_parameter(_param_bend)) {
                index = i;
                break;
            }
        }

        std::rotate(
            boundary.begin(),
            boundary.begin() + index,
            boundary.end()
        );
        // the first vertex is a non contributing vertex

        index = 0;
        for (std::size_t i = 0; i < boundary.size(); i++) {
            const auto& [edge, flip] = boundary[i];
            std::shared_ptr<Vertex> vertex;
            if (not flip) {
                vertex = edge->custom_links().a;
            }
            else {
                vertex = edge->custom_links().b;
            }

            if (not vertex->state.has_parameter(_param_bend)) {
                continue;
            }
            if (vertex->state.get_parameter(_param_bend)) {
                index = i;
                break;
            }
        }

        std::rotate(
            boundary.begin(),
            boundary.begin() + index,
            boundary.end()
        );
        // the first vertex is the contributing vertex

        index = boundary.size() - 1;
        for (int i = boundary.size() - 1; i >= 0; i--) {
            const auto& [edge, flip] = boundary[i];
            std::shared_ptr<Vertex> vertex;
            if (not flip) {
                vertex = edge->custom_links().a;
            }
            else {
                vertex = edge->custom_links().b;
            }

            if (not vertex->state.has_parameter(_param_bend)) {
                continue;
            }
            if (vertex->state.get_parameter(_param_bend)) {
                index = i;
                break;
            }
        }

        for (std::size_t i = 0; i < boundary.size(); i++) {
            const auto& [edge, flip] = boundary[i];
            std::shared_ptr<Vertex> vertex;
            if (not flip) { vertex = edge->custom_links().a; }
            else          { vertex = edge->custom_links().b; }

            if (not vertex->state.has_parameter(_param_bend)) {
                vertex->state.register_parameter(_param_bend, false);
            }
        }

        // remove the trailing vertices that are non contributing
        boundary.erase(boundary.begin() + index, boundary.end());
        for (const auto& [edge, flip] : boundary) {
            edge->custom_links().a->state.update_parameter(_param_bend, true);
            edge->custom_links().b->state.update_parameter(_param_bend, true);
        }

        if (boundary.empty()) {
            throw std::runtime_error("Empty boundary");
        }

        return boundary;
    }

public:
    BoundaryBendElastic (
        std::string name,
        const DataIO::Config& cfg,
        const Model& model
    )
    :
        Base(name, cfg, model),
        _elasticity(get_as<double>("elasticity", cfg)),
        _right(get_as<double>("from_right", cfg, 0.)),
        _top(get_as<double>("from_top", cfg, 0.)),
        _left(get_as<double>("from_left", cfg, 0.)),
        _bottom(get_as<double>("from_bottom", cfg, 0.)),
        _param_bend(this->_name + "_bending"),
        _boundary_condition()
    {    
        const auto bc = get_as<std::string>("boundary", cfg);
        if (bc == "periodic") {
            _boundary_condition = BC::periodic;
        }
        else if (bc == "open") {
            _boundary_condition = BC::open;
        }
        else if (bc == "fixed") {
            _boundary_condition = BC::fixed;

            _fixed_bc = std::make_pair(
                get_as<double>("bc_prior", cfg),
                get_as<double>("bc_post", cfg)
            );
        }
        else {
            throw std::runtime_error(fmt::format(
                "In WF BoundaryBendElastic, unknown boundary condition `{}'. "
                "Please choose one of the following: "
                "- periodic\n"
                "- open\n"
                "- fixed\n",
                bc
            ));
        }

        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.register_parameter(_param_bend, false);
        }

        update_contribution();

        update_parameters(cfg);
    }

    ~BoundaryBendElastic() {
        for (const auto& vertex : this->_am.vertices()) {
            vertex->state.unregister_parameter(_param_bend);
        }
    }

    void compute_and_set_forces () final {
        const auto boundary = get_boundary();
        std::vector<std::shared_ptr<Vertex>> boundary_vs{};
        boundary_vs.reserve(boundary.size() + 1);
        std::vector<SpaceVec> dx{};
        dx.reserve(boundary.size() + 2);
        std::vector<double> l{};
        l.reserve(boundary.size() + 2);
        std::vector<SpaceVec> dx_ds{};
        dx_ds.reserve(boundary.size() + 2);


        // the first vertex
        const auto& [edge, flip] = boundary.front();

        if (not flip) {
            boundary_vs.push_back(edge->custom_links().a);
        }
        else {
            boundary_vs.push_back(edge->custom_links().b);
        }

        /// The junctions
        // the consecutive vertices and junctions
        for (const auto& [edge, flip] : boundary) {
            if (not flip) {
                boundary_vs.push_back(edge->custom_links().b);
            }
            else {
                boundary_vs.push_back(edge->custom_links().a);
            }

            SpaceVec displ = this->_am.displacement(edge, flip);
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
            dx_ds.push_back(displ/arma::norm(displ));
        }


        // the boundary condition
        if (_boundary_condition == BC::periodic) {
            // the prior junction
            auto [edge, flip] = boundary.back();

            SpaceVec displ = this->_am.displacement(edge, flip);
            dx.insert(dx.begin(), displ);
            l.insert(l.begin(), arma::norm(displ));
            dx_ds.insert(dx_ds.begin(), displ/arma::norm(displ));


            // the post junction
            std::tie(edge, flip) = boundary.front();

            displ = this->_am.displacement(edge, flip);
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
            dx_ds.push_back(displ/arma::norm(displ));
        }
        else if (_boundary_condition == BC::open) {
            // the prior junction
            auto [edge, flip] = boundary.front();

            SpaceVec displ = this->_am.displacement(edge, flip);
            dx.insert(dx.begin(), displ);
            l.insert(l.begin(), arma::norm(displ));
            dx_ds.insert(dx_ds.begin(), displ/arma::norm(displ));


            // the post junction
            std::tie(edge, flip) = boundary.back();

            displ = this->_am.displacement(edge, flip);
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
            dx_ds.push_back(displ/arma::norm(displ));
        }
        else if (_boundary_condition == BC::fixed) {
            // the prior junction
            auto [prior, post] = *_fixed_bc__origin;

            SpaceVec displ = this->_am.position_of(boundary_vs.front()) - prior;
            if (arma::norm(displ) < 0.5) {
                prior -= SpaceVec({cos(std::get<0>(_fixed_bc)),
                                   sin(std::get<0>(_fixed_bc))});
                *_fixed_bc__origin = std::make_pair(prior, post);
            }
            dx.insert(dx.begin(), displ);
            l.insert(l.begin(), arma::norm(displ));
            dx_ds.insert(dx_ds.begin(), displ/arma::norm(displ));


            // the post junction
            displ = post - this->_am.position_of(boundary_vs.back());
            if (arma::norm(displ) < 0.5) {
                post += SpaceVec({cos(std::get<1>(_fixed_bc)),
                                  sin(std::get<1>(_fixed_bc))});
                *_fixed_bc__origin = std::make_pair(prior, post);
            }
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
            dx_ds.push_back(displ/arma::norm(displ));
        }


        // The gradients
        const std::size_t N = boundary_vs.size();
        for (std::size_t i = 0; i < N; i++) {
            SpaceVec dE_dx({
                  std::pow(dx_ds[i+1][1] - dx_ds[i][1], 2)
                  * (dx_ds[i+1][0] - dx_ds[i][0])
                  / std::pow(l[i+1] + l[i], 3)
                +   std::pow(dx_ds[i+1][0] - dx_ds[i][0], 3)
                  / std::pow(l[i+1] + l[i], 3)
                + (  dx_ds[i  ][1] * dx_ds[i  ][0] / l[i]
                   + dx_ds[i+1][1] * dx_ds[i+1][0] / l[i+1]  )
                  * (dx_ds[i+1][1] - dx_ds[i][1]) / std::pow(l[i+1] + l[i], 2)
                + (  std::pow(dx_ds[i  ][0], 2) / l[i  ] - 1 / l[i]
                   + std::pow(dx_ds[i+1][0], 2) / l[i+1] - 1 / l[i+1])
                  * (dx_ds[i+1][0] - dx_ds[i][0])
                  / std::pow(l[i+1] + l[i], 2)
                , 
                  std::pow(dx_ds[i+1][0] - dx_ds[i][0], 2) 
                  * (dx_ds[i+1][1] - dx_ds[i][1]) 
                  / std::pow(l[i+1] + l[i], 3)
                + std::pow(dx_ds[i+1][1] - dx_ds[i][1], 3) 
                  / std::pow(l[i+1] + l[i], 3)
                + (  dx_ds[i  ][0] * dx_ds[i  ][1] / l[i]
                   + dx_ds[i+1][0] * dx_ds[i+1][1] / l[i+1])
                  * (dx_ds[i+1][0] - dx_ds[i][0]) / std::pow(l[i+1] + l[i], 2)
                + (  std::pow(dx_ds[i  ][1], 2)/l[i  ] - 1 / l[i] 
                   + std::pow(dx_ds[i+1][1], 2)/l[i+1] - 1 / l[i+1])
                  * (dx_ds[i+1][1] - dx_ds[i][1]) 
                  / std::pow(l[i+1] + l[i], 2)
            });

            boundary_vs[i]->state.add_force(-1. * _elasticity * dE_dx);

            if (i < N - 1 or _boundary_condition == BC::periodic) {
                SpaceVec dE_dx__pre({
                    - dx_ds[i+1][0] * std::pow(dx_ds[i+1][1] - dx_ds[i][1], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i+1][0] * std::pow(dx_ds[i+1][0] - dx_ds[i][0], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i+1][1] * dx_ds[i+1][0] * (dx_ds[i+1][1] - dx_ds[i][1]) 
                      / (l[i+1] * std::pow(l[i+1] + l[i], 2))
                    + (1/l[i+1] - std::pow(dx_ds[i+1][0], 2) / l[i+1]) 
                      * (dx_ds[i+1][0] - dx_ds[i][0])
                      / std::pow(l[i+1] + l[i], 2)
                    , 
                    - dx_ds[i+1][1] * std::pow(dx_ds[i+1][0] - dx_ds[i][0], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i+1][1] * std::pow(dx_ds[i+1][1] - dx_ds[i][1], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i+1][0] * dx_ds[i+1][1] * (dx_ds[i+1][0] - dx_ds[i][0])
                      / (l[i+1] * std::pow(l[i+1] + l[i], 2))
                    + (1/l[i+1] - std::pow(dx_ds[i+1][1], 2) / l[i+1]) 
                      * (dx_ds[i+1][1] - dx_ds[i][1]) 
                      / std::pow(l[i+1] + l[i], 2)
                });
                boundary_vs[(i+1) % N]->state.add_force(
                    -1 * _elasticity * dE_dx__pre
                );
            }
            
            if (i > 0 or _boundary_condition == BC::periodic) {
                SpaceVec dE_dx__post({
                      dx_ds[i][0] * std::pow(dx_ds[i+1][1] - dx_ds[i][1], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    + dx_ds[i][0] * std::pow(dx_ds[i+1][0] - dx_ds[i][0], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i][1]*dx_ds[i][0]*(dx_ds[i+1][1] - dx_ds[i][1]) 
                      / (l[i] * std::pow(l[i+1] + l[i], 2))
                    + (1/l[i] - std::pow(dx_ds[i][0], 2) / l[i])
                      * (dx_ds[i+1][0] - dx_ds[i][0]) 
                      / std::pow(l[i+1] + l[i], 2)
                    ,
                      dx_ds[i][1] * std::pow(dx_ds[i+1][0] - dx_ds[i][0], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    + dx_ds[i][1] * std::pow(dx_ds[i+1][1] - dx_ds[i][1], 2) 
                      / std::pow(l[i+1] + l[i], 3)
                    - dx_ds[i][0]*dx_ds[i][1]*(dx_ds[i+1][0] - dx_ds[i][0]) 
                      / (l[i] * std::pow(l[i+1] + l[i], 2))
                    + (1/l[i] - std::pow(dx_ds[i][1], 2) / l[i]) 
                      * (dx_ds[i+1][1] - dx_ds[i][1]) 
                      / std::pow(l[i+1] + l[i], 2)
                });
                boundary_vs[(i-1+N) % N]->state.add_force(
                    -1 * _elasticity * dE_dx__post
                );
            }
        }
    }

    double compute_energy (
        [[maybe_unused]] const AgentContainer<Vertex>& vertices,
        [[maybe_unused]] const AgentContainer<Edge>& edges,
        [[maybe_unused]] const AgentContainer<Cell>& cells
    ) const final
    {
        double Energy = 0.;

        const auto boundary = get_boundary();
        std::vector<SpaceVec> dx{};
        dx.reserve(boundary.size());
        std::vector<double> l{};
        l.reserve(boundary.size());


        // The junctions
        for (const auto& [edge, flip] : boundary) {
            SpaceVec displ = this->_am.displacement(edge, flip);
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
        }


        // The boundary condition
        if (_boundary_condition == BC::periodic) {
            const auto& [edge, flip] = boundary.front();
            SpaceVec displ = this->_am.displacement(edge, flip);
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
        }
        else if (_boundary_condition == BC::open) {
            // No energy contribution from open boundary condition
        }
        else if (_boundary_condition == BC::fixed) {
            // the prior junction
            auto [prior, post] = *_fixed_bc__origin;

            std::shared_ptr<Vertex> vertex;
            auto [edge, flip] = boundary.front();
            if (not flip) { vertex = edge->custom_links().a; }
            else          { vertex = edge->custom_links().b; }

            SpaceVec displ = this->_am.position_of(vertex) - prior;
            if (arma::norm(displ) < 0.5) {
                prior -= SpaceVec({cos(std::get<0>(_fixed_bc)),
                                   sin(std::get<0>(_fixed_bc))});
                *_fixed_bc__origin = std::make_pair(prior, post);
            }
            dx.insert(dx.begin(), displ);
            l.insert(l.begin(), arma::norm(displ));

            // the post junction
            std::tie(edge, flip) = boundary.back();
            if (not flip) { vertex = edge->custom_links().b; }
            else          { vertex = edge->custom_links().a; }

            displ = post - this->_am.position_of(vertex);
            if (arma::norm(displ) < 0.5) {
                post += SpaceVec({cos(std::get<1>(_fixed_bc)),
                                  sin(std::get<1>(_fixed_bc))});
                *_fixed_bc__origin = std::make_pair(prior, post);
            }
            dx.push_back(displ);
            l.push_back(arma::norm(displ));
        }


        // Calculate the energy
        for (std::size_t i = 1; i < dx.size(); i++) {
            SpaceVec d2x_ds2 = (dx[i]/l[i] - dx[i-1]/l[i-1]) / (l[i] + l[i-1]);

            Energy += 0.5 * _elasticity * std::pow(arma::norm(d2x_ds2), 2);
        }

        return Energy;
    }

    void update_parameters (const DataIO::Config& cfg) final
    {
        _elasticity = get_as<double>("elasticity", cfg, _elasticity);

        if (_boundary_condition == BC::fixed) {
            const auto [prior, post] = _fixed_bc;
            _fixed_bc = std::make_pair(
                get_as<double>("bc_prior", cfg, prior),
                get_as<double>("bc_post", cfg, post)
            );
        }


        _fixed_bc__origin = nullptr;

        if (_boundary_condition == BC::fixed) {
            const auto boundary = get_boundary();
            
            std::shared_ptr<Vertex> vertex__start;
            std::shared_ptr<Vertex> vertex__end;

            auto [edge, flip] = boundary.front();
            if (not flip) {
                vertex__start = edge->custom_links().a;
            }
            else {
                vertex__start = edge->custom_links().b;
            }

            std::tie(edge, flip) = boundary.back();
            if (not flip) {
                vertex__end = edge->custom_links().b;
            }
            else {
                vertex__end = edge->custom_links().a;
            }

            auto [prior, post] = _fixed_bc;
            SpaceVec displ_prior({cos(prior), sin(prior)});
            SpaceVec displ_post ({cos(post) , sin(post) });

            _fixed_bc__origin = std::make_shared<std::pair<SpaceVec, SpaceVec>>(
                this->_am.position_of(vertex__start) - displ_prior,
                this->_am.position_of(vertex__end)   + displ_post
            );
        }
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
