#ifndef UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH
#define UTOPIA_MODELS_PCPVERTEX_BOUNDARY_HH

#include "work_function.hh"


namespace Utopia {
namespace Models {
namespace PCPVertex {
namespace WorkFunction {

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
        const SpaceVec up_corner = ll_corner + get_extent();

        SpaceVec displ({0., 0.});
        for (std::size_t axis = 0; axis <= 1; axis++) {
            if (pos[axis] < ll_corner[axis]) {
                displ[axis] = ll_corner[axis] - pos[axis];
            }
            else if (pos[axis] > up_corner[axis]) {
                displ[axis] = up_corner[axis] - pos[axis];
            }
        }

        return displ;
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

        if (get_as<bool>("deform_plastic", cfg)) {
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
    }
};


} // namespace WorkFunction
} // namespace PCPVertex
} // namespace Models
} // namespace Utopia

#endif
