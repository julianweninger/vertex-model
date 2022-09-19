#ifndef UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH
#define UTOPIA_MODELS_PCPVERTEX_ENTITIES_HH

#include "utils.hh"

namespace Utopia::Models::PCPVertex {

/// The state of a Vertex
class VertexState {
public:
    using SpaceVec = Utopia::SpaceVecType<2>;

private:
    /// The steepest gradient slope
    SpaceVec f;

public:
    /// Whether to fix the position of this vertex, i.e. boundary condition
    bool fix_in_space;

    /// Whether to remove
    bool remove;

    // friend class EntitiesManager;

public:
    /// Constructor
    VertexState ()
    :  
        f(),
        fix_in_space(false),
        remove(false)
    {
        f.zeros();
    }

    inline const SpaceVec& get_force () const {
        return f;
    }

    inline void add_force (const SpaceVec& force) {
        f += force;
    }

    inline void reset_force () {
        f = SpaceVec({0., 0.});
    }
};


/// The Edge object
/** It is defined as the shortest connection between two vertices a and b
 * 
 *  Every edge has two adjoint cells, one to either side, except boundary edges
 *  have only one.
 */
class EdgeState {
private:
    std::map<std::string, std::vector<double>> _parameters;

public:
    double tension;

    /// The time T1 transition was last attempted
    /** 0 if never attempted */
    std::size_t last_T1_attempt;

    /// Whether this object is to be removed 
    bool remove;

    // friend class EntitiesManager;

public:
    /// Constructor
    EdgeState ()
    :
        _parameters({}),
        tension(0.),
        last_T1_attempt(0),
        remove(false)
    { }

public:
    // getters and setters ....................................................
    const std::unordered_set<std::string> list_parameters () const {
        std::unordered_set<std::string> keys({});
        for(const auto& kv : _parameters) {
            keys.insert(kv.first);
        }
        return keys;
    }

    bool has_parameter (const std::string& name) const {
        return _parameters.find(name) != _parameters.end();
    }

    const auto& get_parameter (const std::string& name) const {
        if (_parameters.find(name) == _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot find EdgeState parameter with name `{}`", name));
        }
        return _parameters.at(name);
    }

    /// Register parameter
    void register_parameter (const std::string& name,
                             const std::vector<double>& values) {
        if (_parameters.find(name) != _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot register EdgeState parameter with name `{}`! "
                "A parameter with such a name is already registered.", name));
        }
        _parameters[name] = values;
    }

    /// Remove parameter from register
    void unregister_parameter (const std::string& name) {
        _parameters.erase(name);
    }

    /// Update value of parameter
    void update_parameter (const std::string& name,
                           const std::vector<double>& values) {
        if (_parameters.find(name) == _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot find EdgeState parameter registered with name `{}`. "
                "Please register the parameter first.", name));
        }
        _parameters[name] = values;
    }

    void inherit_from_state (const EdgeState& state) {
        for (const auto& k : state.list_parameters()) {
            this->register_parameter(k, state.get_parameter(k));
        }
    }
};


/// The Cell defined by its id, its vertices and its area
struct CellState {
public:
    using SpaceVec = Utopia::SpaceVecType<2>;

private:
    std::map<std::string, std::vector<double>> _parameters;

public:
    double pressure;

    /// The type of a cell
    std::size_t type;

    /// Whether this object is to be removed 
    bool remove;

    // friend class EntitiesManager;

public:
    /// Constructor of a cell
    /** Construct a cell from a configuration
     * 
     *  \param cfg      The configuration
     *      -  `type` (uint, default: 0): The type of cell.
     */
    CellState (const Config& cfg)
    :
        _parameters({}),
        pressure(0.),
        type(get_as<std::size_t>("type", cfg, 0)),
        remove(false)
    { }
  
public:
    // getters and setters ....................................................
    const std::unordered_set<std::string> list_parameters () const {
        std::unordered_set<std::string> keys({});
        for(const auto& kv : _parameters) {
            keys.insert(kv.first);
        }
        return keys;
    }

    bool has_parameter (const std::string& name) const {
        return _parameters.find(name) != _parameters.end();
    }

    const auto get_parameter (const std::string& name) const {
        if (_parameters.find(name) == _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot find EdgeState parameter with name `{}`", name));
        }
        return _parameters.at(name);
    }

    /// Register parameter
    void register_parameter (const std::string& name,
                             const std::vector<double>& values) {
        if (_parameters.find(name) != _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot register EdgeState parameter with name `{}`! "
                "A parameter with such a name is already registered.", name));
        }
        _parameters[name] = values;
    }

    /// Remove parameter from register
    void unregister_parameter (const std::string& name) {
        _parameters.erase(name);
    }

    /// Update value of parameter
    void update_parameter (const std::string& name,
                           const std::vector<double>& values) {
        if (_parameters.find(name) == _parameters.end()) {
            throw std::runtime_error(fmt::format(
                "Cannot find EdgeState parameter registered with name `{}`. "
                "Please register the parameter first.", name));
        }
        _parameters[name] = values;
    }

    void inherit_from_state (const CellState& state) {
        this->type = state.type;
        for (const auto& k : state.list_parameters()) {
            this->register_parameter(k, state.get_parameter(k));
        }
    }
};

} // namespace Utopia::Models::PCPVertex

#endif