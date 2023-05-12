"""Tests of the configurations for the PCPTopology"""

from copy import deepcopy

import signal
import time

import pytest

from utopya.tools import recursive_update
from utopya.testtools import ModelTest

# Configure the ModelTest class
mtc = ModelTest("PCPTopology", test_file=__file__)

# Fixtures --------------------------------------------------------------------


# Tests -----------------------------------------------------------------------

@pytest.mark.parametrize("cfg_name", mtc.default_config_sets.keys())
def test_run_and_eval_cfgs(cfg_name):
    """Carries out all additional configurations that were specified alongside
    the default model configuration.

    This is done automatically for all run and eval configuration pairs that
    are located in subdirectories of the ``cfgs`` directory (at the same level
    as the default model configuration).
    If no run or eval configurations are given in the subdirectories, the
    respective defaults are used.

    See :py:meth:`~utopya.model.Model.default_config_sets` for more info.
    """
    _params=dict({
        'num_steps': 3,
        'PCPTopology': {
            'initialisation': {
                '1_by_proliferation': {
                    'enabled': False,
                    'num_generations': 0,
                    'iterations_prolog': 0
                }
            },
        }
    })
    debug_params=dict({
        'area_graded': {
            'num_steps': 1,
            'PCPTopology': {
                'PCPVertex': {
                    'agent_manager': {
                        'setup_params': {
                            'hexagonal': {
                                'lattice_rows': 6,
                                'lattice_columns': 8
                            }
                        }
                    }
                }
            }
        },
        'asymmetries': {
            'num_steps': 1,
            'PCPTopology': {
                'PCPVertex': {
                    'agent_manager': {
                        'setup_params': {
                            'hexagonal': {
                                'lattice_rows': 10,
                                'lattice_columns': 20
                            }
                        }
                    }
                }
            }
        },
        'columnar': {
            'num_steps': 1,
            'PCPTopology': {
                'PCPVertex': {
                    'agent_manager': {
                        'setup_params': {'column': {'num_cells': 6}}
                    }
                }
            }
        },
        'evaluate_shear_thinning': {
            'num_steps': 1,
            'PCPTopology': {
                'initialisation': {
                    '1_by_proliferation': {
                        'enabled': True,
                        'num_generations': 1,
                        'iterations_prolog': 1,
                        'minimization': {'num_repeat': 24}
                    }
                },
            }
        },
        'heterogeneities_and_polarity': {
            'PCPTopology': {
                'PCPVertex': {
                    'agent_manager': {
                        'setup_params': {
                            'hexagonal': {
                                'lattice_rows': 4,
                                'lattice_columns': 4
                            }
                        }
                    }
                }
            }
        },
        'polarity': {
            'num_steps': 5,
            'PCPTopology': {
                'PCPVertex': {
                    'agent_manager': {
                        'setup_params': {
                            'hexagonal': {
                                'lattice_rows': 6,
                                'lattice_columns': 6
                            }
                        }
                    }
                }
            }
        },
        'proliferation_minimal__by_generation': {
            'num_steps': 2
        }
    })

    cfg_paths = mtc.default_config_sets[cfg_name]

    params = recursive_update(deepcopy(_params),
                              deepcopy(debug_params.get(cfg_name, dict())))

    print("  With update config {}".format(params))
    
            
    # Run multiverse with a timeout of 120 seconds
    timeout = 120
    def timeout_handler(signum, frame):
        print("Aborting test of '{}' after {} seconds. "
                "Timeout reached!".format(cfg_name, timeout))
        raise RuntimeError(
            "Aborting test of '{}' after {} seconds. "
            "Timeout reached!".format(cfg_name, timeout))

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)

    
    mv, _ = mtc.create_run_load(from_cfg=cfg_paths.get('run'),
                                plot_manager={'raise_exc': True},
                                parameter_space=params
                                )

    # Evaluate multiverse with a timeout of 60 seconds
    timeout = 120
    def timeout_handler(signum, frame):
        print("Aborting test of '{}' after {} seconds. "
                "Timeout reached!".format(cfg_name, timeout))
        raise RuntimeError(
            "Aborting test of '{}' after {} seconds. "
            "Timeout reached!".format(cfg_name, timeout))

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)
    
    mv.pm.plot_from_cfg(plots_cfg=cfg_paths.get('eval'))


    signal.alarm(0)

    print("Succeeded running and evaluating '{}'.\n".format(cfg_name))


def test_disabled_plots():
    mv, dm = mtc.create_run_load(
        plot_only=[
            "equilibrium_cellular_structure__area_elasticity",
            "equilibrium_cellular_structure__cell_contractility",
            "equilibrium_cellular_structure__linetension",
            "equilibrium_cellular_structure__edge_contractility",
            # "equilibrium_cellular_structure__hexatic_order"
        ]
    )