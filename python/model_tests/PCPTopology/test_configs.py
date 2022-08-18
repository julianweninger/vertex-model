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

def test_run_and_eval_cfgs():
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
        'evaluate_diffusion': {
            'num_seeds': 2
        },
        'evaluate_shear_thinning': {
            'num_seeds': 2,
            'num_steps': 1
        },
        'increment_area': {
            'num_seeds': 2,
        },
        'increment_area_graded': {
            'num_seeds': 2,
            'num_steps': 1
        },
        'heterogeneities_and_polarity': {
            'num_seeds': 2,
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
        'proliferation_minmal__by_generation': {
            'num_steps': 2
        }
    })

    print("Running the following tests: ")
    for cfg_name, cfg_paths in mtc.default_config_sets.items():
        print(cfg_name)

    for cfg_name, cfg_paths in mtc.default_config_sets.items():
        print("\nRunning '{}' example ...".format(cfg_name))

        params = recursive_update(deepcopy(_params),
                                  deepcopy(debug_params.get(cfg_name, dict())))

        print("With updated config {}", params)
        
        class AlarmException(Exception):
            pass

        def alarmHandler(signum, frame):
            raise AlarmException

        # Raise timeout after 120 seconds
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(120)
        
        try:
            mv, _ = mtc.create_run_load(from_cfg=cfg_paths.get('run'),
                                        plot_manager={'raise_exc': True},
                                        parameter_space=params
                                        )
        except AlarmException:
            print("Aborting test of '{}' after 2 minutes.")
            raise AlarmException
        
        signal.alarm(60)
        
        try:
            mv.pm.plot_from_cfg(plots_cfg=cfg_paths.get('eval'))
        except AlarmException:
            print("Aborting plotting-test of '{}' after 2 minutes.")
            raise AlarmException

        signal.alarm(0)

        print("Succeeded running and evaluating '{}'.\n".format(cfg_name))


def test_disabled_plots():
    mv, dm = mtc.create_run_load(
        plot_only=[
            "equilibrium_cellular_structure__area_elasticity",
            "equilibrium_cellular_structure__cell_contractility",
            "equilibrium_cellular_structure__linetension",
            "equilibrium_cellular_structure__edge_contractility",
            "equilibrium_cellular_structure__hexatic_order",
            "equilibrium_cellular_structure__rotation"
        ]
    )