"""Tests of the configurations for the PCPTopology"""

import pytest

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
                    'num_generations': 0,
                    'iterations_prolog': 0
                }
            }
        }
    })
    debug_params=dict({
        'polarity': {
            'num_steps': 5,
            'PCPTopology': {
                'initialisation': {
                    '1_by_proliferation': {
                        'enabled': False,
                        'num_generations': 0
                    }
                },
                'setup_params': {
                    'hexagonal': {
                        'lattice_rows': 6,
                        'lattice_columns': 6
                    }
                }
            }
        }
    })

    for cfg_name, cfg_paths in mtc.default_config_sets.items():
        print("\nRunning '{}' example ...".format(cfg_name))

        if cfg_name != 'evaluate_diffusion':
            continue

        params = debug_params.get(cfg_name, _params)

        mv, _ = mtc.create_run_load(from_cfg=cfg_paths.get('run'),
                                    plot_manager={'raise_exc': True},
                                    parameter_space=params
                                    )
        mv.pm.plot_from_cfg(plots_cfg=cfg_paths.get('eval'))

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