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
    for cfg_name, cfg_paths in mtc.default_config_sets.items():
        print("\nRunning '{}' example ...".format(cfg_name))

        mv, _ = mtc.create_run_load(from_cfg=cfg_paths.get('run'),
                                    plot_manager={'raise_exc': True})
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