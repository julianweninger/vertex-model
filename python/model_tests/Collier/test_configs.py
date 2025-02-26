"""Tests of the configurations for the NotchDelta model"""

import pytest

from utopya.testtools import ModelTest

# Configure the ModelTest class
mtc = ModelTest("Collier", test_file=__file__)

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

    cfg_paths = mtc.default_config_sets[cfg_name]
    
    mv, _ = mtc.create_run_load(from_cfg=cfg_paths.get('run'),
                                plot_manager={'raise_exc': True})
    mv.pm.plot_from_cfg(plots_cfg=cfg_paths.get('eval'))

    print("Succeeded running and evaluating '{}'.\n".format(cfg_name))

