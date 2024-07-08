"""Tests of the configurations for the PCPVertex"""

from copy import deepcopy

import signal
import time

import pytest

from utopya.tools import recursive_update
from utopya.testtools import ModelTest

# Configure the ModelTest class
mtc = ModelTest("CopyMeVertex", test_file=__file__)

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
        'num_steps': 3
    })
    debug_params=dict({

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
