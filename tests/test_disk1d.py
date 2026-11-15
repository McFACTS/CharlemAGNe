#!/usr/bin/env python3
######## Imports ########
#### Standard Library ####
import os
import time
#### Local ####
from charlemagne.disk import SG03, TQM05, SG03_DEFAULTS, TQM05_DEFAULTS
######## Tests ########
def test_sg03_obj():
    tic = time.perf_counter()
    sg03 = SG03()
    toc = time.perf_counter()
    print(f"SG03 time: {toc-tic:.6f}s")
    # Test params
    for key, value in SG03_DEFAULTS.items():
        assert key in sg03.params
        assert getattr(sg03, key) == value
    # Test that profiles exist
    assert sg03["r_cyl"].size == SG03_DEFAULTS["npts"]
    assert "kappa" in sg03.fields

def test_tqm05_obj():
    tic = time.perf_counter()
    tqm05 = TQM05()
    toc = time.perf_counter()
    print(f"TQM05 time: {toc-tic:.6f}s")
    # Test params
    for key, value in TQM05_DEFAULTS.items():
        assert key in tqm05.params
        assert getattr(tqm05, key) == value
    # Test that profiles exist
    assert tqm05["r_cyl"].size == TQM05_DEFAULTS["npts"]
    assert "kappa" in tqm05.fields

######## Execution ########
if __name__ == "__main__":
    test_sg03_obj()
    test_tqm05_obj()
