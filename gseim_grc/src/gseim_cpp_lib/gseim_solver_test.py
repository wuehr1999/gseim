import filecmp
import os
import subprocess
import sys

from importlib_resources import files
import pytest

def _test_solver(fname):
    r = subprocess.run(
        [
            str(files('gseim_cpp_lib').joinpath('gseim_solver')),
            str(files('gseim_cpp_lib').joinpath('test_data', 'input', fname)),
        ],
        capture_output=True,
        text=True,
        env={'HOME': 'not set'},
    )

    assert 'Program completed.' in r.stdout
    r.check_returncode()

    print(r.stdout)

    out_fname = os.path.splitext(fname)[0] + '.dat'

    assert filecmp.cmp(
        str(files('gseim_cpp_lib').joinpath('test_data', 'input', out_fname)),
        str(files('gseim_cpp_lib').joinpath('test_data', 'output', out_fname)),
    )

def test_ac_controllerer_1():
    _test_solver('ac_controller_1.in')

def test_buck():
    _test_solver('buck.in')

def test_controlled_rectifier_2():
    _test_solver('controlled_rectifier_2.in')

def test_test_1():
    _test_solver('test_1.in')

def test_test_2():
    _test_solver('test_2.in')

def test_test_3():
    _test_solver('test_3.in')

def test_test_4():
    _test_solver('test_4.in')

def test_test_5():
    _test_solver('test_5.in')

def test_test_6():
    _test_solver('test_6.in')

def test_test_7():
    _test_solver('test_7.in')

def test_test_8():
    _test_solver('test_8.in')

def test_test_9():
    _test_solver('test_9.in')

def test_test_10():
    _test_solver('test_10.in')

def test_test_11():
    _test_solver('test_11.in')

def test_test_12():
    _test_solver('test_12.in')

def test_test_13():
    _test_solver('test_13.in')

def test_test_14():
    _test_solver('test_14.in')

def test_test_15():
    _test_solver('test_15.in')

def test_test_16():
    _test_solver('test_16.in')

def test_test_17():
    _test_solver('test_17.in')

if __name__ == '__main__':
    sys.exit(pytest.main(sys.argv))

