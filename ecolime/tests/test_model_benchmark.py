from __future__ import division, absolute_import, print_function

import pytest

from os.path import dirname, join, abspath

from cobrame.io.jsonme import load_json_me
from ecolime.build_ME_model import return_ME_model
from ecolime.util.ME_model_comparison import find_ME_model_difference

current_dir = dirname(abspath(__file__))
models_dir = current_dir.split('tests')[0] + 'me_models/'

del dirname, abspath


def test_model_benchmark():
    benchmark_model = load_json_me(join(models_dir, 'iLE1678_benchmark.json'))
    test_model = return_ME_model()
    difference = find_ME_model_difference(benchmark_model, test_model, 1e-6)
    print('-----------------------Difference----------------------')
    print(difference)
    assert (len(difference) == 0)

if __name__ == '__main__':
    test_model_benchmark()
