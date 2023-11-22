import os

import pytest

# functions that all pytests have access to
# Not for use in runtime code!

# Usage based on this documentation
# https://docs.pytest.org/en/latest/example/simple.html#pass-different-values-to-a-test-function-depending-on-command-line-options
def pytest_addoption(parser):
	default = os.path.dirname(os.path.abspath(__file__))
	parser.addoption(
			"--datadir",
			default=default,
			help="Directory containing sample data files"
	)

@pytest.fixture()
def vcfdir(request):
	return os.path.join(request.config.getoption("--datadir"), "sample_vcfs")

@pytest.fixture()
def simtrdir(request):
	return os.path.join(request.config.getoption("--datadir"), "sample_simtrdata")

@pytest.fixture()
def regiondir(request):
	return os.path.join(request.config.getoption("--datadir"), "sample_regions")

@pytest.fixture()
def statsdir(request):
	return os.path.join(request.config.getoption("--datadir"),
                     "sample_stats")


