"""Nox sessions."""
import os
from pathlib import Path

import nox
from nox_poetry import Session
from nox_poetry import session


package = "trtools"
python_versions = ["3.7", "3.8", "3.9", "3.10", "3.11"]
nox.needs_version = ">= 2022.11.21"
nox.options.sessions = (
    "tests",
)

# what percentage of the code must be covered?
# set to None to disable coverage testing
COVERAGE_REQUIREMENT = 94

cov_cli_args = []
if COVERAGE_REQUIREMENT is not None:
    cov_cli_args = [
        "--cov=.",
        "--cov-report=term-missing",
        "--cov-fail-under=" + str(COVERAGE_REQUIREMENT),
    ]


def install_handle_python_numpy(session):
    """
    handle incompatibilities with python and numpy versions
    see https://github.com/cjolowicz/nox-poetry/issues/1116
    """
    if session._session.python == "3.11":
        session._session.install(".")
    else:
        session.install(".")


# detect whether conda/mamba is installed
if os.getenv("CONDA_EXE"):
    conda_cmd = "conda"
    if (Path(os.getenv("CONDA_EXE")).parent / "mamba").exists():
        conda_cmd = "mamba"
    conda_args = ["-c", "conda-forge"]

    @session(venv_backend=conda_cmd, venv_params=conda_args, python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.conda_install(
            "numpy",
            "pytest",
            "pytest-cov",
            channel="conda-forge",
        )
        session.conda_install(
            "art",
            "bcftools",
            channel="bioconda",
        )
        install_handle_python_numpy(session)
        session.run(
           "python", "-m", "pytest", *cov_cli_args, *session.posargs
        )
        session.run("./test/cmdline_tests.sh", external=True)

else:

    @session(python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.install("pytest", "pytest-cov")
        install_handle_python_numpy(session)
        session.run(
            "python", "-m", "pytest", *cov_cli_args, *session.posargs
        )
        session.run("./test/cmdline_tests.sh", external=True)
