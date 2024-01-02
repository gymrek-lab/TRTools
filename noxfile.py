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


cov_cli_args = [
    "--cov=.",
    "--cov-report=term-missing",
]


def install_handle_python(session):
    """
    handle incompatibilities between python and other packages
    see https://github.com/cjolowicz/nox-poetry/issues/1116
    """
    # install the latest versions of all dependencies for py3.9+
    # but install the locked versions for < py3.9 to ensure some stability in the CI and
    # help us understand when we need to bump our version constraints
    if session._session.python in ("3.9", "3.10", "3.11"):
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
            "art",
            "bcftools==1.19",
            channel="bioconda",
        )
        session.conda_install(
            "numpy",
            "pytest",
            "pytest-cov",
            channel="conda-forge",
        )
        install_handle_python(session)
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
