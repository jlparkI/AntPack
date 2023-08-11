"""Builds the package."""
import os
import sys

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("Python version >= 3.7 required.")

import platform
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages, Extension
from subprocess import getoutput



def get_version(setup_fpath):
    """Retrieves the version number."""

    os.chdir(os.path.join(setup_fpath, "antpack"))
    with open("__init__.py", "r") as fhandle:
        version_line = [l for l in fhandle.readlines() if
                    l.startswith("__version__")]
        version = version_line[0].split("=")[1].strip().replace('"', "")
    os.chdir(setup_fpath)
    return version


def using_clang():
    """Check to see if we are using Clang."""
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    return "clang" in compiler_ver



if __name__ == "__main__":
    setup_fpath = os.path.dirname(os.path.abspath(__file__))
    setup(
        name="antpack",
        version=get_version(setup_fpath),
        description="A Python package for processing, manipulating and making inferences about antibody sequence data",
        long_description="A Python package for processing, manipulating and making inferences about antibody sequence data",
        packages=find_packages(),
        setup_requires=['pybind11>=2.4'],
        install_requires=['pybind11>=2.4', "numpy", "scipy", "pyhmmer"],
        include_package_data=True,
        license="MIT",
        python_requires=">=3.7",
    )
