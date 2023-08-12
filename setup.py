"""The package setup file for AntPack."""
import os
import sys
from subprocess import getoutput
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages, Extension

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("Python version >= 3.7 required.")




def get_version(setup_fpath):
    """Retrieves the version number."""

    os.chdir(os.path.join(setup_fpath, "antpack"))
    with open("__init__.py", "r", encoding="utf-8") as fhandle:
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
    home_dir = os.path.dirname(os.path.abspath(__file__))
    setup(
        name="antpack",
        version=get_version(home_dir),
        description="A Python package for processing, manipulating and making inferences about antibody sequence data",
        long_description="A Python package for processing, manipulating and making inferences about antibody sequence data",
        packages=find_packages(),
        setup_requires=['pybind11>=2.4'],
        install_requires=['pybind11>=2.4', "numpy", "pyhmmer"],
        include_package_data=True,
        license="MIT",
        python_requires=">=3.7",
    )
