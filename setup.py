"""The package setup file for AntPack."""
import os
import sys
import platform
from subprocess import getoutput
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from setuptools import setup, find_packages, Extension
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext

if sys.version_info[:2] < (3, 7):
    raise RuntimeError("Python version >= 3.7 required.")


def using_clang():
    """Check to see if we are using Clang."""
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    return "clang" in compiler_ver




def get_version(setup_fpath):
    """Retrieves the version number."""

    os.chdir(os.path.join(setup_fpath, "antpack"))
    with open("__init__.py", "r", encoding="utf-8") as fhandle:
        version_line = [l for l in fhandle.readlines() if
                    l.startswith("__version__")]
        version = version_line[0].split("=")[1].strip().replace('"', "")
    os.chdir(setup_fpath)
    return version



def main():
    """Builds the package and extension."""
    home_dir = os.path.dirname(os.path.abspath(__file__))
    read_me = os.path.join(home_dir, "README.md")
    with open(read_me, "r", encoding="utf-8") as fhandle:
        long_description = "".join(fhandle.readlines())

    cpp_extra_link_args = []
    cpp_extra_compile_args = [
        "-std=c++11",
        "-O3"
    ]

    # Mac-specific options
    if platform.system() == "Darwin" and using_clang():
        cpp_extra_compile_args.append("-stdlib=libc++")
        cpp_extra_compile_args.append("-mmacosx-version-min=10.9")
        cpp_extra_link_args.append("-stdlib=libc++")
        cpp_extra_link_args.append("-mmacosx-version-min=10.7")

    extensions=[
        Pybind11Extension("ant_ext",
            sources=[
                "antpack/ext/ant_ext.cpp",
                "antpack/ext/aligners.cpp"
            ],
            include_dirs=[
                "antpack/ext",
                pybind11.get_include(),
            ],
            language="c++",
            extra_compile_args=cpp_extra_compile_args  + ["-fvisibility=hidden"], # needed by pybind
            extra_link_args=cpp_extra_link_args,
        )
    ]

    setup(
        name="antpack",
        version=get_version(home_dir),
        description="A Python package for processing, manipulating and making inferences about antibody sequence data",
        long_description=long_description,
        packages=find_packages(),
        cmdclass={"build_ext":build_ext},
        setup_requires=['pybind11>=2.4'],
        install_requires=['pybind11>=2.4', "numpy", "pyhmmer"],
        include_package_data=True,
        license="MIT",
        ext_modules=extensions,
        python_requires=">=3.7",
    )





if __name__ == "__main__":
    main()
