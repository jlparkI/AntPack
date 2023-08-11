"""The package setup file for AntPack."""
import os
from setuptools import setup, find_packages


setup_fpath = os.path.abspath(os.path.dirname(__file__))


read_me = os.path.join(setup_fpath, "README.md")
with open(read_me, "r", encoding="utf-8") as fhandle:
    long_description = "".join(fhandle.readlines())


setup(name='antpack',
     version='0.0.1',
     description='A toolkit for antibody data processing and inference',
     packages=find_packages(),
     data_files = ['dependencies/muscle', 'dependencies/muscle_macOS'],
     long_description = long_description,
     long_description_content_type="text/markdown",
     include_package_data = True,
     scripts=[],
    )
