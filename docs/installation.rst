Installation
================

Starting with version 0.1.5, AntPack is available on PyPi. To install it,
type:::

  pip install antpack

Also see the `github repo <https://github.com/jlparkI/AntPack>`_
Starting with version 0.27, AntPack is distributed precompiled as a wheel
for some platforms (Windows 64-bit, Linux x86-64). In these cases pip will
install it without any need to compile any code. If your platform is
not one of these, it is still available as an sdist, and pip can still install
easily -- pip will automatically try to find a C++ compiler on your system and
compile from source. In the unlikely event that you encounter any issues please
contact us using the Issues page on Github (see About / Contact section on the main page).

Required dependencies are ``pybind11`` and ``numpy``. These should be
automatically installed but if desired can be installed first. Older versions also require Biopython
but starting with 0.2.0 this is no longer needed. Python3.6 or later is required.

AntPack is also available as a docker image on `DockerHub <https://hub.docker.com/r/jlparkinson1/antpack>`_.
The docker image when run will provide the output of the package unit tests.
We do not encourage the use of the docker image in general, but it is available for
users who are interested. We may discontinue support for the Docker image in
a future release.
