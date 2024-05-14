Installation
================

Starting with version 0.1.5, AntPack is available on PyPi. To install it,
type:::

  pip install antpack

Also see the `github repo <https://github.com/jlparkI/AntPack>`_
Note that since AntPack is distributed as an sdist,
there are portions of the code that have to be compiled on install. Pip should
ordinarily handle this for you automatically. If you encounter any issues please
contact us using the About / Contact section on the main page.

Required dependencies are ``pybind11`` and ``numpy``. These should be
automatically installed but if desired can be installed first. Older versions also require Biopython
but starting with 0.2.0 this is no longer needed.

AntPack is also available as a docker image on `DockerHub <https://hub.docker.com/r/jlparkinson1/antpack>`_.
The docker image when run will provide the output of the package unit tests.
We do not encourage the use of the docker image in general, but it is available for
users who are interested. We may discontinue support for the Docker image in
a future release.
