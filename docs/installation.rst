Installation
================

The AntPack package is not yet provided on PyPi. For now, to install, check
the `github repo <https://github.com/jlparkI/AntPack>`_ and look at the releases
tab. Find the number for the latest release. If it is for example 0.1.1, then
type:::

  pip install git+https://github.com/jlparkI/AntPack@0.1.1

The package contains C++ code which should be automatically compiled on
install. If you encounter any issues please contact us using the About / Contact
section on the main page.

Required dependencies are ``pybind11``, ``numpy``,
``biopython``, ``setuptools``, ``wheel`` and ``scipy``. These should be
automatically installed but if desired can be installed first. To minimize the
number of dependencies as much as possible, we plan to eliminate the
scipy dependency in a future version. To further facilitate ease of installation
we plan to move to distributing on PyPi as well.

AntPack is also available as a docker image on `DockerHub <https://hub.docker.com/r/jlparkinson1/antpack>`_.
The docker image when run will provide the output of the package unit tests.
We do not encourage the use of the docker image in general since the pip installation
procedure shown above is easier and less cumbersome, but it is available for users
who are interested.
