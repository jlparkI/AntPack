Installation
================

The AntPack package is not yet provided on PyPi. For now, to install, check
the `github repo <https://github.com/jlparkI/AntPack>`_ and look at the releases
tab. Find the number for the latest release. If it is for example 0.1.0, then
type:::

  pip install git+https://github.com/jlparkI/AntPack@0.1.0

which will install version 0.0.1. Later, once more tools are added, we may move
to releasing on PyPi as well. Required dependencies are ``pybind11``, ``numpy``,
``biopython``, ``setuptools``, ``wheel`` and ``scipy``. These should be
automatically installed but if desired can be installed first. We may
be able to eliminate the scipy dependency in the near future.
