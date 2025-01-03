Installation
================

The only required dependency is numpy which is automatically installed.
If you want to run the GUI, however, you have to install two additional
optional dependencies:::

  pip install pyside6 qt_material

If you don't want to use the GUI these dependencies are not needed. These
dependencies may no longer be required in future.

To install AntPack, type:::

  pip install antpack

Also see the `github repo <https://github.com/jlparkI/AntPack>`_
AntPack is distributed precompiled as a wheel
for most platforms. In these cases pip will install it without any need to
compile any code. In the event of any issues, antpack is still available as a
source distribution, and pip can still install
easily -- pip will automatically try to find a C++ compiler on your system and
compile from source (C++17 or later only). In the unlikely event that you encounter
any issues please contact us using the Issues page on Github (see About / Contact
section on the main page).

Required dependencies are ``numpy``. These should be
automatically installed but if desired can be installed first.
Python3.8 or later is required.

Licensing
===========

AntPack is licensed under a reasonably permissive open-source license,
so that you are free to use it for your own data analysis in any manner you like
regardless of whether you are working on academic research or industrial R&D.
Note that under the terms of the license however you cannot use AntPack
to build closed-source software that you plan to sell -- if you use AntPack
to build software for distribution the software must also be open-source, must include
the GPL license, and must acknowledge AntPack appropriately. If you are interested
in using AntPack in a closed-source application, please
`contact us <https://mapbioscience.com/contact/>`_ and we can set you up with a version
of AntPack licensed for closed-source distribution.
