Installation
================

Starting with v0.3.9, AntPack can only be used for non-commercial
purposes and is made available under an
ASL license. To use it you'll first need to obtain a license key
`at this link. <https://pwslicensekey.pythonanywhere.com/>`_

The only required dependency is numpy which is automatically installed.
If you want to run the GUI, however, you have to install two additional
optional dependencies:::

  pip install pyside6 qt_material

If you don't want to use the GUI these dependencies are not needed. These
dependencies may no longer be required in future.

To install AntPack, type:::

  pip install antpack

Then:::

  AntPack-setup

AntPack will prompt you for your license key and email then verify these.
These will be saved so you do not need to enter them again, although
you should keep them in case you need to reinstall AntPack.

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
Python3.9 or later is required.
