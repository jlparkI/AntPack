AntPack APIs
================

You can use AntPack from the command line, from Python, or
from C++. The C++ interface is mainly for building software.

The Python interface is extensive, flexible, easy to use
and is the most appropriate way to use AntPack for most data analysis.
If you have a small-medium size dataset, however, and just want some basic
information about your sequences, starting with AntPack
v0.3.5 you can also use the command line interface (CLI).

To run the CLI, type:::

  AntPack

on your command line for a list of options. The CLI generates output csv
files in a format very similar to the ANARCI tool containing all your input
sequences as numbered multiple sequence alignments; it also optionally assigns
VJ genes.

For more extensive analysis we suggest using the Python API, which is
described in the remaining documentation.
