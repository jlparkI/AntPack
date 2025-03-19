Four different ways to use AntPack
=====================================

You can use AntPack from the command line, from Python, from C++,
or using the graphical user interface (GUI).

The Python interface is extensive, flexible, easy to use
and is the most powerful way to use AntPack for most data analysis.
If you have a small-medium size dataset, however, and just want some basic
information about your sequences, starting with AntPack
v0.3.5 you can also use the command line interface (CLI).

To run the CLI, type:::

  AntPack-CLI

on your command line for a list of options. The CLI generates output csv
files in a simple csv file format containing all your input sequences as
numbered multiple sequence alignments; it also optionally assigns VJ genes.

For more extensive analysis we suggest using the Python API, which is
described in the remaining documentation.

Starting in v0.3.7, you can also use the GUI. The GUI makes it easy to quickly
view and compare a few sequences against each other or against their closest
VJ genes in humans or mice. To launch the GUI type:::

  AntPack-GUI
