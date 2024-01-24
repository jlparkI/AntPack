pip uninstall AntPack
rm -rf build antpack.egg-info
rm *.so
python setup.py develop
