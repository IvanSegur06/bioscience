python3 setup.py clean
python3 setup.py sdist bdist_wheel

# Local installation
python3 -m pip install ./dist/bioscience-0.1.4.tar.gz

# Upload to PyPi repository
#twine upload dist/*