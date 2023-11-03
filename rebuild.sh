python3 setup.py clean
python3 setup.py sdist
python3 -m build

# Local installation
python3 -m pip install ./dist/bioscience-0.1.0.tar.gz

# Upload to PyPi repository
#python3 -m twine upload --repository bioscience dist/*
twine upload dist/*