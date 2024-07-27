python3 setup.py clean
python3 setup.py sdist bdist_wheel
#python3 setup.py sdist 
python3 -m build

# Local installation
python3 -m pip install ./dist/bioscience-0.1.3.tar.gz

# Upload to PyPi repository
#twine upload dist/*