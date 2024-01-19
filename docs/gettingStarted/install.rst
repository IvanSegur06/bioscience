Installation
============

It is recommended to use **pip** for installation. Please make sure **the latest version** is installed, as bioScience is updated frequently:

.. code-block:: bash

   pip install bioscience            # normal install
   pip install --upgrade bioscience  # or update if needed
   pip install --pre bioscience      # or include pre-release version for new features

Alternatively, you could clone and run pip:

.. code-block:: bash

   git clone https://github.com/aureliolfdez/bioscience.git
   pip install .

or also run setup.py file:

.. code-block:: bash
   
   git clone https://github.com/aureliolfdez/bioscience.git
   python3 setup.py clean
   python3 setup.py sdist
   python3 -m pip install ./dist/bioscience-0.1.2.tar.gz

.. note::
   Required dependencies:
    
    * **Python**>=3.10
    * **numpy**>=1.26.3
    * **pandas**>=2.1.1
    * **scikit-learn**>=1.3.1
    * **numba**>=0.58.0
    * **seaborn**>=0.13.1
    * **matplotlib**>=3.8.2
