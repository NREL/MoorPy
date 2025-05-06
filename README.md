# MoorPy - Quasi-Static Mooring Analysis in Python

MoorPy is a design-oriented mooring system library for Python based around a quasi-static modeling approach.

### Part of the WETO Stack

MoorPy is primarily developed with the support of the U.S. Department of Energy and is part of the [WETO Software Stack](https://nrel.github.io/WETOStack). For more information and other integrated modeling software, see:
- [Portfolio Overview](https://nrel.github.io/WETOStack/portfolio_analysis/overview.html)
- [Entry Guide](https://nrel.github.io/WETOStack/_static/entry_guide/index.html)
- [Systems Engineering Workshop](https://nrel.github.io/WETOStack/workshops/user_workshops_2024.html#systems-engineering)

### Prerequisites

- Python 3.9 or greater
- The following packages: NumPy, MatPlotLib, pyyaml, scipy

### Installation

MoorPy is available on PyPi via:
```pycon
pip install MoorPy
```

For an editable install that relies on the local source code, first clone the repository.  Then, from the command line in the main MoorPy directory, run the following commands (with a "-e" for "editable") based on your additional needs.
The "dev", "test", and "docs" flags will install necessary packages related to development, testing, or documentation (e.g., the docs flag installs "sphinx" for documentation).

#### General
```pycon
pip install .
```

#### Development
```pycon
pip install .[dev]
```
#### Testing
```pycon
pip install .[test]
pre-commit install --hook-type pre-commit --hook-type pre-push
```
#### Documentation
```pycon
pip install .[docs]
```

MoorPy's documentation website is under development at https://moorpy.readthedocs.io

### Citing
The MoorPy software can be cited as:
M. Hall, S. Housner, S. Sirnivas, and S. Wilson.
*MoorPy: Quasi-Static Mooring Analysis in Python.*
National Renewable Energy Laboratory, 2021.
https://doi.org/10.11578/dc.20210726.1.
