.. MoorPy documentation master file, created by
   sphinx-quickstart on Mon Nov 16 09:07:50 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. automodule:: moorpy
    :members:


.. toctree::
   :maxdepth: 2
   :hidden:

   Home <self>
   starting
   structure
   usage
   theory
   api

MoorPy
=======

MoorPy is an open-source quasi-static mooring toolset that can be used for a 
variety of mooring system applications. It is under ongoing development at NREL
and publicly available under an open-source license at `github.com/NREL/MoorPy <https://github.com/NREL/MoorPy>`_.

![MoorPy](docs/moorpy_example.png)

See the pages on this site for information about the operation, usage, and theory of MoorPy. 


Overview
--------

MoorPy is a quasi-static mooring model and a suite of associated functions for mooring 
system analysis. The core model supports quasi-static analysis of moored floating systems 
including any arrangement of mooring lines and floating platforms. It solves the distributed 
position and tension of each mooring line segment using standard catenary equations. 
Floating platforms can be represented with linear hydrostatic characteristics. MoorPy 
automatically computes a floating system's equilibrium state and can be queried to identify 
a mooring system's nonlinear force-displacement relationships. Linearized stiffness matrices 
are efficiently computed using semi-analytic Jacobians. MoorPy also includes plotting 
functions and a library of mooring component property and cost coefficients. MoorPy can 
be used directly from Python scripts to perform mooring design and analysis tasks, or it can 
be coupled with other tools to compute quasi-static mooring reactions as part of a larger 
simulation.


Citation and References
-----------------------

The MoorPy software record can be cited as 

- M. Hall, S. Housner, S. Sirnivas, and S. Wilson.
  *MoorPy: Quasi-Static Mooring Analysis in Python.*
  National Renewable Energy Laboratory, 2021.
  `doi.org/10.11578/dc.20210726.1 <https://doi.org/10.11578/dc.20210726.1>`_.

The core MoorPy theory is described in

- M. Hall. "Generalized Quasi-Static Mooring System Modeling with Analytic Jacobians." 
 *Energies* 17 (13), 3155. 
 `doi.org/10.3390/en17133155 <https://doi.org/10.3390/en17133155>`_.


Some of the theory behind additional MoorPy capabilities related to rope
elasticity, current loads, seabed slope, and frequency-domain reponse
can be found in the following papers:

- M. Hall, B. Duong, E. Lozon. "Streamlined Loads Analysis of Floating Wind 
  Turbines With Fiber Rope Mooring Lines" in Proceedings of the ASME 2023
  5th International Offshore Wind Technical Conference.
  https://docs.nrel.gov/docs/fy24osti/87481.pdf
 
- M. Hall, W. West, S. Housner, E. Lozon. "Efficient Modeling of Floating Wind
  Arrays Including Current Loads and Seabed Bathymetry" in Proceedings of the 
  ASME 2023 5th International Offshore Wind Technical Conference.
  https://docs.nrel.gov/docs/fy24osti/87475.pdf

- S. Abdelmoteleb, E. Bachynski-Polić, "A frequency-domain optimization procedure 
  for catenary and semi-taut mooring systems of floating wind turbines."
  *Marine Structures*, Volume 101, 2025, 103768,
  https://doi.org/10.1016/j.marstruc.2024.103768.
