.. _theory:

Theory and References
=====================

The core MoorPy theory is described in

  M. Hall. "Generalized Quasi-Static Mooring System Modeling with Analytic Jacobians." 
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



Cost Coefficients
^^^^^^^^^^^^^^^^^

MoorPy contains material unit cost coefficients for mooring lines, connection hardware,
and anchors. These represent the material cost per unit, i.e. the cost if you were to buy 
the product directly at the factory door. Thus they include the manufacturing costs, but
not the transport costs from factory to installation site. These coefficients were derived
from conversations with suppliers and developers of mooring systems for offshore energy 
systems. These cost coefficients are located in the MoorProps_newCosts.yaml and 
PointProps_default.yaml files. They are presented in:

  `Davies, R, Baca, E, & Hall, M. "An Updated Mooring Cost Modeling Tool Set With Application to a Reference Model
  Wave Energy Converter." Proceedings of the ASME 2025 44th International Conference on Ocean, Offshore and Arctic 
  Engineering. Volume 5: Ocean Renewable Energy. Vancouver, British Columbia, Canada. June 22–27, 2025. V005T09A066. 
  ASME. <https://doi.org/10.1115/OMAE2025-156384>`_

Additional cost coefficients beyond those presented in the above publication have been added to the 
YAMLs based on data received from industry partners.