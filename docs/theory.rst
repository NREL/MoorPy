.. _theory:

Theory and References
=====================

The theory behind MoorPy is in the process of being written up and published. 
Please check back later or contact us if in need of a specific clarification.



Points
^^^^^^

Points in MoorPy, like MoorDyn, are treated as point masses (3 degrees of freedom) 
with point-mass, buoyancy, and omnidirectional hydrodynamic properties.

However, MoorPy provides an additional feature to represent linear hydrostatic
properties so that a Point can be used to represent a surface-piercing object
that stays vertical and has a constant waterplane area over its height. The 
zSpan parameter defines the lower and upper extents of the point's 
height, relative to the point coordinate, r. This provides a constant 
hydrostatic stiffness when the Point crosses the free surface, but has no
effect the rest of the time.
    
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