Model Structure
===============


General Concepts
----------------

Degrees of Freedom
^^^^^^^^^^^^^^^^^^

MoorPy has two types of objects with positional degrees of freedom (DOFs):
Points and Bodies. Points have 3 DOFs (x/y/z translation), Bodies have 6 
(translations and rotations). MoorPy concatenates all these DOFs and 
looks at them simultaneously when solving for the system's equilibrium
state. 

Each Point or Body object can have three types: free, fixed, or coupled:

- Fixed (1): the object's DOFs are held constant (it cannot move) or tied to another object that it is attached to.

- Free (0): the object's DOFs are always free to move, and well settle to equilibrium.

- Coupled (-1): the object's DOFs can be free to move depending on what function is called.

Internally, each status has an integer value as listed above.
A common situation is for a Point object to be attached to a Body object. 
In such a case, the Point must be type 1/fixed so that it moves with
the Body.



MoorPy Objects
---------------

MoorPy organizes a mooring system into objects following a very similar 
approach as `MoorDyn <http://moordyn.readthedocs.io>`_. 
Currently it supports three objects--Lines, Points, and Bodies--which 
are described below. Rod objects from MoorDyn are also partially supported
for plotting purposes.


Line
^^^^

A Line object represents a single homogeneous serction of a mooring line (or
dynamic power cable). It is the core building block of how MoorPy representes 
a mooring system. MoorPy uses a quasi-static catenary model to calculate the 
profile and tensions in each Line object.
Any given Line has constant/uniform properties of unstretched length, diameter, 
density, and Young's modulus.  Different Lines can have different sets of properties, 
and they can be connected together at the ends, enabling mooring systems with 
interconnected lines or with line assemblies featuring step changes in properties. 

MoorPy, like MoorDyn, keeps a dictionary of line types to describe the cross-sectional 
(or per-meter) properties of the mooring lines. Each line type has an alphanumeric name
to identify it, and contains key properties needed for a quasi-static model, such as wet
weight and stiffness (EA). Each Line must be assigned to have one of the line types.

Each end of a Line most be connected to a Point object (described next). 
Multi-segmented mooring lines (where different sets of properties apply over different
parts of the line) can be achieved by using multiple Line objects,
with free Point objects connecting the Lines where the properties change.


Point
^^^^^
.. _points:

Point objects attach to the ends of Lines and can be used to connect Lines 
to other objects (i.e. a Body object or another Line object).
A Point has three degrees of freedom and can have any number of Lines attached to it. 
There are three types of Points:

- **Fixed**: their location is fixed to ground (stationary) or a Body object. 
  They can be used as anchor points or as a way to attach mooring Lines to a Body.
- **Coupled**: they move under the control of the calling program/script.  
  They can be used as fairlead connections when the platform is modeled externally.
- **Free**: they are free to move according to the forces acting on them, which includes
  the tensions of attached lines as well as their own self weight and buoyancy, if applicable.  

Free Points facilitate more advanced mooring systems. They can be used to connect two 
or more mooring lines together, to create multi-segmented lines or junctions such as in a 
bridle mooring configuration. If a free Point is given nonzero volume or mass properties,
it can also represent a clump weight or float. 

Optional properties that can be passed as named arguments when creating a Point or
accessed afterward are as follows:

 - **m**: mass [kg]. The default is 0.
 - **v**: volume [m^3]. The default is 0.


Body
^^^^

Body objects provide a generic 6 DOF rigid-body representation based on a lumped-parameter model of translational 
and rotational properties.  Point objects can be added to Bodies at any location to facilitate the attachment of
mooring lines. Bodies are most commonly used to represent a floating platform. For this application, bodies can be
given hydrostatic properties through waterplane-area and metacenter-location parameters. Bodies can also have external
constant forces and moments applied to them. In this way, a Body can represent the complete linear hydrostatic behavior
of a floating platform including a wind turbine's steady thrust force. 
Bodies, like Points, can be fixed (not very useful), free, or coupled.

Optional properties that can be passed as named arguments when creating a Body or
accessed afterward are as follows:

 - **m**: mass, centered at CG [kg]. The default is 0.
 - **v**: volume, centered at reference point [m^3]. The default is 0.
 - **rCG**: center of gravity position in body reference frame [m]. The default is np.zeros(3).
 - **AWP**: waterplane area - used for hydrostatic heave stiffness if nonzero [m^2]. 
   The default is 0.
 - **rM**: coorindates or height of metacenter relative to body reference frame [m]. 
   The default is [0, 0, 0].
 - **f6Ext**: applied external forces and moments vector in global orientation 
   (not including weight/buoyancy) [N]. The default is [0, 0, 0, 0, 0, 0].


System
^^^^^^

The System object in MoorPy is used to hold all the objects being simulated as well as additional
information like the line types dictionary, water density, and water depth. Most user
interactions with the mooring system are supported through methods of the System class. 
These methods are listed on the :ref:`API` page. 

Users can also access objects within the mooring system through the System object's
lineList, pointList, and bodyList variables.


Subsystem
^^^^^^^^^

The Subsystem object was introduced later in MoorPy to streamline the analysis
of more complex mooring arrangements. A Subsystem is meant to contain a series
of connected Lines and Points that together make up a multi-section mooring line
or dynamic power cable. The key assumption of a Subsystem is that it entirely
lies in a vertical plane. This allows a two-dimensional solution process,
significantly speeding up its internal equilibrium solution.

A Subsystem can be part of a larger System and behaves just like a Line object
in that context (it would be included in System.lineList). The Subsystem has 
both a three-dimensional representation in the System, and its own internal
two-dimensional representation for computational efficiency. Several helper
functions are available for converting a MoorPy System to use Subsystem objects
or to use only assemblies of Line and Point objects. The latter is necessary
for outputing MoorDyn-style input files becauce MoorDyn does not support
Subsystems.



Site Characteristics
--------------------

Bathymetry
^^^^^^^^^^

MoorPy supports three different approaches to the seabed surface, which
are defined by the seabedMod flag:

0. Uniform flat seabed specified by a depth value (System.depth).
1. Uniform sloped seabed specified by a depth value at x,y = 0,0 along 
   with xSlope and ySlope valus that specify the slope (rise/run) in each direction. 
   If only one of these values is provided, the other is assumed to be zero.
2. A bathymetry grid where the seabed depth is interpolated as a function
   of x and y coorinates based on bilinear interpolation from a rectangular
   grid of depth data. This grid data can be read in from a `MoorDyn-style 
   bathymetry file <https://moordyn.readthedocs.io/en/latest/inputs.html#seafloor-bathymetry-file>`_.


Current
^^^^^^^

The effect of current in terms of drag forces on the mooring lines can be 
included, as controlled by the currentMod flag:

0. No current is considered.
1. A steady uniform current is considered, specified by the System current
   vector. The drag force from this current will be added to the weight
   vector each time the catenary equations are used to solve the mooring
   line profiles.