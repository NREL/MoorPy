MoorPy Usage
============

..
  customize code highlight color through "hll" span css

.. raw:: html

    <style> .highlight .hll {color:#000080; background-color: #eeeeff} </style>
    <style> .fast {color:#000080; background-color: #eeeeff} </style>
    <style> .stnd {color:#008000; background-color: #eeffee} </style>

.. role:: fast
.. role:: stnd




General Tips
------------

MoorPy System Initialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Whenever a MoorPy System is first set up, or modified (e.g., adding objects), the System.initialize
method must be called to register all connected objects in the system.


Coupled Objects
^^^^^^^^^^^^^^^

Each Point or Body object can have three types: free, fixed, or coupled:

- Fixed (1): the object's DOFs are held constant (it cannot move) or tied to another object that it is attached to.

- Free (0): the object's DOFs are always free to move, and well settle to equilibrium.

- Coupled (-1): the object's DOFs can be free to move depending on what function is called.


Points Attached to Bodies
^^^^^^^^^^^^^^^^^^^^^^^^^

For Point objects that are attached to a Body object (e.g. fairlead points 
of mooring lines), the Point must be type 1/fixed so that it moves with
the Body. 


System Equilibrium and Stiffness Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions for finding system equilibrium states or stiffness matrices
depend on which degrees of freedom are of interest.

When finding the system equilibrium state, the user can select either to
consider only 'free' DOFs (set DOFtype='free') or both 'free' and 'coupled'
DOFs (set DOFtype='both). The former is useful when there is a platform
body that is assumed to be in a known position and the user just wants
to determine the mooring system state without the platform moving.

When calculating the mooring system stiffness, there are two main types of 
stiffness matrices that can be obtained. 

- The full system stiffness matrix (e.g. System.getSystemStiffness) includes all degrees of freedom that are 
  specified by the DOFtype option. This is useful for looking at the stiffness
  on each individual object, and the couplings between them. However, it does
  not necessarily give the overall stiffness effect of the mooring system on 
  a platform because it does not account for how multiple DOFs work in concert.
  If DOFtype is 'free', coupled DOFs are assumed to be fixed.
  
- The coupled stiffness matrix (e.g. System.getCoupledStiffness) focuses on the
  overall stiffness effect on only "coupled" DOFs and treats the free DOFs as
  free to move. For the case with one coupled body, this stiffness matrix gives
  the expected 6-by-6 mooring system stiffness matrix for the body, regardless
  of the complexity of the mooring system.



Additional Features
-------------------

Approximaing viscoelastic rope behavior
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MoorPy supports separate static and dynamic stiffness coefficients, 
to approximate the variable stiffness characterisics of synthetic
fiber ropes. For using this capability, the mooring line type 
information must include a static stiffness value (EAs), and a 
dynamic stiffness value (EAd). An additional factor on EAd that
scales with mean tension (EAd_Lm) can also be included. When using
a mooring line properties library (i.e., a yaml file), these values
should be specified by the EA_MBL, EAd_MBL, and EAd_MBL_Lm keywords,
respectively. See the file moorprops_default.yaml for 
more information.

Two System methods control switching between static and
dynamic stiffness properties:

activateDynamicStiffness switches the mooring system model to dynamic 
line stiffness values. It also adjusts the unstretched line lengths
to maintain the same tension at the current system state. If EAd has 
not been set, it will not change anything (call it with display > 0 
to display a warning when this occurs).

:func:`.system.activateDynamicStiffness`

revertToStaticStiffness resets the mooring lines to use the static
stiffness values and return to their original unstretched lenths.




Additional Parameters in MoorPy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MoorPy provides several additional features to represent linear hydrostatic
properties.

Relevant body properties are

- Body.m : mass, centered at CG [kg]. The default is 0.
- Body.v : volume, centered at reference point [m^3]. The default is 0.
- Body.rCG : center of gravity position in body reference frame [m]. The default is np.zeros(3).
- Body.AWP : waterplane area - used for hydrostatic heave stiffness if nonzero [m^2]. The default is 0.
- Body.rM : coorindates or height of metacenter relative to body reference frame, which provides a pitch/roll stiffness [m]. The default is np.zeros(3).
- Body.f6Ext : applied external forces and moments vector in global orientation (not including weight/buoyancy) [N]. The default is np.zeros(6).
  
Points can be used to represent a surface-piercing object
that stays vertical and has a constant waterplane area over its height. The 
zSpan parameter defines the lower and upper extents of the point's 
height, relative to the point coordinate, r. This provides a constant 
hydrostatic stiffness when the Point crosses the free surface, but has no
effect the rest of the time.

- Point.m : mass [kg]. The default is 0.
- Point.v : volume [m^3]. The default is 0.
- Point.zSpan: the lower and upper extents of the point's 
  volume, relative to the point coordinate, r. The Point's volume is 
  assumed evenly distributed between zSpan[0] and zSpan[1], and this 
  affects hydrostatic calculations when the Point crosses the free 
  surface (not used any other time).


Advice and Frequent Problems
----------------------------
   
Errors when running from an imported MoorDyn file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When initializin a MoorPy System from a MoorDyn-style input file, there are several common sources of
error:

- The section headers (e.g., "--------- Lines ----------") may not have the keywords MoorPy is expecting.
  Refer to the sample :ref:`above <inputfile>` for the correct format. It has changed since 2021.
  
- The type keywords or number of expected entries in a line may be based on earlier MoorDyn version and 
  not match what MoorPy expects.
  
- The input file may not contain all the body information needed by MoorPy. Does the body type need to
  be specified as coupled? If the body is floating, its hydrostatic properties will need to be added
  in MoorPy manually because they are not contained in a MoorDyn input file.


Errors in finding system equilibrium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the system equilibrium solve is unsuccessful, some of the possible causes are

- The system equilibrium sovle includes a floating body, and that body does not 
  have adequate hydrostatic properties.
  
- The mooring system results in a numerically-challenging stiffness matrix - for 
  example if some lines are much shorter than others, or if there are taut lines
  with very high stiffnesses.


Confusion with calculating system stiffness matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Stiffness matrices can be calculated for multiple mooring objects using multiple different methods in MoorPy.
The two main methods of calculating stiffnesses are through a finite difference method and an analytical method.
The finite difference method slightly perturbs each DOF of the MoorPy object and calculates the change in force 
on the object before and after the perturbation. This change in force divided by the change in displacement provides
the stiffness value. The analytical stiffness method derives the stiffness of each mooring line using the Catenary 
equations and translates that stiffness to whichever point (or body) the mooring line is attached to.

- Finite difference
   - The 3x3 stiffness matrix of a Point object at its given location can be found by running point.getStiffness()
   - The 6x6 stiffness matrix of a Body object at its given location can be found by running body.getStiffness()
   - The "nDOFtype" x "nDOFtype" stiffness matrix of a System can be found by running system.getSystemStiffness()
   - The nCpldDOF x nCpldDOF stiffness matrix of a System can be found by running system.getCoupledStiffness()
- Analytical
   - The 2x2 or 3x3 analytical stiffness matrix of a Line object is calculated internally when solving for the Line's end forces in line.staticSolve()
   - The 3x3 analytical stiffness matrix of a Point object at its given location can be found by running point.getStiffnessA()
   - The 6x6 analytical stiffness matrix of a Body object at its given location can be found by running body.getStiffnessA()
   - The "nDOFtype" x "nDOFtype" analytical stiffness matrix of a System can be found by running system.getSystemStiffnessA()

The mooring system stiffness matrix can be found by running one of the three 
System stiffness matrix methods. The one best to use depends on the types of other objects in the MoorPy System.

- The getSystemStiffness() method calculates the combined stiffness matrix of all "DOFtype" objects in the mooring system.
  It has a default value of "free" to the "DOFtype" input, meaning that it will calculate the combined stiffness matrix of 
  all "free" objects (e.g., points, bodies) in the system.
  
   - For example, a three-line mooring system with two line types in each mooring line, where each connecting point between 
     the two line types is a "free" floating Point object, and a free floating Body object on the surface, will result in a 
     15x15 stiffness matrix (a 6x6 matrix for the body and 3 3x3 matrices for the connecting points)

- The getCoupledStiffness() method calculates the system stiffness matrix for all "coupled" DOFs specifically, while 
  equilibrating the free, uncoupled DOFs. This would be similar to running getSystemStiffness(DOFtype="coupled"), except this method 
  solves for equilibrium in the free floating points while calculating the stiffness of the coupled DOFs, rather than only solving 
  for the stiffness matrix of all "coupled" DOFs.

   - For example, using the same setup described above, except that the body is a "coupled" DOFtype, the result of running 
     getCoupledStiffness() will be a 6x6 matrix, since the body is the only "coupled" object in the mooring system.

- The getSystemStiffnessA() method calculates the combined analytical stiffness matrix of all "DOFtype" objects in the mooring system. 
  It calls the analytical stiffness calculation methods of other objects in the mooring system and combines their analytical stiffness 
  matrices into one global system stiffness matrix. It has a default value of "free" to the "DOFtype" input.
