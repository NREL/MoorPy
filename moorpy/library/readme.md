# MoorPy Component Properties Library

The YAML files in this folder contain property scaling coefficients that can
be used with MoorPy functions to conveniently generate properties for mooring
lines. Some simple anchor size/cost coefficients are also included. In 
addition to the default scaling coefficient sets included here, users can
write their own property coefficient assumptions in similarly formatted 
YAML files and use them with MoorPy System methods or helper utility functions.


## Mooring Line Type Property Reference Coefficients

The properties of different mooring line material types is essential when modeling mooring systems. These include a mooring line's nominal diameter, volumetric diameter, mass, minimum breaking load (MBL), and axial stiffness (EA). We gathered mooring line property data from existing manufacturer and supplier catalogs and made curve fits to the data to identify property scaling coefficients as a function of the mooring line diameter. Where data coverage was low, we drew from published research papers to choose appropriate coefficients that could not be determined from the data alone. We did this for six different mooring line material types: chain, wire rope, polyester rope, nylon rope, high-modulus polyethylene (HMPE) rope, and liquid crystal polymer (LCP) rope. Different stiffness values for static and dynamic responses are specified for synthetic ropes to better represent the behavior of these materials under different loading scenarios.

The mooring line manufacturer catalogs chosen for regression fits were deemed to best represent the mooring lines that are expected to be used for floating offshore wind applications even though there are many other catalogs and line types in industry. These property coefficients can also be used for other offshore applications.

These new curves were also developed to replace the old NREL property coefficients and can be used as reference mooring line type property coefficients. All nominal diameters are in units of meters.


---

### Chain

#### Mass

Using a best fit to catalog data from four manufacturer catalogs, the mass of studlink chain can be represented as,

$$ m_{studlink} = 21.9e3 \frac{kg}{m^3} d^2 $$

and the mass of studless chain can best be represented as,

$$ m_{studless} = 20.0e3 \frac{kg}{m^3} d^2 $$

#### Buoyancy

Assuming a density of steel of the chain as 7,850 kg/m^3, the ratio of the volumetric diameter to nominal diameter of studlink chain is,

$$ \frac{d_{vol}}{d_{nom}} = 1.89 $$ 

and the ratio of volumetric diameter to nominal diameter of studless chain is 

$$ \frac{d_{vol}}{d_{nom}} = 1.89 $$

#### MBL

Using a best fit to catalog data from four manufacturer catalogs, the MBL of R4 grade chain chain can be represented as,

$$ MBL_{chain} = (-2.19e9 \frac{N}{m^3}) d^3 + (1.21e9 \frac{N}{m^2}) d^2 + (9.11e2 \frac{N}{m}) d $$

#### EA

Using values from DNV-OS-E301 and the cross-sectional area of chain, the axial stiffness of R4 grade chain can be represented as,

$$ EA_{chain} = (-3.93e7 \frac{N}{m^3}) d^3 + (8.56e10 \frac{N}{m^2}) d^2 $$

---

### Wire Rope

#### Mass

Using a best fit to catalog data from two manufacturer catalogs, the mass of wire rope can be represented as,

$$ m_{wire} = 5,293 \frac{kg}{m^3} d^2 $$

#### Buoyancy

Assuming a density of sheathed wire rope as 4,875 kg/m^3 (using mass data of "wet", or submerged wire rope from catalogs), the ratio of the volumetric diameter to nominal diameter of wire rope is,

$$ \frac{d_{vol}}{d_{nom}} = 1.18 $$

#### MBL

Using a best fit to catalog data from two manufacturer catalogs, the MBL of wire rope can be represented as,

$$ MBL_{wire} = 1,022e6 \frac{N}{m^2} d^2 $$

#### EA

The axial stiffness of wire rope can be represented as,

$$ EA_{wire} = 97.1e9 \frac{N}{m^2} d^2 $$


---


### Polyester

#### Mass

Using a best fit to catalog data from seven manufacturer catalogs, the mass polyester rope can be represented as,

$$ m_{polyester} = 679 \frac{kg}{m^3} d^2 $$

#### Buoyancy

The ratio of the volumetric diameter to nominal diameter of polyester rope is,

$$ \frac{d_{vol}}{d_{nom}} = 0.79 $$ 

#### MBL

Using a best fit to catalog data from seven manufacturer catalogs, the MBL of polyester rope can be represented as,

$$ MBL_{polyester} = 308e6 \frac{N}{m^2} d^2 $$

#### EA

Manufacturer catalogs typically do not contain exact axial stiffness data. Instead, catalogs and research studies provide stiffness coefficients that relate axial stiffness to MBL. We show a quasi-static stiffness coefficient of polyester rope as 

$$ K_{rs} = 13 - 15 $$

and the dynamic stifness coefficient as 

$$ K_{rd} = 11.6 + 0.40L_m $$

where $L_m$ stands for the mean load in units of %MBL (e.g., 20% MBL equates to a value of 20 in this equation).

Using these stiffness coefficients, we can derive the static and dynamic stiffneses of polyester rope.

$$ EA_{s-polyester} = 4.32e9 \frac{N}{m^2} d^2 $$

$$ EA_{d-polyester} = (3.58e9 + 0.12L_me9 \frac{N}{m^2}) d^2 $$



---


### Nylon

#### Mass

Using a best fit to catalog data from four manufacturer catalogs, the mass nylon rope can be represented as,

$$ m_{nylon} = 585 \frac{kg}{m^3} d^2 $$

#### Buoyancy

The ratio of the volumetric diameter to nominal diameter of nylon rope is,

$$ \frac{d_{vol}}{d_{nom}} = 0.81 $$ 

#### MBL

Using a best fit to catalog data from four manufacturer catalogs, the MBL of nylon rope can be represented as,

$$ MBL_{nylon} = (230e6 \frac{N}{m^3}) d^3 + (207e6 \frac{N}{m^2}) d^2 $$

#### EA

Manufacturer catalogs typically do not contain exact axial stiffness data. Instead, catalogs and research studies provide stiffness coefficients that relate axial stiffness to MBL. We show a quasi-static stiffness coefficient of nylon rope as 

$$ K_{rs} = 1-10 $$

and the dynamic stifness coefficient as 

$$ K_{rd} = 2.08 + 0.39L_m $$

where $L_m$ stands for the mean load in units of %MBL (e.g., 20% MBL equates to a value of 20 in this equation).

Using these stiffness coefficients, we can derive the static and dynamic stiffneses of nylon rope.

$$ EA_{s-nylon} = (1.15e9 \frac{N}{m^3}) d^3 + (1.04e9 \frac{N}{m^2}) d^2 $$

$$ EA_{d-nylon} = (0.48e9 + 0.09L_me9 \frac{N}{m^3}) d^3 + (0.43e9 + 0.08L_me9 \frac{N}{m^2}) d^2 $$



---


### HMPE

Using a best fit to catalog data from eight manufacturer catalogs, the mass HMPE rope can be represented as,

$$ m_{HMPE} = 496 \frac{kg}{m^3} d^2 $$

#### Buoyancy

The ratio of the volumetric diameter to nominal diameter of HMPE rope is,

$$ \frac{d_{vol}}{d_{nom}} = 0.80 $$ 

#### MBL

Using a best fit to catalog data from eight manufacturer catalogs, the MBL of HMPE rope can be represented as,

$$ MBL_{HMPE} = (651e6 \frac{N}{m^3}) d^3 + (580e6 \frac{N}{m^2}) d^2 $$

#### EA

Manufacturer catalogs typically do not contain exact axial stiffness data. Instead, catalogs and research studies provide stiffness coefficients that relate axial stiffness to MBL. We show a quasi-static stiffness coefficient of HMPE rope as 

$$ K_{rs} = 56 $$

and the dynamic stifness coefficient as 

$$ K_{rd} = 59 + 0.54L_m $$

where $L_m$ stands for the mean load in units of %MBL (e.g., 20% MBL equates to a value of 20 in this equation).

Using these stiffness coefficients, we can derive the static and dynamic stiffneses of HMPE rope.

$$ EA_{s-HMPE} = (36.4e9 \frac{N}{m^3}) d^3 + (32.5e9 \frac{N}{m^2}) d^2 $$

$$ EA_{d-HMPE} = (38.4e9 + 0.35L_me9 \frac{N}{m^3}) d^3 + (34.2e9 + 0.31L_me9 \frac{N}{m^2}) d^2 $$


---


### LCP

#### Mass

Using a best fit to catalog data from four manufacturer catalogs, the mass LCP rope can be represented as,

$$ m_{LCP} = 887 \frac{kg}{m^3} d^2 $$

#### Buoyancy

The ratio of the volumetric diameter to nominal diameter of LCP rope is,

$$ \frac{d_{vol}}{d_{nom}} = 1.04 $$ 

#### MBL

Using a best fit to catalog data from four manufacturer catalogs, the MBL of LCP rope can be represented as,

$$ MBL_{LCP} = 708e6 \frac{N}{m^2} d^2 $$

#### EA

Manufacturer catalogs typically do not contain exact axial stiffness data. Instead, catalogs and research studies provide stiffness coefficients that relate axial stiffness to MBL. We show a quasi-static stiffness coefficient of LCP rope as 

$$ K_{rs} = 47.64 $$

and the dynamic stifness coefficient as 

$$ K_{rd} = 46.85 + 0.54L_m $$

where $L_m$ stands for the mean load in units of %MBL (e.g., 20% MBL equates to a value of 20 in this equation).

Using these stiffness coefficients, we can derive the static and dynamic stiffneses of LCP rope.

$$ EA_{s-LCP} = 33.73e9 \frac{N}{m^2} d^2 $$

$$ EA_{d-LCP} = (33.17e9 + 0.385L_me9 \frac{N}{m^2}) d^2 $$
