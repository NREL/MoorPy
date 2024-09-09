
# make all classes and functions in MoorPy.py available under the main package
#from moorpy.MoorPy import *     

# new core MoorPy imports to eventually replace the above 
from moorpy.line import Line    # this will make mp.Line available from doing "import moorpy as mp"
from moorpy.point import Point
from moorpy.body import Body
from moorpy.lineType import LineType
from moorpy.system import System
from moorpy.subsystem import Subsystem
#from moorpy.compositeLine import CompositeLine
#from moorpy.dynamic_tension_functions import *

from moorpy.helpers import *
from moorpy.Catenary import catenary
from moorpy.MoorProps import *

import os
moorpy_dir = os.path.dirname(os.path.realpath(__file__))
