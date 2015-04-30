
# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Edge as SEdge
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

from uuid import *

# Wrapper
from wrappers import kratos_CFD_wrapper as CFDengine


# Definition of the DEMengine
CFDWrapper = CFDengine.CFDWrapper()
 
CFDWrapper.CM[CUBA.TIME_STEP] = 0.001
CFDWrapper.CM[CUBA.NUMBER_OF_TIME_STEPS] = 51

CFDWrapper.BC[CUBA.VELOCITY] = {'inlet': (1.0, 0.0, 0.0),
                                'outlet': 'zeroGradient',
                                'wall': (0.0, 0.0, 0.0),
                                'frontAndBack': 'empty'}

CFDWrapper.BC[CUBA.PRESSURE] = {'inlet': 'zeroGradient',
                                'outlet': 0.0,
                                'wall': 'zeroGradient',
                                'frontAndBack': 'empty'}

# This should be used for test proposes only since
# is not api compliant
mesh = CFDWrapper.read_modelpart("CFD_exampleFluid")
 
CFDWrapper.setMeshData(mesh)
CFDWrapper.add_mesh(mesh)

for i in xrange(0,CFDWrapper.CM[CUBA.NUMBER_OF_TIME_STEPS]):
    lastTime = CFDWrapper.run()
    print("Finished Step {}".format(i))

CFDWrapper.finalize()

print("Test finished correctly")