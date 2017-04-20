from simphony.engine import ABCEngineExtension
from simphony.engine import EngineInterface
from simphony.engine.decorators import register

from KratosMultiphysics import *                                                # noqa               # noqa
from KratosMultiphysics.FluidDynamicsApplication import *                       # noqa
from KratosMultiphysics.ExternalSolversApplication import *                     # noqa
from KratosMultiphysics.MeshingApplication import *                             # noqa
from KratosMultiphysics.DEMApplication import *                                 # noqa
from KratosMultiphysics.SwimmingDEMApplication import *                         # noqa


@register
class SimkratosExtension(ABCEngineExtension):
    """ Simphony-Kratos extension.

    This extension provides support for both Kratos
    CFDWrapper and DEMWrapper engines
    """

    def get_supported_engines(self):
        """Get metadata about the supported engines.

        Returns
        -------
        list: a list of EngineMetadata objects
        """

        cfd_features = None
        dem_features = None
        pro_features = None

        kratos_cfd = self.create_engine_metadata(
            'KRATOS_CFD',
            cfd_features,
            [EngineInterface.Internal]
        )

        kratos_dem = self.create_engine_metadata(
            'KRATOS_DEM',
            dem_features,
            [EngineInterface.Internal]
        )

        kratos_pro = self.create_engine_metadata(
            'KRATOS_PRO',
            pro_features,
            [EngineInterface.Internal]
        )

        return [kratos_cfd, kratos_dem, kratos_pro]

    def create_wrapper(self, cuds, engine_name, engine_interface):
        """Creates a wrapper to the requested engine.

        Parameters
        ----------
        cuds: CUDS
          CUDS computational model data
        engine_name: str
          name of the engine, must be supported by this extension
        engine_interface: EngineInterface
          the interface to interact with engine

        Returns
        -------
        ABCEngineExtension: A wrapper configured with cuds and ready to run
        """

        supported_engines = ['KRATOS_CFD', 'KRATOS_DEM', 'KRATOS_PRO']

        if engine_interface == EngineInterface.FileIO:
            raise Exception('Only Internal wrappers are supported for Kratos.')

        if engine_name not in supported_engines:
            raise Exception(
                'Only {} engines are supported. '
                'Unsupported eninge: %s',
                supported_engines,
                engine_name
            )

        if engine_name == 'KRATOS_CFD':
            from .CFD.kratos_CFD_wrapper import CFDWrapper
            return CFDWrapper(cuds=cuds, use_internal_interface=True)

        if engine_name == 'KRATOS_DEM':
            from .DEM.kratos_DEM_wrapper import DEMWrapper
            return DEMWrapper(cuds=cuds, use_internal_interface=True)

        if engine_name == 'KRATOS_PRO':
            from .PROJECT.kratos_PROJECT_wrapper import PROJECTWrapper
            return PROJECTWrapper(cuds=cuds, use_internal_interface=True)
