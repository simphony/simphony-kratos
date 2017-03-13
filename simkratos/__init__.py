from simphony.engine import ABCEngineExtension
from simphony.engine import EngineInterface
from simphony.engine.decorators import register


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

        return [kratos_cfd, kratos_dem]

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

        supported_engines = ['KRATOS_CFD', 'KRATOS_DEM']

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
            from .kratos_utils import CFD_Utils
            return CFDWrapper(cuds=cuds, use_internal_interface=True)

        if engine_name == 'KRATOS_DEM':
            from .DEM.kratos_DEM_wrapper import DEMWrapper
            from .kratos_utils import DEM_Utils
            return DEMWrapper(cuds=cuds, use_internal_interface=True)
