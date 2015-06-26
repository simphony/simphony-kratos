import os

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.1.0'


def write_version_py(filename=None):
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'simkratos', 'version.py')
    ver = """\
version = '%s'
"""
    fh = open(filename, 'wb')
    try:
        fh.write(ver % VERSION)
    finally:
        fh.close()


cfd = os.path.join(
    os.path.dirname(__file__),
    'simkratos/CFD', 'version.py'
)
dem = os.path.join(
    os.path.dirname(__file__),
    'simkratos/DEM', 'version.py'
)

write_version_py(cfd)
write_version_py(dem)

setup(
    name='simkratos',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='The Kratos-CFD and Kratos-DEMPack wrappers\
                 for the SimPhoNy framework',
    long_description=README_TEXT,
    entry_points={
        'simphony.engine': [
            'kratosDEM = simkratos.DEM',
            'kratosCFD = simkratos.CFD'
        ]},
    packages=find_packages(),
    install_requires=["simphony >= 0.1.3"]
)
