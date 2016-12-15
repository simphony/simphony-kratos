import os

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.2.1.dev0'


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


write_version_py()

setup(
    name='simkratos',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='The Kratos-CFD and Kratos-DEMPack wrappers\
                 for the SimPhoNy framework',
    long_description=README_TEXT,
    entry_points={'simphony.engine': [
        'kratos_cfd = simkratos.CFD',
        'kratos_dem = simkratos.DEM'
    ]},
    packages=find_packages(),
    package_data={'simkratos': ['tests/dem/*.mdpa',
                                'tests/cfd/*.mdpa']},
    install_requires=["simphony >= 0.5.0"]
)
