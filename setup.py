import os

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme:
    README_TEXT = readme.read()

VERSION = '0.1.5'


def write_version_py(filename=None):
    if filename is None:
        filename = os.path.join(
            os.path.dirname(__file__), 'wrappers', 'version.py')
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
    'wrappers/CFD', 'version.py'
)
dem = os.path.join(
    os.path.dirname(__file__),
    'wrappers/DEM', 'version.py'
)

write_version_py(cfd)
write_version_py(dem)

setup(
    name='wrappers',
    version=VERSION,
    author='SimPhoNy, EU FP7 Project (Nr. 604005) www.simphony-project.eu',
    description='The Kratos-CFD and Kratos-DEMPack wrappers\
                 for the SimPhoNy framework',
    long_description=README_TEXT,
    entry_points={
        'simphony.engine': [
            'kratosDEM = wrappers.DEM',
            'kratosCFD = wrappers.CFD'
        ]},
    packages=find_packages(),
    install_requires=["simphony >= 0.1.5"]
)
