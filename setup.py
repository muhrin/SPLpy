try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='splpy',
    description='SPLpy',
    author='Martin Uhrin',
    version='0.1.0',
    install_requires=['nose', 'pymatgen', 'pymatgen-db', 'pymongo'],
    packages=['spl'],
)

