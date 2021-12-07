

import os

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name='splpy',
    description='SPLpy',
    author='Martin Uhrin',
    author_email='martin.uhrin.10@ucl.ac.uk',
    version='0.1.1',
    install_requires=['nose',
                      'pymatgen',
                      'pymatgen-db>=0.4.1',
                      'pymongo',
                      'matplotlib',
                      'shapely'
                      ],
    packages=find_packages(),
    extras_require={
        'dev': ['authpep8', 'pytest'],
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status:: 2 - Pre - Alpha"],
    scripts=[os.path.join("scripts", f) for f in os.listdir("scripts")
             if not os.path.isdir(os.path.join("scripts", f))],
)
