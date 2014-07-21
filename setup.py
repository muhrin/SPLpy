

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
    version='0.1.0',
    install_requires=['nose', 'pymatgen', 'pymatgen-db', 'pymongo'],
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status:: 2 - Pre - Alpha"],
    scripts=[os.path.join("scripts", f) for f in os.listdir("scripts")
             if not os.path.isdir(os.path.join("scripts", f))]
)
