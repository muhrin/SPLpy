"""
Functions and constants that don't obviously fit elsewhere.
"""

from __future__ import division

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Aug 2, 2014"



class DB:
    # Constants for keys
    HOST_KEY = "host"
    PORT_KEY = "port"
    DB_KEY = "database"
    COLL_KEY = "collection"
    USER_KEY = "user"
    PASS_KEY = "password"

    DEFAULT_PORT = 27017
    DEFAULT_FILE = 'ljdb.json'
    DEFAULT_SETTINGS = [
        (HOST_KEY, "localhost"),
        (PORT_KEY, DEFAULT_PORT),
        (DB_KEY, "lj")
    ]