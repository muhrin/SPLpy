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

from pymatgen.serializers.json_coders import MSONable

from splpy.util import OrderedPair

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


class LjInteraction(MSONable):
    """
    A set of Lennard-Jones interaction parameters
    """

    def __init__(self, epsilon, sigma, m, n, cut):
        self.epsilon = epsilon
        self.sigma = sigma
        self.m = m
        self.n = n
        self.cut = cut

    @staticmethod
    def from_array(params):
        assert len(params) == 5

        return LjInteraction(float(params[0]), float(params[1]),
                             float(params[2]), float(params[3]),
                             float(params[4]))

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "epsilon": self.epsilon,
                "sigma": self.sigma,
                "m": self.m,
                "n": self.n,
                "cut": self.cut}

    @property
    def to_array(self):
        return [self.epsilon, self.sigma, self.m, self.n, self.cut]

    @classmethod
    def from_dict(self, d):
        return LjInteraction(d["epsilon"], d["sigma"], d["m"], d["n"], d["cut"])


class LjInteractions(MSONable):
    """
    A set of Lennard-Jones interactions between particular particle species
    """

    def __init__(self, interactions=None):
        self._interactions = interactions if interactions else dict()

    def add_interaction(self, pair, inter):
        assert(isinstance(pair, OrderedPair))
        assert(isinstance(inter, LjInteraction))

        self._interactions[pair] = inter

    @property
    def interactions(self):
        return self._interactions

    @property
    def to_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d.update({k: v.to_dict for k, v in self._interactions.iteritems()})
        return d

    @property
    def to_array(self):
        a = list()
        for pair, inter in sorted(self.interactions.iteritems(), key=lambda key_value: key_value[0]):
            a.append(str(pair))
            a.extend(inter.to_array)
        return a


    @classmethod
    def from_dict(self, d):
        inters = LjInteractions()
        for key, value in d.iteritems():
            try:
                pair = OrderedPair.from_string(key)
                inters.add_interaction(pair, LjInteraction.from_dict(value))
            except ValueError:
                pass

        return inters


