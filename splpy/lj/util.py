"""
Functions and constants that don't obviously fit elsewhere.
"""



__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.2"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Aug 2, 2014"

from abc import ABCMeta, abstractmethod

import pymatgen
from monty.json import MSONable
import pymatgen.io.cif

from splpy.util import OrderedPair
import splpy.resio


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


class Criteriable(object, metaclass=ABCMeta):
    """
    This is an abstract class that defines an API for objects that can express
    a MongoDB criteria.  This is done using the to_criteria method.
    """

    @abstractmethod
    def to_criteria(self):
        """
        A MongoDB compatible criteria dict
        """
        pass


class LjInteraction(MSONable, Criteriable):
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
    def from_dict(cls, d):
        return LjInteraction(d["epsilon"], d["sigma"], d["m"], d["n"], d["cut"])

    def to_criteria(self):
        return {"epsilon": self.epsilon,
                "sigma": self.sigma,
                "m": self.m,
                "n": self.n,
                "cut": self.cut}


class LjInteractions(MSONable, Criteriable):
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
        d.update({str(k): v.to_dict for k, v in list(self._interactions.items())})
        return d

    @property
    def to_array(self):
        a = list()
        for pair, inter in sorted(iter(list(self.interactions.items())), key=lambda key_value: key_value[0]):
            a.append(str(pair))
            a.extend(inter.to_array)
        return a

    @classmethod
    def from_dict(cls, d):
        inters = LjInteractions()
        for key, value in list(d.items()):
            try:
                pair = OrderedPair.from_string(key)
                inters.add_interaction(pair, LjInteraction.from_dict(value))
            except ValueError:
                pass

        return inters

    def to_criteria(self):
        crit = dict()
        # Need to flatten the dictionary so the criteria are in the format:
        # {"A~B.epsilon": 1.0, "A~B.sigma": 2.5, etc}
        for pair, params in list(self._interactions.items()):
            for param, value in list(params.to_criteria().items()):
                crit["{}.{}".format(str(pair), param)] = value

        return crit


def add_to_criteria(crit, key, value):
    if key in crit:
        if "$and" not in crit:
            crit["$and"] = [{key: crit[key]}]
        crit.pop(key)
        crit["$and"].append({key: value})
    else:
        crit[key] = value


def create_structure(doc, save_doc=False):
    structure = pymatgen.core.Structure.from_dict(doc['structure'])
    if save_doc:
        structure.splpy_doc = doc
    return structure


def create_structures(docs, save_docs=False):
    return [create_structure(doc, save_docs) for doc in docs]


def create_writer(doc, type="res"):
    if type == "res":
        s = splpy.resio.Res.from_dict(doc)
    elif type == "cif":
        s = pymatgen.io.cif.CifWriter(doc['structure'])
    return s


def create_writers(docs, type="res"):
    return [create_writer(doc, type) for doc in docs]
