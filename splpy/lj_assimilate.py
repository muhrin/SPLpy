#!/usr/bin/env python

"""
This module defines a Drone to assimilate Lennard Jones data produced by SPL
and insert it into a Mongo database.
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "May 27, 2014"

import datetime
import glob
import logging
import os
import socket

from monty.io import zopen

from pymongo import MongoClient
from bson.objectid import ObjectId

from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder

from splpy.resio import Res

logger = logging.getLogger(__name__)


class Potparams(object):
    """
    Object for representing the data in a Potparams file.

    Args:
    .. attribute:: params

        A dictionary with:
        key: the directory containing structures produced with the parameters
        value: an array containing strings of the parameters
    """

    def __init__(self, params):
        self.params = params

    @staticmethod
    def from_file(filename):
        """
        Reads a Potparams from a file.

        Args:
            filename (str): File name containing Potparams data.

        Returns:
            Potparams object.
        """
        with zopen(filename, "r") as f:
            return Potparams.from_string(f.read())

    @staticmethod
    def from_string(data):
        """
        Reads a Potparams from a string.

        Args:
            data (str): String containing Potparams data.

        Returns:
            Potparams object.
        """
        params = dict()
        num_params = 0

        for line in data.splitlines():
            tokens = line.split(' ')
            if num_params != 0 and len(tokens) > 1:
                # Save the parameters in the dictionary
                # with the first entry (the directory)
                # as the key
                params[tokens[0]] = tokens[1:num_params + 1]
            elif len(tokens) > 2:
                for tok in tokens:
                    if len(tok) > 2 and tok[:2] == 'p_':
                        # Params count up from 0
                        num_params = int(tok[2:]) + 1

        return Potparams(params)


class LjInteraction(object):
    """
    LjParams represents a set of Lennard-Jones potential parameters
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


class LjToDbTaskDrone(AbstractDrone):
    """
    LjToDbTaskDrone assimilates directories containing lj input to
    inserted db tasks. This drone is meant ot be used with pymatgen's
    BorgQueen to assimilate entire directory structures and insert them into
    a database using Python's multiprocessing.
    """

    #Version of this lj assimilate document.
    __version__ = "0.0.1"

    def __init__(self, host="127.0.0.1", port=27017, database="lj",
                 user=None, password=None, collection="structures",
                 tags=None, simulate_mode=False, additional_fields=None,
                 update_duplicates=True, use_full_uri=True):
        """
        Args:
            host:
                Hostname of database machine. Defaults to 127.0.0.1 or
                localhost.
            port:
                Port for db access. Defaults to mongo's default of 27017.
            database:
                Actual database to access. Defaults to "lj".
            user:
                User for db access. Requires write access. Defaults to None,
                which means no authentication.
            password:
                Password for db access. Requires write access. Defaults to
                None, which means no authentication.
            collection:
                Collection to query. Defaults to "structures".
            tags:
                A list of tags to tag the assimilated structures with.
            simulate_mode:
                Allows one to simulate db insertion without actually performing
                the insertion.
            additional_fields:
                Dict specifying additional fields to append to each doc
                inserted into the collection. For example, allows one to add
                an author or tags to a whole set of runs for example.
            update_duplicates:
                If True, if a duplicate path exists in the collection, the
                entire doc is updated. Else, duplicates are skipped.
            use_full_uri:
                Whether to use full uri path (including hostname) for the
                directory name. Defaults to True. If False, only the abs
                path will be used.
        """
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.collection = collection
        self.port = port
        self.tags = tags
        self.simulate = simulate_mode
        self.additional_fields = {} if not additional_fields \
            else additional_fields
        self.update_duplicates = update_duplicates
        self.use_full_uri = use_full_uri
        if not simulate_mode:
            conn = MongoClient(self.host, self.port, j=True)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)

    def assimilate(self, path):
        """
        Parses spl lj runs. Then insert the result into the db. and return the
        ObjectId or doc of the insertion.

        Returns:
            If in simulate_mode, the entire doc is returned for debugging
            purposes. Else, only the ObjectId of the inserted doc is returned.
        """
        tids = []
        try:
            for paramfile in glob.glob(os.path.join(path, "*.potparams")):
                potparams = Potparams.from_file(paramfile)
                for dir, params_entry in potparams.params.iteritems():
                    abs_dir = os.path.join(path, dir)
                    params = get_lj_interactions(params_entry)
                    for resfile in glob.glob(os.path.join(abs_dir, "*.res")):
                        d = self.generate_doc(resfile)
                        tid = self._insert_doc(d, params)
                        tids.append(tid)

            return tids
        except Exception as ex:
            import traceback

            print
            traceback.format_exc(ex)
            logger.error(traceback.format_exc(ex))
            return False

    @classmethod
    def generate_doc(cls, resfile):
        """
        Read the resfile and populate with the information we want to store in
        the db.
        """
        d = cls.process_res(resfile)

        d["schema_version"] = LjToDbTaskDrone.__version__

        cls.post_process(resfile, d)

        return d

    @classmethod
    def process_res(cls, resfile):
        fullpath = os.path.abspath(resfile)
        res     = Res.from_file(resfile)
        d = res.to_dict
        d["file_name"] = fullpath

        s = res.structure
        # Set the composition and formulas for the system
        comp = s.composition
        el_amt = s.composition.get_el_amt_dict()
        d.update({"unit_cell_formula": comp.to_dict,
                  "reduced_cell_formula": comp.to_reduced_dict,
                  "elements": list(el_amt.keys()),
                  "nelements": len(el_amt),
                  "pretty_formula": comp.reduced_formula,
                  "anonymous_formula": comp.anonymized_formula,
                  "nsites": comp.num_atoms,
                  "chemsys": "-".join(sorted(el_amt.keys()))})
        #d["density"] = s.density

        # Figure out the symmetry group
        sg = SymmetryFinder(s, cls.normalised_symmetry_precision(s))
        d["spacegroup"] = {"symbol": unicode(sg.get_spacegroup_symbol(),
                                             errors="ignore"),
                           "number": sg.get_spacegroup_number(),
                           "point_group": unicode(sg.get_point_group(),
                                                  errors="ignore"),
                           "source": "spglib",
                           "crystal_system": sg.get_crystal_system(),
                           "hall": sg.get_hall()}

        return d

    @classmethod
    def post_process(cls, resfile, d, use_full_uri=True):

        fullpath = os.path.abspath(resfile)

        #Convert to full uri path.
        if use_full_uri:
            d["file_name"] = get_uri(resfile)

        logger.info("Post-processed " + fullpath)


    @classmethod
    def length_per_site(cls, structure):
        return (structure.volume / structure.num_sites) ** 0.5

    @classmethod
    def normalised_symmetry_precision(cls, structure, precision=0.01):
        return precision * cls.length_per_site(structure)


    def _insert_doc(self, d, params):
        if not self.simulate:
            # Perform actual insertion into db. Because db connections cannot
            # be pickled, every insertion needs to create a new connection
            # to the db.
            conn = MongoClient(self.host, self.port)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            coll = db[self.collection]

            #params_id = self._get_params_id(db, params)
            #if params_id is None:
            #    logger.info("Failed to read params for {}, skipping.".format(d["file_name"]))
            #    return
            d["potential"] = {"name": "lennard_jones", "params": params}

            d["last_updated"] = datetime.datetime.today()
            result = coll.find_one({"file_name": d["file_name"]},
                                   fields=["file_name", "_id"])
            if result is None or self.update_duplicates:
                #d["params_id"] = params_id
                d["tags"] = self.tags

                if result is not None and "_id" in result:
                    id = result["_id"]
                else:
                    id = ObjectId()

                coll.update({"file_name": d["file_name"]}, {"$set": d},
                            upsert=True)

                logger.info("Inserting {} with _id = {}"
                            .format(d["file_name"], id))

                return str(id)
            else:
                logger.info("Skipping duplicate {}".format(d["file_name"]))
        else:
            logger.info("Simulated insert into database for {}"
                        .format(d["file_name"]))
            return d

    def _get_params_id(self, db, params):
        if not self.simulate:
            params_collection = db["params"]
            entry = params_collection.find_one({"params": params})
            if entry is None:
                id = params_collection.insert({"params": params})
            else:
                id = entry["_id"]
            return id
        else:
            return 0

    def get_valid_paths(self, path):
        """
        A Lennard-Jones run can be assimilated from a potparam file.
        """
        (parent, subdirs, files) = path
        if len(glob.glob(os.path.join(parent, "*.potparams"))) > 0:
            return [parent]
        return []

    def __str__(self):
        return "LjToDbTaskDrone"

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])

    @property
    def to_dict(self):
        init_args = {"host": self.host, "port": self.port,
                     "database": self.database, "user": self.user,
                     "password": self.password,
                     "collection": self.collection,
                     "simulate_mode": self.simulate,
                     "additional_fields": self.additional_fields,
                     "update_duplicates": self.update_duplicates}
        output = {"name": self.__class__.__name__,
                  "init_args": init_args, "version": __version__}
        return output


def get_uri(path):
    """
    Returns the URI path for a directory. This allows files hosted on
    different file servers to have distinct locations.

    Args:
        path:
            A relative path on the filesystem.

    Returns:
        Full URI path, e.g., fileserver.host.com:/full/path/of/path.
    """
    fullpath = os.path.abspath(path)
    try:
        hostname = socket.gethostbyaddr(socket.gethostname())[0]
    except:
        hostname = socket.gethostname()
    return "{}:{}".format(hostname, fullpath)


def get_lj_interactions(params):
    """
    Returns the dictionary representing the Lennard-Jones interactions from an
    array of strings.

    Args:
        params:
            An array containing Lennard-Jones parameters in the format:
             [ A~A, epsilon, sigma, m, n, cut,
              A~B, epsilon, sigma, m, n, cut,
              etc...
            ]

    Returns:
        A dictionary containing the Lennard-Jones parameters
    """
    # The old version of the lj potparams file used to have two more
    # entries that are superfluous
    leftover = len(params) % 6
    if not (leftover == 0 or leftover == 2):
        return

    lj = dict()
    for i in range(0, len(params) - leftover, 6):
        lj[params[i]] = LjInteraction.from_array(params[i + 1:i + 6]).to_dict
    return lj