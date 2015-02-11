"""
This module contains functions used for creating maps from the database
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jan 26, 2015"

import splpy.lj.db_hulls as db_hulls


def create_property_getter(property_name):
    if property_name == "prototype":
        return Prototype()
    elif property_name == "stable_stoichiometries":
        return StableStoichiometries()
    else:
        return SimpleProperty(property_name)


def create_mask_property_getter(property_name):
    if property_name == "stable":
        return StableStructures()

    return None


class SimpleProperty(object):
    def __init__(self, property_name):
        self._property = property_name

    @property
    def properties(self):
        return [self._property]

    def get_values(self, structure_docs, query_engine):
        return {doc['_id']: doc[self._property] for doc in structure_docs}


class Prototype(object):
    @property
    def properties(self):
        return ["prototype_id", "spacegroup.symbol"]

    def get_values(self, structure_docs, query_engine):
        protos = query_engine.db["prototypes"]

        docs = list(structure_docs)
        proto_ids = set(doc["prototype_id"] for doc in docs)
        proto_ids.discard(None)

        proto_values = dict()
        for proto_id in proto_ids:
            value = None
            doc = protos.find({"_id": proto_id})[0]
            if "structure_type" in doc:
                value = doc["structure_type"]
            else:
                if "strukturbericht" in doc:
                    value = doc["structure_type"]
                elif "pearson" in doc:
                    value = "{} ({})".format(doc["spacegroup"]["symbol"], doc["pearson"])

            if value:
                proto_values[proto_id] = value

        # Set the values for each of the parameter points
        values = dict()
        for doc in docs:
            # Use the spacegroup as a fallback
            values[doc["_id"]] = proto_values[doc["prototype_id"]] \
                if doc["prototype_id"] in proto_values else doc["spacegroup.symbol"]
        return values


class StableStoichiometries(object):
    @property
    def properties(self):
        return ["potential.params_id"]

    def get_values(self, structure_docs, query_engine):
        values = dict()
        for doc in structure_docs:
            stoichs = db_hulls.get_stable_stoichiometries(query_engine.db, doc["potential.params_id"])
            if stoichs:
                values[doc['_id']] = "+".join(sorted(stoichs))
        return values


class StableStructures(object):
    @property
    def properties(self):
        return []

    def get_values(self, structure_docs, query_engine):
        values = dict()
        for doc in structure_docs:
            stable = db_hulls.is_stable(query_engine.db, doc["_id"])
            if stable is not None:
                values[doc["_id"]] = stable
        return values