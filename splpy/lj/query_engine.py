"""
This module provides a QueryEngine that simplifies queries for Mongo databases
generated using hive.
"""

from __future__ import division

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.2"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "July 23, 2014"


import matgendb as mgdb
from matgendb.query_engine import QueryEngine
from matgendb.query_engine import QueryListResults

from pymatgen import Structure, Composition
from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry


class LjQueryEngine(QueryEngine):
    """This class defines a QueryEngine interface to a Mongo Collection based on
    a set of aliases. This query engine also provides convenient translation
    between various pymatgen objects and database objects.

    This engine has special functions to make interfacing with the
    Lennard-Jones collection easier.

    The major difference between the QueryEngine's query() method and pymongo's
    find() method is the treatment of nested fields. QueryEngine's query
    will map the final result to a root level string, while pymmongo will
    return the doc as is. For example, let's say you have a document
    that is of the following form::

        {"a": {"b" : 1}}

    Using pymongo.find({}, fields=["a.b"]), you will get a doc where you need
    to do doc["a"]["b"] to access the final result (1). Using
    QueryEngine.query(properties=["a.b"], you will obtain a result that can be
    accessed simply as doc["a.b"].
    """

    def __init__(self, host="127.0.0.1", port=27017, database="lj",
                 user=None, password=None, collection="structures",
                 aliases_config=None, default_properties=None,
                 connection=None, **ignore):
        """Constructor.

        Args:
            host:
                Hostname of database machine. Defaults to 127.0.0.1 or
                localhost.
            port:
                Port for db access. Defaults to mongo's default of 27017.
            database:
                Actual database to access. Defaults to "lj".
            user:
                User for db access. Defaults to None, which means no
                authentication.
            password:
                Password for db access. Defaults to None, which means no
                authentication.
            collection:
                Collection to query. Defaults to "structures".
            connection:
                If given, ignore 'host' and 'port' and use existing connection.
            aliases_config:
                An alias dict to use. Defaults to None, which means the default
                aliases defined in "aliases.json" is used. The aliases config
                should be of the following format::

                    {
                        "aliases": {
                            "e_above_hull": "analysis.e_above_hull",
                            "energy": "output.final_energy",
                            ....
                        },
                        "defaults": {
                            "state": "successful"
                        }
                    }

                The "aliases" key defines mappings, which makes it easier to
                query for certain nested quantities. While Mongo does make it
                easy to map collections, it is sometimes beneficial to
                organize the doc format in a way that is different from the
                query format.
                The "defaults" key specifies criteria that should be applied
                by default to all queries. For example, a collection may
                contain data from both successful and unsuccessful runs but
                for most querying purposes, you may want just successful runs
                only. Note that defaults do not affect explicitly specified
                criteria, i.e., if you suppy a query for {"state": "killed"},
                this will override the default for {"state": "successful"}.
            default_properties:
                List of property names (strings) to use by default, if no
                properties are given to the 'properties' argument of
                query().
        """
        super(LjQueryEngine, self).__init__(host, port, database, user, password,
                                       collection, aliases_config,
                                       default_properties, connection)
        self.params = self.db["params"]

    def get_entries(self, criteria, inc_structure=False, optional_data=None):
        """
        Get ComputedEntries satisfying a particular criteria.

        Args:
            criteria:
                Criteria obeying the same syntax as query.
            inc_structure:
                Optional parameter as to whether to include a structure with
                the ComputedEntry. Defaults to False. Use with care - including
                structures with a large number of entries can potentially slow
                down your code to a crawl.
            optional_data:
                Optional data to include with the entry. This allows the data
                to be access via entry.data[key].

        Returns:
            List of pymatgen.entries.ComputedEntries satisfying criteria.
        """
        all_entries = list()
        optional_data = [] if not optional_data else list(optional_data)
        fields = [k for k in optional_data]
        fields.extend(["_id", "unit_cell_formula", "energy"])
        for c in super(LjQueryEngine, self).query(fields, criteria):
            optional_data = {k: c[k] for k in optional_data}
            if inc_structure:
                struct = Structure.from_dict(c["structure"])
                entry = ComputedStructureEntry(struct, c["energy"],
                                               data=optional_data,
                                               entry_id=c["_id"])
            else:
                entry = ComputedEntry(Composition(c["unit_cell_formula"]),
                                      c["energy"], data=optional_data,
                                      entry_id=c["_id"])
            all_entries.append(entry)

        return all_entries

    def get_param_ids(self, crit, properties=None):
        cur = self.params.find(crit, fields=properties)
        return mgdb.query_engine.QueryListResults(None, cur)

    def get_param_id_criteria(self, params_range):
        # Add the set of parameter ids to look for
        criteria = dict()
        param_ids = self.get_param_ids(params_range.to_criteria())

        # TODO: Put in special case for just 1 parameter
        criteria["potential.params_id"] =\
            {"$in": [k["_id"] for k in param_ids]}

        return criteria

    def query_in_params_range(self, params_range, properties=None, criteria=None):
        if criteria is None:
            criteria = dict()

        # Add the set of parameter ids to look for
        criteria.update(self.get_param_id_criteria(params_range))

        return super(LjQueryEngine, self).query(properties, criteria)

    def query_at_each_param_point(self, params_range, properties=None, criteria=None):

        points = list()
        for value in self.get_param_ids(params_range.to_criteria()):
            crit = criteria if criteria is not None else dict()
            crit["potential.params_id"] = value["_id"]
            points.append((value["_id"], super(LjQueryEngine, self).query(properties, crit)))

        return points
