"""
This module provides a QueryEngine that simplifies queries for Mongo databases
generated using hive.
"""


from __future__ import division

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "July 23, 2014"


import matgendb as mgdb


class LjQueryEngine(object):
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
        self._query = mgdb.QueryEngine(host, port, database, user, password,
                                       collection, aliases_config,
                                       default_properties, connection)

    def set_collection(self, collection):
        """
        Switch to another collection. Note that you may have to set the
        aliases and default properties via set_aliases_and_defaults if the
        schema of the new collection differs from the current collection.

        Args:
            collection:
                Name of collection.
        """
        self._query.set_collection(collection)

    def get_entries_in_system(self, species, inc_structure=False,
                              optional_data=None, additional_criteria=None):
        """
        Gets all entries in a chemical system, e.g. A-B-C will return all
        A-B, B-C, A-C, A-B-C compounds.

        Args:
            elements:
                Sequence of species symbols, e.g. ['A','B','C']
            inc_structure:
                Optional parameter as to whether to include a structure with
                the ComputedEntry. Defaults to False. Use with care - including
                structures with a large number of entries can potentially slow
                down your code to a crawl.
            optional_data:
                Optional data to include with the entry. This allows the data
                to be access via entry.data[key].
            additional_criteria:
                Added ability to provide additional criteria other than just
                the chemical system.

        Returns:
            List of ComputedEntries in the chemical system.
        """
        chemsys_list = []
        for i in range(len(species)):
            for combi in itertools.combinations(species, i + 1):
                chemsys = "-".join(sorted(combi))
                chemsys_list.append(chemsys)
        crit = {"chemsys": {"$in": chemsys_list}}
        if additional_criteria is not None:
            crit.update(additional_criteria)
        return self.get_entries(crit, inc_structure,
                                optional_data=optional_data)