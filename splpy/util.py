#!/usr/bin/env python

"""
Module for helper functions and other things that don't fit neatly elsewhere.
"""


def find_or_create(collection, query, update, **kwargs):
    return collection.find_and_modify(query, {"$setOnInsert": update}, upsert=True, new=True, **kwargs)


class OrderedPair(object):
    """
    Object for representing a pair of objects (usually numbers) that are
    accessed in order using the .first and .second attributes

    """
    def __init__(self, x0=None, x1=None):
        self._first = None
        self._second = None
        self.set(x0, x1)

    @property
    def first(self):
        return self._first

    @property
    def second(self):
        return self._second

    def set(self, x0, x1):
        if x0 < x1:
            self._first = x0
            self._second = x1
        else:
            self._first = x1
            self._second = x0

    def __str__(self):
        return "{}~{}".format(self.first, self.second)

    def __hash__(self):
        return hash(self.first) ^ hash(self.second)

    def __cmp__(self, other):
        """Compare."""
        if other is None or not isinstance(other, OrderedPair):
            return 1
        start_cmp = cmp(self.first, other.first)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.second, other.second)

    @classmethod
    def from_string(cls, s):
        values = s.split("~")
        if len(values) != 2:
            raise ValueError("OrderedPair expects format A~B")
        return OrderedPair(values[0], values[1])