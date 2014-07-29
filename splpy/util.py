#!/usr/bin/env python

"""
Module for helper functions and other things that don't fit neatly elsewhere.
"""


def find_or_create(collection, query, update, **kwargs):
    return collection.find_and_modify(query, {"$setOnInsert": update}, upsert=True, new=True, **kwargs)