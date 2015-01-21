"""
This module contains code relating to convex hulls and their storage
in the database
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.1.0"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jan 20, 2015"

import splpy.util as util


def _regenerate_hull(db, hull_id):
    # Get all the structures at the parameter point of this hull

    # Save them to a temporary directory

    # Calculate formation enthalpies and distance above hull

    # Save the new entries

    pass

def _update_hull(db, params_id):
    hulls = db['hulls']

    # Find or create the hull entry and find the id
    hull = util.find_or_create(hulls, {'params_id': params_id}, {'params_id': params_id},
                               fields={'last_updated': 1, 'entries.structure_id': 1})

    # Check if we need to remove any structures from the hull because they're no longer
    # in the structures collection
    structures = db['structures']
    hull_changed = False
    structure_ids = list()
    for entry in hull['entries']:
        cur = structures.find({'_id': entry['_id']}).limit(1)
        if cur.size() == 0:
            hulls.update({'_id': hull['_id']},
                         {'$pull': {'entries.structure_id': entry['_id']}})
            hull_changed = True
        else:
            structure_ids.append(entry['_id'])

    if not hull_changed and structures.find({'potential.params_id': params_id, '_id': {'$nin': structure_ids}}).limit(
            1).size() > 0:
        hull_changed = True

    if hull_changed:
        _regenerate_hull(db, hull['_id'])
