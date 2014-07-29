"""
This module defines a set of functions that make querying the Lennard Jones
database easier and groups common tasks.
"""

__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Jul 16, 2014"


class HalfOpenInterval(object):
    """
    Represents a half-open interval [start,end).
    Start and end do not have to be numeric types.

    This class was taken from here:

    http://code.activestate.com/recipes/576816-interval/

    and adapted.  It is covered by the MIT license.

    The MIT License (MIT)

    Copyright (c) 2014 John Reid

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
    """

    def __init__(self, start, end):
        """Construct, start must be <= end."""
        if start > end:
            raise ValueError('Start (%s) must not be greater than end (%s)' % (start, end))
        self._start = start
        self._end = end

    start = property(fget=lambda self: self._start, doc="The interval's start")
    end = property(fget=lambda self: self._end, doc="The interval's end")

    def __str__(self):
        """As string."""
        return '[%s,%s)' % (self.start, self.end)

    def __repr__(self):
        """String representation."""
        return '[%s,%s)' % (self.start, self.end)

    def __cmp__(self, other):
        """Compare."""
        if None == other:
            return 1
        start_cmp = cmp(self.start, other.start)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.end, other.end)

    def __hash__(self):
        """Hash."""
        return hash(self.start) ^ hash(self.end)

    def intersection(self, other):
        """Intersection. @return: An empty intersection if there is none."""
        if self > other:
            other, self = self, other
        if self.end <= other.start:
            return HalfOpenInterval(self.start, self.start)
        return HalfOpenInterval(other.start, self.end)

    def hull(self, other):
        """@return: Interval containing both self and other."""
        if self > other:
            other, self = self, other
        return HalfOpenInterval(self.start, other.end)

    def overlap(self, other):
        """@return: True iff self intersects other."""
        if self > other:
            other, self = self, other
        return self.end > other.start

    def __contains__(self, item):
        """@return: True iff item in self."""
        return self.start <= item < self.end

    def zero_in(self):
        """@return: True iff 0 in self."""
        return self.start <= 0 < self.end

    def subset(self, other):
        """@return: True iff self is subset of other."""
        return self.start >= other.start and self.end <= other.end

    def proper_subset(self, other):
        """@return: True iff self is proper subset of other."""
        return self.start > other.start and self.end < other.end

    def empty(self):
        """@return: True iff self is empty."""
        return self.start == self.end

    def singleton(self):
        """@return: True iff self.end - self.start == 1."""
        return self.end - self.start == 1

    def separation(self, other):
        """@return: The distance between self and other."""
        if self > other:
            other, self = self, other
        if self.end > other.start:
            return 0
        else:
            return other.start - self.end


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
        if other is None:
            return 1
        start_cmp = cmp(self.first, other.first)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.second, other.second)


class LennardJonesSearchRange(object):
    class InteractionRange(object):
        """
        Object for representing a range of parameters in the Lennard Jones
        parameter space.  Each of the parameters can either be:
        HalfOpenInterval    - match parameter values in this interval
        None                - this parameter can take any value
        (string convertible)- match exactly this value
        """

        def __init__(self, species_pair, epsilon=None, sigma=None, m=None, n=None, cut=None):
            self.species_pair = species_pair
            self.epsilon = epsilon
            self.sigma = sigma
            self.m = m
            self.n = n
            self.cut = cut

        def to_dict(self):
            if self.species_pair is None:
                return

            criteria = dict()

            # Epsilon
            condition = self._param_to_condition(self.epsilon)
            if condition is not None:
                criteria["potential.params.{}.epsilon".format(self.species_pair)] = condition

            # Sigma
            condition = self._param_to_condition(self.sigma)
            if condition is not None:
                criteria["potential.params.{}.sigma".format(self.species_pair)] = condition

            # M
            condition = self._param_to_condition(self.m)
            if condition is not None:
                criteria["potential.params.{}.m".format(self.species_pair)] = condition

            # N
            condition = self._param_to_condition(self.n)
            if condition is not None:
                criteria["potential.params.{}.n".format(self.species_pair)] = condition

            # Cut
            condition = self._param_to_condition(self.cut)
            if condition is not None:
                criteria["potential.params.{}.cut".format(self.species_pair)] = condition

            return criteria

        def _param_to_condition(self, param):
            if param is None:
                return None
            elif isinstance(param, HalfOpenInterval):
                return {"$gte": param.start, "$lt": param.end}
            else:
                return param

    def __init__(self):
        self.interactions = dict()

    def add_interaction(self, inter):
        self.interactions[inter.species_pair] = inter

    def to_dict(self):
        c = dict()
        for interaction in self.interactions.itervalues():
            criteria = interaction.to_dict()
            if criteria is not None:
                c.update(criteria)

        return c
