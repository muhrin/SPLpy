#!/usr/bin/env python

"""
Module implementing an Res file object class.
"""



__author__ = "Martin Uhrin"
__copyright__ = "Copyright 2014, Martin Uhrin"
__version__ = "0.0.1"
__maintainer__ = "Martin Uhrin"
__email__ = "martin.uhrin.10@ucl.ac.uk"
__date__ = "Aug 12, 2014"

from abc import ABCMeta
from abc import abstractmethod
from abc import abstractproperty
import re


class Interval(object, metaclass=ABCMeta):

    def __init__(self, start, end):
        """Construct, start must be <= end."""
        if start > end:
            raise ValueError('Start (%s) must not be greater than end (%s)' % (start, end))
        self._start = start
        self._end = end

    start = property(fget=lambda self: self._start, doc="The interval's start")
    end = property(fget=lambda self: self._end, doc="The interval's end")

    @classmethod
    def from_string(cls, s):
        try:
            return Closed.from_string(s)
        except ValueError:
            pass
        try:
            return LeftClosedRightOpen.from_string(s)
        except ValueError:
            pass
        raise ValueError("Passed invalid interval format: {}.".format(s))


    @abstractproperty
    def left_closed(self):
        pass

    @abstractproperty
    def right_closed(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __hash__(self):
        pass


class Closed(Interval):
    """
    Represents a closed interval [start,end].
    Start and end do not have to be numeric types.
    """

    def __init__(self, start, end):
        super(Closed, self).__init__(start, end)

    @classmethod
    def from_string(cls, s):
        inter_patt = re.compile("\[(.+)\s*,\s*(.+)\]")
        match = inter_patt.match(s)
        if not match:
            raise ValueError("Interval string {} is not in format: (start, end)".format(s))

        return Closed(float(match.group(1)), float(match.group(2)))

    @property
    def left_closed(self):
        return False

    @property
    def right_closed(self):
        return False

    def __str__(self):
        """As string."""
        return '[%s,%s]' % (self.start, self.end)

    def __repr__(self):
        """String representation."""
        return '[%s,%s]' % (self.start, self.end)

    def __cmp__(self, other):
        """Compare."""
        if None == other or isinstance(other, LeftClosedRightOpen):
            return 1
        start_cmp = cmp(self.start, other.start)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.end, other.end)

    def __hash__(self):
        """Hash."""
        return hash("[") ^ hash(self.start) ^ hash(self.end) ^ hash("]")

class LeftClosedRightOpen(Interval):
    """
    Represents a half-open interval [start,end).
    Start and end do not have to be numeric types.
    """

    def __init__(self, start, end):
        super(LeftClosedRightOpen, self).__init__(start, end)

    @classmethod
    def from_string(cls, s):
        inter_patt = re.compile("\[(.+)\s*,\s*(.+)\)")
        match = inter_patt.match(s)
        if not match:
            raise ValueError("Interval string {} is not in format: [start, end)".format(s))

        return LeftClosedRightOpen(float(match.group(1)), float(match.group(2)))

    @property
    def left_closed(self):
        return False

    @property
    def right_closed(self):
        return True

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
        return hash("[") ^ hash(self.start) ^ hash(self.end) ^ hash(")")

    def intersection(self, other):
        """Intersection. @return: An empty intersection if there is none."""
        if self > other:
            other, self = self, other
        if self.end <= other.start:
            return LeftClosedRightOpen(self.start, self.start)
        return LeftClosedRightOpen(other.start, self.end)

    def hull(self, other):
        """@return: Interval containing both self and other."""
        if self > other:
            other, self = self, other
        return LeftClosedRightOpen(self.start, other.end)

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