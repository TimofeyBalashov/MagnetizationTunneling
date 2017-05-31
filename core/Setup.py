from functools import wraps

def setupmethod(f):
    """Use this decorator for methods that change the internal state / parameters
       of the system, such that the eigenstates have to be recalculated.

       This method will reset the 'data was changed' flag"""
    @wraps(f)
    def decorated(self, *ops, **kwops):
        if self.ready:
            self.makeNotReady()
        return f(self, *ops, **kwops)
    return decorated


def buildmethod(f):
    """Use this decorator for methods that only need to be run once for a certain
       set of parameters.

       This method will only be run if data was changed.
       The 'data was changed' flag will remain unchanged."""
    @wraps(f)
    def decorated(self, *ops, **kwops):
        if not self.ready:
            return f(self, *ops, **kwops)
    return decorated


def resultmethod(f):
    """Use this decorator for methods that use the eigenstates to
       calculate their resutls.

       This method will make the system ready and set the
       'data was changed'"""
    @wraps(f)
    def decorated(self, *ops, **kwops):
        self.makeReady()
        return f(self, *ops, **kwops)
    return decorated


class SetupClass:
    """A class for building object trees, where a change to a leaf
    leads to a change in all parent nodes up to the root.

    The change is normally postponed to avoid repeated cascading
    update propagation until access to the state-dependent data
    of a parent node is required. The instance variable *ready*
    is used to indicate pending updates. Changes to this variable
    should be made through instance methods `makeReady()` and
    `makeNotReady()`.

    Subclasses should implement the `_build()` method that updates
    the object state as necessary based on children state.
    """
    def __init__(self, parent=None):
        self.ready = False
        self.parent = parent

    def makeReady(self):
        """Ensures that the object state is up-to-date.
        Will call `_build()` if necessary.
        """
        if not self.ready:
            self._build()
            self.ready = True

    def makeNotReady(self):
        """Sets 'update pending' status for this node
        and all parent nodes.
        """
        self.ready = False
        if self.parent is not None:
            self.parent.makeNotReady()
