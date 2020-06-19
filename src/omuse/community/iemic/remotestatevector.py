from amuse.support.exceptions import CodeException

class RemoteStateVector(object):
    def __init__(self, interface):
        self.interface=interface
        self._id=interface._new_state()
        self.grid=interface.get_grid(self)
    def copy_to(self, state):
        self.interface._copy_state(self._id, state._id)
    def copy(self):
        new=RemoteStateVector(self.interface)
        self.copy_to(new)
        return new
    def norm(self):
        return self.interface._get_state_norm(self._id)
    def __del__(self):
        try:
          self.interface._remove_state(self._id)
        except CodeException:
          pass
        except AttributeError:
          pass
    def __mul__(self, other):
        new=self.copy()
        self.interface._mul_state(new._id, other)
        return new
    def __rmul__(self, other):
        return self.__mul__(other)
    def __truediv__(self, other):
        return self.__mul__(1./other)
    def __add__(self, other):
        new = self.copy()
        self.interface._update_state(new._id, other._id, 1.0)
        return new
    def __sub__(self, other):
        new = self.copy()
        self.interface._update_state(new._id, other._id, -1.0)
        return new
    def __neg__(self):
        new = self.copy()
        self.interface._mul_state(new._id, -1.0)
        return new
    def __getitem__(self, index):
        return self.grid[index]
