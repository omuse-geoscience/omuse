import logging
import re
import numpy
import netCDF4

# File readers for DALES experiments. We support single-timestep prifile ASCII
# file and netcdf files for time-dependent output profiles or time series.

log = logging.getLogger(__name__)


# Factory method for creating file readers.
def make_file_reader(filename):
    if filename.endswith(".nc"):
        result = NetcdfReader(filename)
    else:
        result = AsciiReader(filename)
    log.info("Detected %d variables in file %s" % (len(result.variables), filename))
    return result


# Base class for the DALES file readers.
class DalesReader(object):

    def __init__(self, filepath):
        self.filepath = filepath
        self.variables = []
        self.dimensions = []
        self.units = []
        self.data = []

    # Returns the heights read from the file. If the file only contains time series,
    # returns a list with a single zero height.
    def get_heights(self, var):
        return []

    # Returns the times read from the file in seconds since the starttime.
    # If the file only contains time-independent profiles, returns a list with
    # a single zero time step.
    def get_times(self, var):
        return []

    # Returns the unit of the given variable.
    def get_unit(self, var):
        if var in self.variables:
            index = self.variables.index(var)
            return self.units[index]
        else:
            return "1"

    # Returns the shape of the given variable array.
    def get_shape(self, var):
        return []

    # Returns the missing value for the givan variable
    def get_missval(self, var):
        return None


# ASCII file reader. Currently, we only support single time step text files
# (initial profiles) with all text commented.
class AsciiReader(DalesReader):
    comment_char = '#'

    def __init__(self, filepath):
        super(AsciiReader, self).__init__(filepath)
        with open(filepath, 'r') as f:
            f.readline()
            first_row = f.readline()
        if first_row.startswith(AsciiReader.comment_char):
            cols = first_row[1:].replace("(z)", "").replace(',', ' ').split()
        else:
            cols = first_row[:].replace("(z)", "").replace(',', ' ').split()
        for col in cols:
            tokens = [s for s in re.split("\(|\)|\[|\]", col)]
            self.variables.append(tokens[0])
            self.units.append(tokens[1] if len(tokens) > 1 else "1")
        for z in ["z", "zf", "zc", "zm"]:
            if z in self.variables:
                index = self.variables.index(z)
                self.variables[index] = "height"
                break
        if "height" in self.variables:
            self.dimensions.append("height")
        self.data = numpy.loadtxt(self.filepath, skiprows=2).transpose()

    def __getitem__(self, varname, *args):
        if varname in self.variables:
            index = self.variables.index(varname)
            new_args = (index,) + args
            return self.data.__getitem__(new_args)
        else:
            return []

    def get_heights(self, var):
        return self.__getitem__("height")

    def get_times(self, var):
        return [0]

    def get_shape(self, var):
        if var in self.variables:
            return self.data.shape[1:]
        else:
            return super(AsciiReader, self).get_shape(var)


# Netcdf file reader implementation.
class NetcdfReader(DalesReader):

    def __init__(self, filepath):
        super(NetcdfReader, self).__init__(filepath)
        self.dataset = netCDF4.Dataset(filepath, 'r')
        self.variables = [v for v in self.dataset.variables if v not in self.dataset.dimensions]
        self.units = [getattr(self.dataset.variables[v], "units", None) for v in self.variables]

    def __getitem__(self, varname, *args):
        v = self.dataset.variables.get(varname, None)
        if not v:
            return []
        return v.__getitem__(*args)

    def get_heights(self, var):
        v = self.dataset.variables.get(var, None)
        if not v:
            return []
        for lev_type in ["zt", "zm"]:
            if lev_type in v.dimensions:
                dim = self.dataset.variables.get(lev_type, None)
                return dim[:] if dim else [0.]

    def get_times(self, var):
        v = self.dataset.variables.get(var, None)
        if not v:
            return []
        if "time" in v.dimensions:
            dim = self.dataset.variables.get("time", None)
            return dim[:] if dim else []

    def get_shape(self, var):
        if var in self.variables:
            return self.dataset.variables[var].shape
        else:
            return super(NetcdfReader, self).get_shape(var)

    def get_missval(self, var):
        if var in self.dataset.variables:
            return getattr(self.dataset.variables[var], "_FillValue", None)
        else:
            return super(NetcdfReader, self).get_missval(var)
