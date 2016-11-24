import os
import numpy
import f90nml
import shutil

from omuse.units import units
from amuse.community.interface.common import CommonCodeInterface,CommonCode
from amuse.support.interface import InCodeComponentImplementation
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.rfi.core import CodeInterface,LegacyFunctionSpecification
from amuse.rfi.core import legacy_function,remote_function
from amuse.units.core import system,no_system
from amuse.community.interface.stopping_conditions import StoppingConditionInterface,StoppingConditions
from amuse import datamodel

class OpenIFSInterface(CodeInterface,
                       CommonCodeInterface,
                       StoppingConditionInterface,
                       LiteratureReferencesMixIn):
    """

    OpenIFS: The Open Integral Forecasting System

    """
    use_modules=["StoppingConditions","openifs_interface"]

    # Constructor
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    # Determines worker binary name
    def name_of_the_worker(self):
        return "openifs_worker"

    # Returns model internal time in sec.
    @remote_function
    def get_model_time():
        returns (time = 0. | units.s)

    # Returns the model time step in sec.
    @remote_function
    def get_timestep():
        returns (dt = 0. | units.s)

    # Inserts grid, not yet implemented
    @remote_function
    def commit_grid():
        pass

    # Utility method returning the full u-component flow field
    @remote_function(must_handle_array=True)
    def get_field_U_(i = 0,k = 0):
        returns (out = 0. | units.m / units.s)

    # Utility method returning the full v-component flow field
    @remote_function(must_handle_array=True)
    def get_field_V_(i = 0,k = 0):
        returns (out = 0. | units.m / units.s)

    # Utility method returning the full temperature field
    @remote_function(must_handle_array=True)
    def get_field_T_(i = 0,k = 0):
        returns (out = 0. | units.K)

    # Utility method returning the specific humidity field
    @remote_function(must_handle_array=True)
    def get_field_SH_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the liquid water content field
    @remote_function(must_handle_array=True)
    def get_field_QL_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the ice water content field
    @remote_function(must_handle_array=True)
    def get_field_QI_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the ozone field
    @remote_function(must_handle_array=True)
    def get_field_O3_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the gridpoint latitudes and longitudes.
    @remote_function(must_handle_array=True)
    def get_gridpoints(i = 0):
        returns (lats = 0. | units.rad,lons = 0. | units.rad)

    # Returns the horizontal and vertical index ranges
    @remote_function
    def get_grid_sizes():
        returns (i=1,k=1)

    # Model evolution function
    @remote_function
    def evolve_model(tend = 0. | units.s):
        pass

# OpenIFS class implementation
class OpenIFS(CommonCode):

    inputfile = "fort.4"
    backupfile = "fort.4.bkp"

    # Constructor
    def __init__(self,**options):
        self.stopping_conditions = StoppingConditions(self)
        self.itot = 0
        self.ktot = 0
        self.latitudes = None
        self.longitudes = None
        if(not os.path.exists(OpenIFS.inputfile)):
            raise Exception("File fort.4 not found. Creating an openIFS model from scratch is not supported yet.")
        else:
            os.rename(OpenIFS.inputfile,OpenIFS.backupfile)
            self.params = f90nml.read(OpenIFS.backupfile)
            self.patch = {"NAMPAR0":{"NPROC":options.get("number_of_workers",1)}}
            f90nml.patch(OpenIFS.backupfile,self.patch,OpenIFS.inputfile)
        CommonCode.__init__(self,OpenIFSInterface(**options),**options)

    # Commit run parameters, not implemented yet
    def commit_parameters(self):
        self.overridden().commit_parameters()

    # Commit grid, reads grid size and fills geometry cache
    def commit_grid(self):
        self.itot,self.ktot = self.get_grid_sizes()
        indices = numpy.array([i for i in numpy.arange(0,self.itot)])
        self.latitudes,self.longitudes = self.get_gridpoints(indices)
        self.overridden().commit_grid()

    def cleanup_code(self):
        if(os.path.exists(OpenIFS.backupfile) and os.path.exists(OpenIFS.inputfile)):
            os.remove(OpenIFS.inputfile)
        if(os.path.exists(OpenIFS.backupfile)):
            os.rename(OpenIFS.backupfile,OpenIFS.inputfile)
        self.overridden().cleanup_code()

    # State machine definition
    def define_state(self, object):
        object.set_initial_state("UNINITIALIZED")

        object.add_transition("UNINITIALIZED","INITIALIZED","initialize_code")
        object.add_transition("!UNINITIALIZED!STOPPED","END","cleanup_code")
        object.add_transition("END","STOPPED","stop",False)
        object.add_transition("INITIALIZED","EDIT","commit_parameters")
        object.add_transition("EDIT","RUN","commit_grid")
        object.add_transition("RUN","EVOLVED","evolve_model",False)

        object.add_method("INITIALIZED","before_set_interface_parameter")
        object.add_method("UNINITIALIZED","before_get_parameter")
        object.add_method("!UNINITIALIZED","before_set_parameter")
        object.add_method("END","before_get_parameter")
        object.add_method("STOPPED","stop")
        for state in ["RUN","EDIT","EVOLVED"]:
            object.add_method(state,"get_model_time")
            object.add_method(state,"get_timestep")
            object.add_method(state,"get_volume_field")
            object.add_method(state,"get_profile_field")
        object.add_method("EVOLVED","evolve_model")

    def get_field(self,fid,i,k):
        if(fid == "U"):
            return self.get_field_U_(i,k)
        elif(fid == "V"):
            return self.get_field_V_(i,k)
        elif(fid == "T"):
            return self.get_field_T_(i,k)
        elif(fid == "SH"):
            return self.get_field_SH_(i,k)
        elif(fid == "QL"):
            return self.get_field_QL_(i,k)
        elif(fid == "QI"):
            return self.get_field_QI_(i,k)
        elif(fid == "O3"):
            return self.get_field_O3_(i,k)
        else:
            raise Exception("Unknown atmosphere prognostic field identidier:",fid)

    def get_volume_field(self,fid):
        pts = numpy.mgrid[0:self.itot,0:self.ktot]
        i,k = pts.reshape(2,-1)
        return self.get_field(fid,i,k).reshape((self.itot,self.ktot))

    def get_layer_field(self,fid,layerindex):
        i,k = numpy.arange(0,self.itot),numpy.full([self.itot],layerindex)
        return self.get_field(fid,i,k)

    def get_profile_field(self,fid,colindex):
        i,k = numpy.full([self.ktot],colindex),numpy.arange(0,self.ktot)
        return self.get_field(fid,i,k)
