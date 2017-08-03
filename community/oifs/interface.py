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

    # set the experiment name, used as part of the input file names
    @remote_function
    def set_exp_name(str='a'):
        pass

    # Inserts grid, not yet implemented
    @remote_function
    def commit_grid():
        pass

    # Utility method returning the full u-component flow field
    @remote_function(must_handle_array=True)
    def get_field_Pfull_(i = 0,k = 0):
        returns (out = 0. | units.Pa)

    # Utility method returning the full u-component flow field
    @remote_function(must_handle_array=True)
    def get_field_Phalf_(i = 0,k = 0):
        returns (out = 0. | units.Pa)

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

    # Utility method returning the cloud coverage field
    @remote_function(must_handle_array=True)
    def get_field_A_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the ozone field
    @remote_function(must_handle_array=True)
    def get_field_O3_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the full u-component tendency
    @remote_function(must_handle_array=True)
    def get_tendency_U_(i = 0,k = 0):
        returns (out = 0. | units.m / units.s)

    # Utility method returning the full v-component tendency
    @remote_function(must_handle_array=True)
    def get_tendency_V_(i = 0,k = 0):
        returns (out = 0. | units.m / units.s)

    # Utility method returning the full temperature tendency
    @remote_function(must_handle_array=True)
    def get_tendency_T_(i = 0,k = 0):
        returns (out = 0. | units.K)

    # Utility method returning the specific humidity tendency
    @remote_function(must_handle_array=True)
    def get_tendency_SH_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the liquid water content tendency
    @remote_function(must_handle_array=True)
    def get_tendency_QL_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the ice water content tendency
    @remote_function(must_handle_array=True)
    def get_tendency_QI_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the cloud coverage tendency
    @remote_function(must_handle_array=True)
    def get_tendency_A_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method returning the ozone tendency
    @remote_function(must_handle_array=True)
    def get_tendency_O3_(i = 0,k = 0):
        returns (out = 0.)

    # Utility method setting the u-component tendency
    @remote_function(must_handle_array=True)
    def set_tendency_U_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the v-component tendency
    @remote_function(must_handle_array=True)
    def set_tendency_V_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the temperature tendency
    @remote_function(must_handle_array=True)
    def set_tendency_T_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the specific humidity tendency
    @remote_function(must_handle_array=True)
    def set_tendency_SH_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the liquid water tendency
    @remote_function(must_handle_array=True)
    def set_tendency_QL_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the ice water tendency
    @remote_function(must_handle_array=True)
    def set_tendency_QI_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the cloud coverage tendency
    @remote_function(must_handle_array=True)
    def set_tendency_A_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting the ozone tendency
    @remote_function(must_handle_array=True)
    def set_tendency_O3_(i = 0,k = 0,v = 0.):
        pass

    # Utility method setting a superparametrization mask
    @remote_function(must_handle_array=True)
    def set_mask(i = 0):
        pass

    # Utility method resetting a superparametrization mask
    @remote_function(must_handle_array=True)
    def reset_mask(i = 0):
        pass

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

    # Single step model evolution
    @remote_function
    def evolve_model_single_step():
        pass

    # Single step model evolution, returns when cloud scheme is reached
    @remote_function
    def evolve_model_until_cloud_scheme():
        pass

    # Single step model evolution, only executes the IFS cloud scheme
    @remote_function
    def evolve_model_cloud_scheme():
        pass

    # Single step model evolution, enters time step function after cloud scheme
    @remote_function
    def evolve_model_from_cloud_scheme():
        pass

# OpenIFS class implementation
class OpenIFS(CommonCode):

    inputfile = "fort.4"
    backupfile = "fort.4.bkp"

    # Constructor
    def __init__(self,**options):
        print('__init__')

        self.stopping_conditions = StoppingConditions(self)
        self.itot = 0
        self.ktot = 0
        self.latitudes = None
        self.longitudes = None

        if(not os.path.exists(OpenIFS.inputfile)):
            print('inputfile:' + OpenIFS.inputfile)
            print('cwd:' + os.getcwd())
            raise Exception("File fort.4 not found. Creating an openIFS model from scratch is not supported yet.")
        else:
            os.rename(OpenIFS.inputfile,OpenIFS.backupfile)
            self.params = f90nml.read(OpenIFS.backupfile)
            self.patch = {"NAMPAR0":{"NPROC":options.get("number_of_workers",1)}}
            f90nml.patch(OpenIFS.backupfile,self.patch,OpenIFS.inputfile)
            print('***** Done patching ', OpenIFS.inputfile)
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
        object.add_transition("RUN","EVOLVED","evolve_model_single_step",False)

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
        object.add_method("EVOLVED","evolve_model_single_step")
        object.add_method("EVOLVED","evolve_model_until_cloud_scheme")
        object.add_method("EVOLVED","evolve_model_cloud_scheme")
        object.add_method("EVOLVED","evolve_model_from_cloud_scheme")

    def get_field(self,fid,i,k):
        if(fid == "Pfull"):
            return self.get_field_Pfull_(i,k)
        elif(fid == "Phalf"):
            return self.get_field_Phalf_(i,k)
        elif(fid == "U"):
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
        elif(fid == "A"):
            return self.get_field_A_(i,k)
        elif(fid == "O3"):
            return self.get_field_O3_(i,k)
        else:
            raise Exception("Unknown atmosphere prognostic field identidier:",fid)

    def set_tendency(self,fid,i,k,v):
        if(fid == "U"):
            return self.set_tendency_U_(i,k,v)
        elif(fid == "V"):
            return self.set_tendency_V_(i,k,v)
        elif(fid == "T"):
            return self.set_tendency_T_(i,k,v)
        elif(fid == "SH"):
            return self.set_tendency_SH_(i,k,v)
        elif(fid == "A"):
            return self.set_tendency_A_(i,k,v)
        elif(fid == "O3"):
            return self.set_tendency_O3_(i,k,v)
        else:
            raise Exception("Unknown atmosphere tendency field identidier:",fid)

    def get_volume_field(self,fid):
        ktop = (self.ktot + 1) if fid == "Phalf" else self.ktot
        pts = numpy.mgrid[0:self.itot,0:ktop]
        i,k = pts.reshape(2,-1)
        return self.get_field(fid,i,k).reshape((self.itot,ktop))

    def get_layer_field(self,fid,layerindex):
        i,k = numpy.arange(0,self.itot),numpy.full([self.itot],layerindex)
        return self.get_field(fid,i,k)

    def get_profile_field(self,fid,colindex):
        ktop = (self.ktot + 1) if fid == "Phalf" else self.ktot
        i,k = numpy.full([ktop],colindex),numpy.arange(0,ktop)
        return self.get_field(fid,i,k)

    # like get_profile_field, but gets several columns at once
    def get_profile_fields(self,fid,colindex):
        ktop = (self.ktot + 1) if fid == "Phalf" else self.ktot

        i,k = numpy.meshgrid(colindex,numpy.arange(ktop),indexing='ij')
        i,k = i.reshape(-1), k.reshape(-1) #i = a,a,a,...,b,b,b,...   if colindex = [a, b, ...] 
                                           #k = 1,2,3,...,1,2,3,...

        #print('get_profile_fields', colindex)
        #print(i)
        #print(k)
        f = self.get_field(fid,i,k)
        f = f.reshape((len(colindex),-1))
        #print (f)
        return f

    def get_tendency(self,fid,i,k):
        if(fid == "U"):
            return self.get_tendency_U_(i,k)
        elif(fid == "V"):
            return self.get_tendency_V_(i,k)
        elif(fid == "T"):
            return self.get_tendency_T_(i,k)
        elif(fid == "SH"):
            return self.get_tendency_SH_(i,k)
        elif(fid == "QL"):
            return self.get_tendency_SH_(i,k)
        elif(fid == "QI"):
            return self.get_tendency_SH_(i,k)
        elif(fid == "A"):
            return self.get_tendency_A_(i,k)
        elif(fid == "O3"):
            return self.get_tendency_O3_(i,k)
        else:
            raise Exception("Unknown atmosphere tendency field identidier:",fid)

    def get_volume_tendency(self,fid):
        ktop = self.ktot
        pts = numpy.mgrid[0:self.itot,0:ktop]
        i,k = pts.reshape(2,-1)
        return self.get_tendency(fid,i,k).reshape((self.itot,ktop))

    def get_layer_tendency(self,fid,layerindex):
        i,k = numpy.arange(0,self.itot),numpy.full([self.itot],layerindex)
        return self.get_tendency(fid,i,k)

    def get_profile_tendency(self,fid,colindex):
        ktop = (self.ktot + 1) if fid == "Phalf" else self.ktot
        i,k = numpy.full([ktop],colindex),numpy.arange(0,ktop)
        return self.get_tendency(fid,i,k)

    def set_profile_tendency(self,fid,colindex, v):
        ktop = (self.ktot + 1) if fid == "Phalf" else self.ktot
        i,k = numpy.full([ktop],colindex),numpy.arange(0,ktop)
        self.set_tendency(fid,i,k,v)
