import os.path
import numpy
from amuse.test.amusetest import TestWithMPI
import time

from omuse.community.adcirc.interface import AdcircInterface,Adcirc

from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

default_options={}
#default_options = dict(redirection="none")
#default_options=dict(debugger="gdb")

import logging
logging.basicConfig(level=logging.DEBUG)
logging.getLogger("code").setLevel(logging.DEBUG)

from omuse.community.adcirc.read_grid import adcirc_grid_reader, adcirc_parameter_reader
from omuse.community.adcirc.write_grid import adcirc_grid_writer, adcirc_parameter_writer

AMIG=0.000140525700000 | units.s**-1
PER=2*numpy.pi/AMIG
EMO=0.3048 | units.m
DRAMP=2 | units.day

def ramp(t):
    return numpy.tanh(2*t/DRAMP)

def tidal_force_function(t):
    ncyc=numpy.floor(t/PER)
    aj=AMIG*(t-ncyc*PER)
    return ramp(t)*EMO*numpy.cos(aj)


def run(tend=5. | units.day, state=None):

    param=adcirc_parameter_reader("data/test/2d/fort.15")
    param.read_parameters(NETA=9)
    param.parameters['NBFR']=-1

    gr=adcirc_grid_reader("data/test/2d/fort.14")
    gr.read_grid()
    nodes,elements,elev_boundary,flow_boundary=gr.get_sets()

    code=Adcirc()

    code._parameters=param.parameters        
    code.assign_grid_and_boundary(nodes,elements,elev_boundary, flow_boundary)

    code.parameters.use_interface_elevation_boundary=True
    code.parameters.use_interface_parameters=True
    code.parameters.use_interface_grid=True
    code.parameters.A_H=param.parameters["ESLM"] | units.m**2/units.s
    code.parameters.timestep=abs(param.parameters["DTDP"]) | units.s
    code.parameters.bottom_friction_law=["linear","quadratic","hybrid"][param.parameters["NOLIBF"]]
    try:
      code.parameters.linear_bottom_friction_coeff=param.parameters["TAU"]| units.s**-1
    except:
      pass
    try:
      code.parameters.quadratic_bottom_friction_coeff=param.parameters["CF"]
    except:
      pass
    code.parameters.use_predictor_corrector=param.parameters["DTDP"]<0
    code.parameters.use_interface_met_forcing=False

    if state:
      nodes,elements=state

    if state:
      code.parameters.begin_time=nodes.collection_attributes.time

    if state:
      channel=nodes.new_channel_to(code.nodes)
      channel.copy_attributes(["eta","deta_dt","status","vx","vy"])
      channel=elements.new_channel_to(code.elements)
      channel.copy_attributes(["status"])

    tnow=code.model_time
    dt=code.parameters.timestep

    elev_boundaries= list(code.elevation_boundaries())

    eta61=[]
    time=[]
    forcing=[]

    while tnow<tend-dt/2:
        elev_boundaries[0].eta=tidal_force_function(tnow+dt/2)
        code.evolve_model(tnow+dt)
        tnow=code.get_model_time()

        eta=code.nodes[60].eta.number
        time.append(tnow.number)
        eta61.append(eta)  
        forcing.append(elev_boundaries[0].eta[0].number)
  
    state=code.nodes.copy(),code.elements.copy()
    state[0].collection_attributes.time=code.model_time
    print "done at:", code.model_time.in_(units.day)
    
    code.stop()
    
    return state,time,eta61,forcing

class TestAdcircRestart(TestWithMPI):

    def test1(self):
        from matplotlib import pyplot

        pyplot.ion()
        f=pyplot.figure(figsize=(8,6))
        pyplot.show()

        t1=time.time()
        ref_state,ref_time,ref_eta61,ref_forcing=run(tend=5. | units.day)
        t2=time.time()
    
        t3=time.time()
        state,stime,eta61,forcing=run(tend=4.5 | units.day)
        t4=time.time()
        
        t5=time.time()
        state,time_,eta61_,forcing_=run(tend=5. | units.day,state=state)
        t6=time.time()
        
        print "timing:", t2-t1,t6-t5+t4-t3
        
        stime.extend(time_)
        eta61.extend(eta61_)
        forcing.extend(forcing_)

        print ref_eta61[-1]
        print eta61[-1]

        pyplot.plot(ref_time,ref_eta61,'r')
        pyplot.plot(ref_time,ref_forcing,'g:')
        pyplot.plot(ref_time,tidal_force_function((ref_time| units.s)).number)    
        pyplot.plot(stime,eta61,'r+')
        pyplot.plot(stime,forcing,'g+')
        pyplot.draw()

        deta=(ref_state[0].eta-state[0].eta)

        self.assertTrue(abs(deta.min().number) < 1.e-10)
        self.assertTrue(abs(deta.max().number) < 1.e-10)

        print deta.min(),deta.max()

        raw_input()
