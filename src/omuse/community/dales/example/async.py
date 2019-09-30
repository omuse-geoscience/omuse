"""A demonstration of asynchronous interface functions with the Dales model

Asynchronous function calls and requests is a mechanism in OMUSE to
let multiple models or code instances do work concurrently.

Any method call on an OMUSE object can be made asynchronous, by
appending ``.asynchronous`` to the method's name.
When an asynchronous call is made to a code, the Python script
continues to run while the call is being made.  The Python script can
then perform computations or calls to other codes in OMUSE.
Subsequent calls to the same code instance can be made, but they
will be performed in sequence.


An asynchronous OMUSE function call returns  immediately,
i.e. it does not wait for the worker code to finish the function.
Instead of the normal return value, the asynchronous function returns a request
object, which can be used to access the result of the function call later.
The request can also be added to a request pool, which manages
several concurrent requests.

Caveats:

Mixing asynchronous and synchronous calls produces correct results,
but has consequences for the sequencing of the work: Performing a
normal, synchronous call to a code after one or more asynchronous
calls, causes an implicit wait for the asynchronous calls to complete.

Accessing result() on a request causes a wait for the request to
complete.

When making several calls to the same code instance, the first call
begins immediately. The second call begins only when the code is waited on,
not automatically when the first call completes.

"""

from __future__ import division
from __future__ import print_function

from omuse.community.dales.interface import Dales
from omuse.units import units
from amuse.rfi.async_request import AsyncRequestsPool
import time


# create Dales objects
d1 = Dales(workdir='dales1', channel_type='sockets', number_of_workers=1, case='bomex')
d2 = Dales(workdir='dales2', channel_type='sockets', number_of_workers=1, case='bomex')

# explicitly initialize the codes
# otherwise implicitly done when calling evolve_model
d1.commit_parameters()
d2.commit_parameters()

# add parameter redirection='none' to see DALES diagnostics output

target_time = 120 | units.s # target time

# create a pool for managing asynchronous requests
t = time.time()
pool = AsyncRequestsPool()

# add requests to the two codes to the pool
request1 = d1.evolve_model.asynchronous(target_time, exactEnd=True)
pool.add_request(request1)

request2 = d2.evolve_model.asynchronous(target_time, exactEnd=True)
pool.add_request(request2)

print ('Generating asynchronous requests  %f s'%(time.time()-t))

# wait for the requests to finish
print('Calling pool.waitall()')
t = time.time()
pool.waitall()
print('pool.waitall() returned %f s'%(time.time()-t))


# setting grid data

# normal synchronous call
t = time.time()
d1.fields[:,:,3:6].U = 1 | units.m / units.s
print ('Synchronous setting %f s'%(time.time()-t))

# asynchronous call
t = time.time()
d2.fields[:,:,3:6].request.U = 1 | units.m / units.s
print ('Asynchronous setting %f s'%(time.time()-t))

# getting grid data

t = time.time()
# synchronous call, returns an array
uprofile = d1.fields[5,5,1:10].U
print ('Synchronous getting %f s'%(time.time()-t))
print(uprofile)

# asynchronous call, returns a request.
# the request.result() gives the result of the call, an array.
t = time.time()
request = d2.fields[5,5,1:10].request.U
print ('Making asynchronous getting call %f s'%(time.time()-t))
t = time.time()
uprofile = request.result() # implicit wait for the request to finish
print ('Retreiving asynchronous result %f s'%(time.time()-t))
print(uprofile)
 

