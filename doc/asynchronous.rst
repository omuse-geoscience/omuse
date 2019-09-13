.. _asynchronous:

Asynchronous calls
==================

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
::

   request1 = d1.evolve_model.asynchronous(target_time, exactEnd=True)
   # ... do something else ...
   print(request1.result()) 

Example
-------

Running two models simultaneously::

    from omuse.community.dales.interface import Dales
    from omuse.units import units
    from amuse.rfi.async_request import AsyncRequestsPool

    # create Dales objects
    d1 = Dales(workdir='dales1', channel_type='sockets', number_of_workers=1, case='bomex')
    d2 = Dales(workdir='dales2', channel_type='sockets', number_of_workers=1, case='bomex')

    # create a pool for managing asynchronous requests
    pool = AsyncRequestsPool()

    # add requests to the two codes to the pool
    request1 = d1.evolve_model.asynchronous(target_time, exactEnd=True)
    pool.add_request(request1)

    request2 = d2.evolve_model.asynchronous(target_time, exactEnd=True)
    pool.add_request(request2)

    # wait for the requests to finish
    pool.waitall()


Asynchronous variable access::
  
    # setting grid data
    # normal synchronous call
    d1.fields[:,:,3:6].U = 1 | units.m / units.s

    # asynchronous call
    d2.fields[:,:,3:6].request.U = 1 | units.m / units.s

    # getting grid data
    # synchronous call, returns an array
    uprofile = d1.fields[5,5,1:10].U

    # asynchronous call, returns a request.
    request = d2.fields[5,5,1:10].request.U
    uprofile = request.result() # retrieve result. Implicit wait for the request to finish

 

Caveats
-------

* Mixing asynchronous and synchronous calls produces correct results,
  but has consequences for the sequencing of the work: Performing a
  normal, synchronous call to a code after one or more asynchronous
  calls, causes an implicit wait for the asynchronous calls to complete.

* Accessing ``result()`` on a request causes a wait for the request to
  complete.

* When making several calls to the same code instance, the first call
  begins immediately. The second call begins only when the code is waited on,
  not automatically when the first call completes.
