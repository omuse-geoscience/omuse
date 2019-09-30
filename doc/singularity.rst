.. _Singularity-section:

Singularity image
=================

A Singularity recipe is included with OMUSE. When building it, the result is a Singularity image
which is a portable way to test OMUSE and the included models.

The image includes Jupyter for interactive Python notebooks and the following OMUSE models:

 * cdo
 * dales
 * qgcm (hogg2006,ocean_only,hogg2006_20km)
 * qgmodel
 * swan

   

Building the image::

    sudo singularity build omuse.img Singularity 

Create a directory which the container can use at runtime for storing notebooks::

    mkdir run

Launch the container - start Jupyter server inside::

    singularity run -B src/omuse/community/dales/example/:/opt/notebooks/,run:/run/user/ ./omuse.img 
    
Then visit localhost:8888 with a browser.
In the singularity command, the -B option specifies paths that will be mounted inside the container.
This example shows the path of the dales examples folder, which contains the notebook `bubble-notebook.ipynb`
demonstrating a warm air bubble simulation.

It is also possible to launch a shell inside the container::

    singularity shell ./omuse.img

Inside the container, one can then run the example programs::
  
    cd src/omuse/community/dales/example/
    python bubble.py
    

    
