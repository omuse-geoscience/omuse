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

    singularity run -B notebooks:/opt/notebooks/,run:/run/user/ ./omuse.img 
    
Then visit localhost:8888 with a browser.
