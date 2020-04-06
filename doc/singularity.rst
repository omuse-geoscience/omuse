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

Launch the container to start Jupyter server inside::

    singularity run --contain -B examples:/opt/notebooks,run:/run/user omuse.img 
    
Then visit the reported localhost:8888 with a browser.

In the singularity command above, the --contain option disables mounting your home directory whereas the -B option specifies paths that will be mounted inside the container, in this case the examples distributed with OMUSE. If you have a folder with notebooks using OMUSE, you can mount this instead of the examples directory to execute them with the singularity container. The container above can also be used to run a regular python script that uses OMUSE functionality by piping the contents to the container python command::

    cat myscript.py | singularity exec --contain omuse.img python3 

Finally, it is also possible to launch a shell inside the container::

    singularity shell --contain omuse.img

to execute your python code with all OMUSE dependencies findable. We advise to use the --contain option whenever you
have OMUSE installed on your host system in $HOME/.local 
