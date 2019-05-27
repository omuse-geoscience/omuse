#!/usr/bin/env python

import os
import sys

#purpose of this script:
#
#calling:
#diff -u cdo-1.7.0rc5/src/Makefile src/src/Makefile > makefile.patch
# or
#patch cdo-1.7.0rc5/src/Makefile < makefile.patch
#
#for a bunch of files

code = "cdo-1.7.0rc5" 

source_dir = "src_original/src/"
target_dir = code + "/src/" 
patch_dir = "patches/"

files = ["Makefile", "remap_scrip_io.c", "remap_conserv_scrip.c", "EOFs.c" ]

def make_patches(filename):
    command = 'diff' + " -I [C/]" + ' -u ' + target_dir + filename + " " + source_dir + filename + ' > ' + patch_dir + filename + ".patch"
    print command
    os.system(command)

def apply_patches(filename):
    command =  "patch" + " " + target_dir + filename + " < " + patch_dir + filename + ".patch" 
    os.system(command)

def print_usage():
    print "Usage: make_patches.py [diff|patches]"
    print "       option diff generates patch files and stores them in directory patch"
    print "       option patch applies the patches to the code"

if __name__ == "__main__": 
    if len(sys.argv) != 2 or (len(sys.argv) == 2 and not (sys.argv[1] == "diff" or sys.argv[1] == "patch")):
        print_usage()
        exit()

    if len(sys.argv) > 1 and sys.argv[1] == "diff":
        if not os.path.exists(patch_dir):
            os.makedirs(patch_dir)
        for file in files:
            make_patches(file)

    if len(sys.argv) > 1 and sys.argv[1] == "patch":
        for file in files:
            apply_patches(file)




