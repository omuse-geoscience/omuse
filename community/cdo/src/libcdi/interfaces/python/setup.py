import os, string

INCFLAGS = []
CFLAGS   = []
LDFLAGS  = []
LIBS     = []
if 'INCFLAGS' in os.environ:
   INCFLAGS = os.environ['INCFLAGS'].split()

if 'CFLAGS' in os.environ:
   CFLAGS  = os.environ['CFLAGS'].split()

if 'LDFLAGS' in os.environ:
   LDFLAGS  = os.environ['LDFLAGS'].split('-L')

if 'LIBS' in os.environ:
   LIBS     = os.environ['LIBS'].split()
