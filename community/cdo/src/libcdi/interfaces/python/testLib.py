import CdiLib
ifile = "../testdata/mulval.nc"
streamID = CdiLib.streamOpenRead(ifile)
vlistID  = CdiLib.streamInqVlist(streamID)
nvars    = CdiLib.vlistNvars(vlistID)
for i in range(0,nvars):
  print CdiLib.vlistInqVarCode(vlistID, i),

CdiLib.streamClose(streamID)
