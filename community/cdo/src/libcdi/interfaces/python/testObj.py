import CdiObj

ifile = "../testdata/mulval.grb"

cdi = CdiObj.Cdi(ifile)

print('Stream: ',cdi.streamID,' vlistID:',cdi.vlistID,' nvars:{d}', cdi.nvars)

print('#========== TAXES ====================================#')
for k in range(cdi.taxes.size()):
  print(k,": ", cdi.taxes[k].unitname)

print('#========== GRIDS ====================================#')
for k in range(cdi.grids.size()):
  print(k,": ", cdi.grids[k].size,' ', cdi.grids[k].xname,' ', cdi.grids[k].yname,' ', cdi.grids[k].ylongname) 

print("#========== ZAXES ====================================#")
for k in range(cdi.zaxes.size()):
  print(k,": ", cdi.zaxes[k].size,' ', cdi.zaxes[k].name,' ', cdi.zaxes[k].units)

print("#========== VARIABLES ================================#")
for k in range(cdi.variables.size()):
  v = cdi.variables[k]
  print(v.name," ",v.size, " ", v.missval)

print("#========== VARIABLEcdi.NAMES =================================#")
for k in range(cdi.variables.size()):
  print(cdi.variables[k].longname,' ',cdi.variables[k].units)

print("#========== VAR by index ======================================#")
var = cdi.variables[1]
var.getValues()
val = var.values
i=0; print('val[',i,'] = ',val[i])
i=1; print('val[',i,'] = ',val[i])
i=2; print('val[',i,'] = ',val[i])
i=3; print('val[',i,'] = ',val[i])
i=4; print('val[',i,'] = ',val[i])
i=5; print('val[',i,'] = ',val[i])
print("#========= Var by name ===============================#")
name ="tsurf"
newvar = cdi.var[name]
print("name ",name," var.name: ", newvar.name, " var.grids.xsize: " , newvar.grid.xsize)
print("#========= Var by code ===============================#")
code = 169
newvar = cdi.varByCode[code]
newvar.sinfo()
