# ptype:
# simple: rw, scalar value and implementation generated by interface
# normal: rw, scalar value, custom implementation 
# ro: read only, scalar value, generated 
# vector: rw, generated vector valued

class FortranCodeGenerator(object):
    _getter_string="""
      function get_{0}(x) result(ret)
        integer :: ret
        {1} :: x
        x={0}
        ret=0
      end function
                   """
    _setter_string="""
      function set_{0}(x) result(ret)
      integer :: ret
        {1} :: x
        {0}=x
        ret=0
      end function
                   """
    
    _vector_getter_string="""
      function get_{0}(i,x) result(ret)
        integer :: i,ret
        {1} :: x
        x={0}(i)
        ret=0
      end function
                   """
    _vector_setter_string="""
      function set_{0}(i,x) result(ret)
      integer :: i,ret
        {1} :: x
        {0}(i)=x
        ret=0
      end function
                   """

    _grid_getter_template="""
      function get_{0}({2},{0}_out_,n) result(ret)
      integer :: n,{3},k,ret
      {1} :: {0}_out_(n)
      ret=0
      do k=1,n
        {5}
        {0}_out_(k)={0}({4})
      enddo
      end function
      """

    _grid_setter_template="""
      function set_{0}({2},{0}_in_,n) result(ret)
      integer :: n,{3},k,ret
      {1} :: {0}_in_(n)
      ret=0
      do k=1,n
        {5}
        {0}({4})={0}_in_(k)
      enddo
      end function
      """    
    
    datatypedict={"string" : "character(len=*) ", 
                  "float64" : "real*8",
                  "float32" : "real",
                  "int32" : "integer",
                  "bool" : "logical",
                  }

    def __init__(self,parameter_definition=None, grid_variable_definition=None):
        self.parameter_definition=parameter_definition
        self.grid_variable_definition=grid_variable_definition

    def _grid_format_arg(self,name, dtype,ndim=1, index_ranges=None):
      arg=','.join(['i{0}'.format(i) for i in range(ndim)]) 
      dec=','.join(['i{0}(n)'.format(i) for i in range(ndim)])
      assign=','.join(['i{0}(k)'.format(i) for i in range(ndim)])
      if index_ranges is None:
        check="! no check on index range"
      else:
        checks=(".OR. &\n"+" "*11).join(["({0}.LT.{1}.OR.{0}.GT.{2})".format('i{0}(k)'.format(i), 
                              index_ranges[i][0],index_ranges[i][1]) for i in range(ndim)]) 
        check="if({0}) then\n".format(checks)+ \
              " "*10+"ret=-1\n"+ \
              " "*10+"exit\n"+ \
              " "*8+"endif"              
      dtype=self.datatypedict[dtype]
      return name,dtype,arg,dec,assign,check
    
    def grid_getter(self, name, dtype,ndim=1, index_ranges=None):
      return self._grid_getter_template.format(*self._grid_format_arg(name,dtype,ndim,index_ranges))

    def grid_setter(self, name, dtype,ndim=1, index_ranges=None):
      return self._grid_setter_template.format(*self._grid_format_arg(name,dtype,ndim,index_ranges))

    def parameter_getter_setters(self):
        filestring=""
        py_to_f=self.datatypedict
        for par,d in self.parameter_definition.items():
          if d["ptype"] in ["ro"]:
            filestring+=self._getter_string.format(d["short"],py_to_f[d["dtype"]])
          if d["ptype"] in ["simple"]:
            filestring+=self._setter_string.format(d["short"],py_to_f[d["dtype"]])
            filestring+=self._getter_string.format(d["short"],py_to_f[d["dtype"]])
          if d["ptype"] in ["vector"]:
            filestring+=self._vector_setter_string.format(d["short"],py_to_f[d["dtype"]])
            filestring+=self._vector_getter_string.format(d["short"],py_to_f[d["dtype"]])
    
        return filestring

    def grid_getter_setters(self):
        string=""
        for var, d in self.grid_variable_definition.items():
            vartype=d.get("vartype", None)
            dtype=d.get("dtype", "float64")
            ndim=d.get("ndim", 1)
            index_ranges=d.get("index_ranges",None)
            for forvar in d["forvar"]:
              string+=self.grid_getter(forvar, dtype, ndim, index_ranges)
              if vartype!="ro":
                string+=self.grid_setter(forvar, dtype, ndim, index_ranges)
        return string

    def generate_getters_setters(self,filename="getter_setters.f90"):
        filestring=""
        #~ filestring+=input_grid_string(_unstructured_input_grid_template)
        #~ filestring+=input_grid_string(_regular_input_grid_template)
        filestring+=self.parameter_getter_setters()
        filestring+=self.grid_getter_setters()
        with open(filename,"w") as f:
            f.write(filestring)

    def generate_parameter_interface_functions(self):
        output=""
        for par,d in self.parameter_definition.items():
            dtype=d["dtype"]
            if hasattr(d["default"],"unit"):
              unit=d["default"].unit.reference_string()
            else:
              unit="None"
            short=d["short"]
            ptype=d["ptype"]
            if ptype in ["ro"]:
              output+=("@legacy_function\ndef get_"+short+"():\n  function = LegacyFunctionSpecification()\n"
                   "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
                   "  function.result_type = 'int32'\n  return function\n")
            if ptype in ["simple"]:
              output+=("@legacy_function\ndef get_"+short+"():\n  function = LegacyFunctionSpecification()\n"
                   "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
                   "  function.result_type = 'int32'\n  return function\n")
              output+=("@legacy_function\ndef set_"+short+"():\n  function = LegacyFunctionSpecification()\n"
                   "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n"
                   "  function.result_type = 'int32'\n  return function\n")
            if ptype in ["vector"]:
              output+=("@legacy_function\ndef get_"+short+"():\n  function = LegacyFunctionSpecification()\n"
                   "  function.addParameter('i', dtype='i', direction=function.IN)\n"
                   "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n"
                   "  function.can_handle_array=True\n"
                   "  function.result_type = 'int32'\n  return function\n")
              output+=("@legacy_function\ndef set_"+short+"():\n  function = LegacyFunctionSpecification()\n"
                   "  function.addParameter('i', dtype='i', direction=function.IN)\n"
                   "  function.addParameter('"+short+"', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n"
                   "  function.can_handle_array=True\n"
                   "  function.result_type = 'int32'\n  return function\n")
              length=d["length"]
              output+=( "def get_"+short+"_range(self):\n" + (
                    ("  return 1," + str(length)) if isinstance(length, int) else
                    ("  return 1, self.get_"+length+"()['"+length+"']\n")  )
                    )
        return output

    def generate_grid_interface_functions(self):
        output=""
        for var, d in self.grid_variable_definition.items():
            vartype=d.get("vartype", None)
            dtype=d.get("dtype", "float64")
            dtype=dtype.__name__ if isinstance(dtype, type) else str(dtype)
            ndim=d.get("ndim", 1)
            index_ranges=d.get("index_ranges",None)
            unit=d.get("unit",None)
            unit="None" if unit is None else unit.reference_string() 
            for pyvar, forvar in zip(d["pyvar"], d["forvar"]):
              if vartype!="ro":
                output+=("@legacy_function\ndef set_"+forvar+"():\n  function = LegacyFunctionSpecification()\n" + \
                     "".join(["  function.addParameter('index{0}', dtype='i', direction=function.IN)\n".format(i) for i in range(ndim)]) + \
                     "  function.addParameter('"+pyvar+"', dtype='"+dtype+"', direction=function.IN, unit="+unit+")\n" + \
                     "  function.addParameter('n', direction=function.LENGTH)\n" + \
                     "  function.must_handle_array=True\n" + \
                     "  function.result_type = 'int32'\n  return function\n")
              output+=("@legacy_function\ndef get_"+forvar+"():\n  function = LegacyFunctionSpecification()\n" + \
                   "".join(["  function.addParameter('index{0}', dtype='i', direction=function.IN)\n".format(i) for i in range(ndim)]) + \
                   "  function.addParameter('"+pyvar+"', dtype='"+dtype+"', direction=function.OUT, unit="+unit+")\n" + \
                   "  function.addParameter('n', direction=function.LENGTH)\n" + \
                   "  function.must_handle_array=True\n" + \
                   "  function.result_type = 'int32'\n  return function\n")
        return output

    def generate_interface_functions(self):
        output=""
        output+=self.generate_parameter_interface_functions()
        output+=self.generate_grid_interface_functions()
        return output

    def generate_parameter_definitions(self, object):      
        for name,d in self.parameter_definition.items():
            short=d["short"]
            ptype=d["ptype"]
            dtype=d["dtype"]
            getter="get_"+short
            if ptype in ["simple","normal","vector"]:
              setter="set_"+short
            else:
              setter=None
            range_method="get_"+short+"_range"
            if ptype in ["simple", "normal", "ro"]:
              if dtype!='bool':
                  object.add_method_parameter(
                      getter,
                      setter,
                      name,
                      d["description"], 
                      d["default"]
                  )
              else:
                  object.add_boolean_parameter(
                      getter,
                      setter,
                      name,
                      d["description"], 
                      d["default"]
                  )
            else:
              object.add_array_parameter(
                      getter,
                      setter,
                      range_method, 
                      name,
                      d["description"]
                  )

if __name__=="__main__":
  grid_var={
    "pressure" : dict( pyvar=["pressure"], forvar=["pom"], dtype="float64", ndim=1, index_ranges=[(1,10)],vartype="ro"),
    "test2" : dict( pyvar=["y"], forvar=["y"], dtype="float64", ndim=2, index_ranges=[(1,"nla"),(1,"nla")])
  }
  f=FortranCodeGenerator(grid_variable_definition=grid_var)
  print(f.grid_getter_setters())
  print(f.generate_grid_interface_functions())
  
