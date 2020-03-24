module dflowfm_omuse

    character(len=256) :: input_configfile = "amuse.mdu"


contains

  function initialize() result(ret)
    integer :: ret
    ret=0
  end function

  function commit_parameters() result(ret)
    use dflowfm_omuse_lib
    integer :: ret
  
    ret=initialize_dflowfm(input_configfile)

    if(ret.NE.0) return
    
    ret=initialize_index_map()

  end function
  
  function get_model_time(t_) result(ret)
    use dflowfm_omuse_lib
    integer :: ret
    double precision :: t_
    
    ret=get_current_time(t_)
  
  end function
  
  function evolve_model(t_) result(ret)
    use dflowfm_omuse_lib, evolve_model_ => evolve_model
    integer :: ret
    double precision :: t_
    
    ret=evolve_model_(t_)
  
  end function

 ! Flow node numbering:
 ! 1:ndx2D, ndx2D+1:ndxi, ndxi+1:ndx1Db, ndx1Db:ndx
 ! ^ 2D int ^ 1D int      ^ 1D bnd       ^ 2D bnd ^ total

  function get_2d_flow_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=ndx2d
    ret=0
  end function

  function get_1d_flow_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=ndx2d+1
    imax=ndxi
    ret=0
  end function

  function get_1d_boundary_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=ndxi+1
    imax=ndx1Db
    ret=0
  end function

  function get_2d_boundary_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=ndx1Db+1
    imax=ndx
    ret=0
  end function

! the following global indexing assumes no 1D nodes

  function get_global_2d_flow_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=ndxi_glob
    ret=0
  end function

  function get_global_2d_boundary_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=ndxi_glob+1
    imax=ndx_glob
    ret=0
  end function


  function get_x_position(i, x, n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n,i_
    double precision :: x(n)
    ret=get_x_position_(i, x,n)
  end function

  function get_y_position(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_y_position_(i, x,n)
  end function

  function get_water_level(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_water_level_(i,x,n)
  end function

end module

