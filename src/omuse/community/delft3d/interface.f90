module dflowfm_omuse
!~   use hash

! global node counts

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

  function get_x_position(i, x, n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    x(1:n)=xz(i(1:n))
    ret=0
  end function

  function get_y_position(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    x(1:n)=yz(i(1:n))
    ret=0
  end function

  function get_water_level(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    x(1:n)=s1(i(1:n))
    ret=0
  end function

end module

