module dflowfm_omuse

    character(len=256) :: input_configfile = "omuse.mdu"

    logical :: use_wind=.FALSE., use_patm=.FALSE. 
    logical :: use_waterlevel=.FALSE. 

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

    if(ret.NE.0) return

    ret=initialize_link_map()

    if(ret.NE.0) return

    ret=initialize_interface_forcings( use_wind, use_patm)
    
    if(ret.NE.0) return

    ret=initialize_interface_boundaries(use_waterlevel)    
        
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

  function get_use_wind(x) result(ret)
    integer :: ret
    logical :: x
    x=use_wind
    ret=0
  end function

  function set_use_wind(x) result(ret)
    integer :: ret
    logical :: x
    use_wind=x
    ret=0
  end function

  function get_use_patm(x) result(ret)
    integer :: ret
    logical :: x
    x=use_patm
    ret=0
  end function

  function set_use_patm(x) result(ret)
    integer :: ret
    logical :: x
    use_patm=x
    ret=0
  end function

  function get_use_waterlevel(x) result(ret)
    integer :: ret
    logical :: x
    x=use_waterlevel
    ret=0
  end function

  function set_use_waterlevel(x) result(ret)
    integer :: ret
    logical :: x
    use_waterlevel=x
    ret=0
  end function


 ! Flow node numbering:
 ! 1:ndx2D, ndx2D+1:ndxi, ndxi+1:ndx1Db, ndx1Db:ndx
 ! ^ 2D int ^ 1D int      ^ 1D bnd       ^ 2D bnd ^ total

  function get_flow_nodes_range(imin,imax) result(ret)
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

  function get_boundary_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=ndx1Db+1
    imax=ndx
    ret=0
  end function

! the following global indexing assumes no 1D nodes

  function get_global_flow_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=ndx_glob
    ret=0
  end function

  function get_global_internal_flow_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=ndxi_glob
    ret=0
  end function

  function get_global_boundary_flow_nodes_range(imin,imax) result(ret)
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

  function get_ucx(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_ucx_(i,x,n)
  end function

  function get_ucy(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_ucy_(i,x,n)
  end function

  function get_patm(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_patm_(i,x,n)
  end function

  function set_patm(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=set_patm_(i,x,n)
  end function



! for net nodes:

  function get_net_nodes_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=numk
    ret=0
  end function

  function get_x_position_net_nodes(i, x, n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n,i_
    double precision :: x(n)
    do i_=1, n
      x(i_)=xk(i(i_))
    enddo    
    ret=0
  end function

  function get_y_position_net_nodes(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    do i_=1, n
      x(i_)=yk(i(i_))
    enddo    
    ret=0
  end function

! for flow links

  function get_global_flow_links_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=lnx_glob
    ret=0
  end function
  
  function get_global_internal_flow_links_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=1
    imax=lnxi_glob
    ret=0
  end function  

  function get_global_boundary_flow_links_range(imin,imax) result(ret)
    use dflowfm_omuse_lib
    integer :: ret,imin,imax
    imin=lnxi_glob+1
    imax=lnx_glob
    ret=0
  end function

  function get_x_position_flow_links(i, x, n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n,i_
    double precision :: x(n)
    ret=get_x_position_flow_links_(i, x,n)
    ret=0
  end function

  function get_y_position_flow_links(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_y_position_flow_links_(i, x,n)
    ret=0
  end function

  function get_wx(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_wx_(i, x,n)
    ret=0
  end function

  function get_wy(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_wy_(i, x,n)
    ret=0
  end function

  function set_wx(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=set_wx_(i, x,n)
    ret=0
  end function

  function set_wy(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=set_wy_(i, x,n)
    ret=0
  end function

! waterlevel boundaries

  function get_zbndz(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_zbndz_(i, x,n)
    ret=0
  end function

  function set_zbndz(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=set_zbndz_(i, x,n)
    ret=0
  end function

  function get_is_waterlevel_bnd(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    logical :: x(n)
    integer :: x_(n)
    ret=get_is_waterlevel_bnd_(i, x_,n)
    where(x_.EQ.1) x=.TRUE.
    ret=0
  end function

  function get_xbndz(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_xbndz_(i, x,n)
    ret=0
  end function

  function get_ybndz(i, x,n) result (ret)
    use dflowfm_omuse_lib
    integer :: ret,i(n),n
    double precision :: x(n)
    ret=get_ybndz_(i, x,n)
    ret=0
  end function


end module

