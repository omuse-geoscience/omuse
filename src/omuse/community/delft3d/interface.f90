module dflowfm_omuse


contains

  function initialize() result(ret)
    use dflowfm_omuse_lib
    integer :: ret
    character(len=256) :: somefile = "westerscheldt.mdu"
  
    ret=initialize_dflowfm(somefile)
  
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


end module

