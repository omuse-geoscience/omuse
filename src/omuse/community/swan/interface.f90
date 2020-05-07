module swan_interface
  use amuse_swan
  implicit none

  real :: begin_time=0.

  character*12 :: grid_type="regular"
  character*12 :: calc_mode="stationary"
  character*12 :: coordinates="cartesian"
  character*12 :: projection_method="quasi-cart."
  logical :: wrap_x = .FALSE.
  integer :: number_dimensions=2
  integer :: grid_mxc,grid_myc
  real :: grid_xpc,grid_ypc,grid_xlenc,grid_ylenc,grid_alpc
    
  character*12 :: input_grid_type="regular"
  real :: input_xp,input_yp,input_dx,input_dy,input_alp
  integer :: input_mx,input_my
  logical ::  use_input_depth=.TRUE., &
              use_input_water_level=.FALSE., &
              use_input_current=.FALSE., &
              use_input_air_sea_temp_diff=.FALSE., &
              use_input_friction=.FALSE., &
              use_input_wind=.FALSE., &
              use_input_plant_density=.FALSE., &
              use_input_turbulent_visc=.FALSE., &
              use_input_mud_layer=.FALSE.

  character*80 :: north_boundary_spec_file="none", &
                  south_boundary_spec_file="none", &
                  east_boundary_spec_file="none", &
                  west_boundary_spec_file="none", &
                  unstructured_boundary_spec_file="none"

  integer :: boundary_marker

  logical ::  use_gen3=.FALSE., &
              use_breaking=.FALSE., &
              use_friction=.FALSE., &
              use_triads=.FALSE.

  logical ::  use_csigma_cfl_limiter=.FALSE., &
              use_ctheta_cfl_limiter=.FALSE.

  real :: umin=1.
  real :: rho_air=1.28
  real :: ctb=0.01
  real :: under_relaxation_factor=0.0

contains

include "getter_setters.f90"

function initialize_code(coord_, mode_, grid_,input_grid_) result(ret)
  integer :: ret
  character*12 :: coord_, mode_, grid_,input_grid_
  
  grid_type=grid_
  input_grid_type=input_grid_
  calc_mode=mode_
  coordinates=coord_

  ret=swan_entry()
  if(ret.NE.0) then
    return
  endif
  
  ret=swan_init()
  if(ret.NE.0) then
    ret=ret-10
    return
  endif
  
  ret=swan_init_mode(calc_mode, number_dimensions)
  if(ret.NE.0) then
    ret=ret-100
    return
  endif

  ret=swan_init_coord(coordinates, projection_method)
  if(ret.NE.0) then
    ret=ret-1000
    return
  endif
  
  PROJID='AMUSE'
  PROJNR='1'
  
  if(grid_type.EQ."regular") OPTG=1
  if(grid_type.EQ."curvilinear") OPTG=3
  if(grid_type.EQ."unstructured") OPTG=5

  umin=PWIND(12)
  rho_air=PWIND(16)
  ctb=PTURBV(1)
  
end function

function commit_grid_positions() result(ret)
  integer :: ret
  ret=0
  if(grid_type.EQ."unstructured") ret=swan_init_unstructured_comp_grid()
  
end function

function commit_grids() result(ret)
  integer :: ret
  ret=0
  if(grid_type.EQ."regular") ret=swan_init_regular_comp_grid()
  if(grid_type.EQ."curvilinear") ret=swan_init_curvilinear_comp_grid()
end function


function initialize_grid() result(ret)
  integer :: ret
  
  if(grid_type.EQ."regular") then
    ret=swan_init_regular_grid(grid_mxc,grid_myc, &
      grid_xpc,grid_ypc,grid_xlenc,grid_ylenc,grid_alpc, wrap_x)
    if(ret.NE.0) return
  endif
  if(grid_type.EQ."curvilinear") then
    ret=swan_init_curvilinear_grid()
    if(ret.NE.0) return
  endif
  if(grid_type.EQ."unstructured") then
    ret=swan_init_unstructured_grid()
    if(ret.NE.0) return
  endif
  ret=swan_init_freq_grid()
  if(ret.NE.0) return
  ret=swan_report_grids()  

end function

function initialize_input_grids() result(ret)
  integer :: ret
  ret=0
  
  if(input_grid_type.EQ."curvilinear") then
    ret=-2
    return
  endif
  
  if(use_input_depth) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(1, input_mx, input_my, input_xp,  &
                                      input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(1)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(DEPTH)) DEALLOCATE(DEPTH)
    ALLOCATE(DEPTH(MXG(1)*MYG(1)))
    DEPTH=1
    LEDS(1)=2
  endif

  if(use_input_current) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(2 , input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
      ret=swan_init_regular_input_grid(3 , input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(2)
      ret=swan_init_unstructured_input_grid(3)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(UXB)) DEALLOCATE(UXB)
    IF (ALLOCATED(UYB)) DEALLOCATE(UYB)
    ALLOCATE(UXB(MXG(2)*MYG(2)))      
    ALLOCATE(UYB(MXG(3)*MYG(3)))      
    UXB=0
    UYB=0
    LEDS(2)=2
    LEDS(3)=2
  endif
  
  if(use_input_friction) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(4, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(4)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(FRIC)) DEALLOCATE(FRIC)
    ALLOCATE(FRIC(MXG(4)*MYG(4)))
    FRIC=0
    LEDS(4)=2
  endif
  
  if(use_input_wind) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(5, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
      ret=swan_init_regular_input_grid(6, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(5)
      ret=swan_init_unstructured_input_grid(6)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(WXI)) DEALLOCATE(WXI)
    IF (ALLOCATED(WYI)) DEALLOCATE(WYI)
    ALLOCATE(WXI(MXG(5)*MYG(5)))
    ALLOCATE(WYI(MXG(6)*MYG(6)))
    WXI=0
    WYI=0
    LEDS(5)=2
    LEDS(6)=2
  endif

  if(use_input_water_level) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(7, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(7)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(WLEVL)) DEALLOCATE(WLEVL)
    ALLOCATE(WLEVL(MXG(7)*MYG(7)))
    WLEVL=0
    LEDS(7)=2
  endif
  
  if(use_input_air_sea_temp_diff) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(10, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(10)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(ASTDF)) DEALLOCATE(ASTDF)
    ALLOCATE(ASTDF(MXG(10)*MYG(10)))
    ASTDF=0
    LEDS(10)=2
  endif      
    
  if(use_input_plant_density) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(11, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(11)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(NPLAF)) DEALLOCATE(NPLAF)
    ALLOCATE(NPLAF(MXG(11)*MYG(11)))
    NPLAF=0
    LEDS(11)=2
  endif
  
  if(use_input_turbulent_visc) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(12, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(12)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(TURBF)) DEALLOCATE(TURBF)
    ALLOCATE(TURBF(MXG(12)*MYG(12)))
    TURBF=0
    LEDS(12)=2
  endif
  
  if(use_input_mud_layer) then
    if(input_grid_type.EQ."regular") then
      ret=swan_init_regular_input_grid(13, input_mx, input_my, input_xp,  &
                                        input_yp, input_alp, input_dx, input_dy)
    else if(input_grid_type.EQ."unstructured") then
      ret=swan_init_unstructured_input_grid(13)
    endif
    if(ret.NE.0) return
    IF (ALLOCATED(MUDLF)) DEALLOCATE(MUDLF)
    ALLOCATE(MUDLF(MXG(13)*MYG(13)))
    MUDLF=0
    LEDS(13)=2
  endif    
  
end function


function commit_parameters() result(ret)
  integer :: ret
  ret=0
  TIMCO=begin_time

  if(.not.use_input_depth) then
    ret=-1
    return
  endif
  if(use_input_current) ICUR=1
  if(use_input_air_sea_temp_diff) then
    VARAST = .TRUE.
    IF (JASTD2.LE.1) THEN
      MCMVAR = MCMVAR + 2
      JASTD2 = MCMVAR - 1
      JASTD3 = MCMVAR
      ALOCMP = .TRUE.
    ENDIF
  endif
  if(use_input_friction) then
    VARFR  = .TRUE.
    MCMVAR = MCMVAR + 2
    JFRC2  = MCMVAR - 1
    JFRC3  = MCMVAR
    ALOCMP = .TRUE.
  endif
  if(use_input_wind) then
    IWIND  = 3 ! note: fix
    VARWI  = .TRUE.
  endif
  if(use_input_water_level) VARWLV=.TRUE.
  if(use_input_plant_density) then
    VARNPL = .TRUE.
    IF (JNPLA2.LE.1) THEN
      MCMVAR = MCMVAR + 2
      JNPLA2 = MCMVAR - 1
      JNPLA3 = MCMVAR
      ALOCMP = .TRUE.
    ENDIF
  endif
  if(use_input_turbulent_visc) then
    ITURBV = 1
    VARTUR = .TRUE.
    PTURBV(2) = -1
    IF (JTURB2.LE.1) THEN
      MCMVAR = MCMVAR + 2
      JTURB2 = MCMVAR - 1
      JTURB3 = MCMVAR
      ALOCMP = .TRUE.
    ENDIF
  endif
  if(use_input_mud_layer) then
    VARMUD = .TRUE.
    IMUD   = 1
    IF (JMUDL2.LE.1) THEN
      MCMVAR = MCMVAR + 3
      JMUDL1 = MCMVAR - 2
      JMUDL2 = MCMVAR - 1
      JMUDL3 = MCMVAR
      ALOCMP = .TRUE.
    ENDIF
  endif
  
  if(.not.use_input_wind) then
    IWIND=3
    VARWI=.FALSE.
  endif
  
  if(use_gen3) then
    ret=swan_physics_gen3()
    if(ret.NE.0) return
  endif
  if(use_breaking) then
    ret=swan_physics_breaking()
    if(ret.NE.0) return
  endif
  if(use_friction) then
    ret=swan_physics_friction()
    if(ret.NE.0) return
  endif
  if(use_triads) then
    ret=swan_physics_triads()
    if(ret.NE.0) return
  endif
  
  if(use_csigma_cfl_limiter) then 
    PNUMS(33)=1.
  else
    PNUMS(33)=0
  endif
  if(use_ctheta_cfl_limiter) then
    PNUMS(35)=1.
  else
    PNUMS(35)=0
  endif

  PWIND(12)=umin
  PWIND(16)=rho_air
  PWIND(17)=RHO
  PWIND(9)=rho_air/RHO

  PTURBV(1)=ctb
  
  PNUMS(30)=under_relaxation_factor

  ret=0
end function

function initialize_boundary() result(ret)
  integer :: ret
  CHARACTER(len=255) :: cwd
  ret=0

  call getcwd(cwd)

  if(grid_type.EQ."regular") then
    if(north_boundary_spec_file.NE."none") then
      ret=swan_regular_add_boundary_from_file("N",north_boundary_spec_file)
      if(ret.NE.0) return
    endif
    if(south_boundary_spec_file.NE."none") then
      ret=swan_regular_add_boundary_from_file("S",south_boundary_spec_file)
      if(ret.NE.0) return
    endif
    if(west_boundary_spec_file.NE."none") then
      ret=swan_regular_add_boundary_from_file("W",west_boundary_spec_file)
      if(ret.NE.0) return
    endif
    if(east_boundary_spec_file.NE."none") then
      ret=swan_regular_add_boundary_from_file("S",east_boundary_spec_file)
      if(ret.NE.0) return
    endif    
  endif

  if(grid_type.EQ."unstructured") then
    ret=swan_init_unstructured_boundary()
    if(ret.NE.0) return
    if(unstructured_boundary_spec_file.NE."none") then
      ret=swan_unstructured_add_boundary_from_file(boundary_marker,unstructured_boundary_spec_file)
      if(ret.NE.0) return
    endif
!    ret=swan_compute_wave_induced_force_unstructured() ! maybe move later on
    if(ret.NE.0) return
  endif

  ret=0
end function

function evolve_model(tend) result(ret)
  integer :: ret
  real*8 :: tend
  real :: dt_save
  character*4 :: comp

  ret=0
  dt_save=DT
  if(calc_mode.EQ."stationary") then
    TINIC=TIMCO
    TFINC=tend
    TIMCO=tend
    DT = 1.E10
    RDTIM = 0.
    NSTATC = 0
    MTC = 1
    comp="COMP"
  else
    TINIC=TIMCO
    if(DT.LE.0.) then
      ret=-1
      return
    endif
    MTC = CEILING ((tend - TINIC)/DT)
    if(MTC.LE.0) then
      MTC=0
      TFINC=TINIC
      RDTIM=0.
      comp="NOCO"
    else
      TFINC=TINIC+MTC*DT
      RDTIM=1./DT
      NSTATC=1
    endif
  endif
  IF (NSTATC.EQ.0) THEN
    ITERMX = MXITST
  ELSE
    ITERMX = MXITNS
  ENDIF
  ret=swan_compute("COMP")
  DT=dt_save
  
  if(grid_type.EQ."unstructured") ret=swan_compute_wave_induced_force_unstructured()
  
end function


function recommit_parameters() result(ret)
  integer :: ret
! as yet nothing tbd
  ret=0
end function
function cleanup_code() result(ret)
  integer :: ret
  ret=swan_cleanup()
end function

function get_depth_regular(i,j,x,n) result(ret)
  integer :: ret,i(n),j(n),m,n,ii
  real*8 :: x(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.MXC.OR. &
        j(m).LT.1.OR.j(m).GT.MYC) then
      x(m)=0
      ret=-1
      cycle
    endif
    ii=KGRPNT(i(m),j(m))
    if(ii.EQ.1) then
      x(m)=0
    else if(ii.LE.1.OR.ii.GT.MCGRD) then
      x(m)=0
      ret=-1
    else
      x(m)=COMPDA(ii,JDP2)
    endif
  enddo
end function

function get_wave_stress_regular(i,j,x,y,n) result(ret)
  integer :: ret,i(n),j(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=-2
end function

function get_ac2_regular(i,j,k,l,x,n) result(ret)
  integer :: ret,i(n),j(n),k(n),l(n),m,n,ii
  real*8 :: x(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.MXC.OR. &
        j(m).LT.1.OR.j(m).GT.MYC.OR. &
        k(m).LT.1.OR.k(m).GT.MDC.OR. &
        l(m).LT.1.OR.l(m).GT.MSC ) then
      x(m)=0
      ret=-1
      cycle
    endif
    ii=KGRPNT(i(m),j(m))
    if(ii.EQ.1) then
      x(m)=0
    else if(ii.LE.1.OR.ii.GT.MCGRD) then
      x(m)=0
      ret=-1
    else
      x(m)=AC2(k(m),l(m),ii)
    endif
  enddo
end function

function get_ac2_unstructured(ii,k,l,x,n) result(ret)
  integer :: ret,k(n),l(n),m,n,ii(n)
  real*8 :: x(n)
  ret=0
  do m=1,n
    if( ii(m).LT.1.OR.ii(m).GT.nvertsg.OR. &
        k(m).LT.1.OR.k(m).GT.MDC.OR. &
        l(m).LT.1.OR.l(m).GT.MSC ) then
      x(m)=0
      ret=-1
      cycle
    endif
    x(m)=AC2(k(m),l(m),ii(m))
  enddo
end function

function get_depth_unstructured(ii,x,n) result(ret)
  integer :: ret,m,n,ii(n)
  real*8 :: x(n)
  ret=0
  do m=1,n
    if( ii(m).LT.1.OR.ii(m).GT.nvertsg) then
      x(m)=0
      ret=-1
      cycle
    endif
    x(m)=COMPDA(ii(m),JDP2)
  enddo
end function

function get_wave_stress_unstructured(ii,x,y,n) result(ret)
  integer :: ret,m,n,ii(n)
  real*8 :: x(n),y(n)
  ret=0
  do m=1,n
    if( ii(m).LT.1.OR.ii(m).GT.nvertsg) then
      x(m)=0
      y(m)=0
      ret=-1
      cycle
    endif
    x(m)=FORCE(ii(m),1)
    y(m)=FORCE(ii(m),2)
  enddo
end function

function get_grid_position_regular(i,j,x,y,n) result(ret)
  integer :: ret,i(n),j(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.MXC.OR. &
        j(m).LT.1.OR.j(m).GT.MYC ) then
      x(m)=0
      y(m)=0
      ret=-1
      cycle
    endif
    x(m)=XCGRID(i(m),j(m))
    y(m)=YCGRID(i(m),j(m))
  enddo
end function  

function get_grid_lonlat_regular(i,j,x,y,n) result(ret)
  integer :: ret,i(n),j(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=get_grid_position_regular(i,j,x,y,n)
end function


function set_grid_position_unstructured(i,x,y,n) result(ret)
  integer :: ret,i(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.nvertsg ) then
      ret=-1
      cycle
    endif
    xcugrd(i(m))=x(m)
    ycugrd(i(m))=y(m)
  enddo
  if(LXOFFS) ret=-2
end function  

function set_grid_lonlat_unstructured(i,x,y,n) result(ret)
  integer :: ret,i(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=set_grid_position_unstructured(i,x,y,n)
end function

function get_grid_position_unstructured(i,x,y,n) result(ret)
  integer :: ret,i(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.nvertsg ) then
      x(m)=0
      y(m)=0
      ret=-1
      cycle
    endif
    x(m)=xcugrd(i(m))
    y(m)=ycugrd(i(m))
  enddo
  if(LXOFFS) then
    x=x+XOFFS
    y=y+YOFFS
  endif
end function

function get_grid_lonlat_unstructured(i,x,y,n) result(ret)
  integer :: ret,i(n),m,n,ii
  real*8 :: x(n),y(n)
  ret=get_grid_position_unstructured(i,x,y,n)
end function

function set_grid_vmark_unstructured(i,vm,n) result(ret)
  integer :: ret,i(n),m,n,ii,vm(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.nvertsg ) then
      ret=-1
      cycle
    endif
    vmark(i(m))=vm(m)
  enddo
end function  

function get_grid_vmark_unstructured(i,vm,n) result(ret)
  integer :: ret,i(n),m,n,ii,vm(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.nvertsg ) then
      vm(m)=0
      ret=-1
      cycle
    endif
    vm(m)=vmark(i(m))
  enddo
end function

function get_input_grid_position_regular(i,j, x,y,n) result(ret)
    integer :: ret,n,i(n),j(n),k,ii, igrid=1
    real*8 :: x(n),y(n),dx,dy
    ret=0
    do k=1,n
      if(i(k).LT.1.OR.i(k).GT.MXG(igrid).OR. &
         j(k).LT.1.OR.j(k).GT.MYG(igrid) ) THEN
        ret=-1
      else
        dx=DXG(igrid)*(i(k)-1)
        dy=DYG(igrid)*(j(k)-1)
        x(k)=XPG(igrid)+(COSPG(igrid)*dx-SINPG(igrid)*dy)
        y(k)=YPG(igrid)+(SINPG(igrid)*dx+COSPG(igrid)*dy)
      endif
    enddo
end function

function get_input_grid_lonlat_regular(i,j, x,y,n) result(ret)
    integer :: ret,n,i(n),j(n),k,ii, igrid=1
    real*8 :: x(n),y(n),dx,dy
    ret=get_input_grid_position_regular(i,j,x,y,n)
end function

! note the offset -1  
function set_element_nodes(i,n1,n2,n3,n) result(ret)
  integer :: ret,i(n),m,n,n1(n),n2(n),n3(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.ncellsg ) then
      ret=-1
      cycle
    endif
    kvertc(1,i(m))=n1(m)+1
    kvertc(2,i(m))=n2(m)+1
    kvertc(3,i(m))=n3(m)+1
  enddo
end function

function get_element_nodes(i,n1,n2,n3,n) result(ret)
  integer :: ret,i(n),m,n,n1(n),n2(n),n3(n)
  ret=0
  do m=1,n
    if( i(m).LT.1.OR.i(m).GT.ncellsg ) then
      n1(m)=0
      n2(m)=0
      n3(m)=0
      ret=-1
      cycle
    endif
    n1(m)=kvertc(1,i(m))-1
    n2(m)=kvertc(2,i(m))-1
    n3(m)=kvertc(3,i(m))-1
  enddo
end function

function get_number_of_unstructured_boundary_segments(n) result(ret)
  integer :: n,ret
  n=nbpol
  ret=0
end function

function get_number_of_nodes_in_unstructured_boundary_segment(i,n) result(ret)
  integer :: i,n,ret
  if(i.LT.1.OR.i.GT.nbpol) then
    ret=-1
    return
  endif
  n=nbpt(i)
  ret=0
end function

function get_unstructured_boundary_node(i,j,n) result(ret)
  integer :: i,j,n,ret
  if(j.LT.1.OR.j.GT.nbpol) then
    ret=-1
    return
  endif
  if(i.LT.1.OR.i.GT.nbpt(j)) then
    ret=-2
    return
  endif
  n=blist(i,j)
  ret=0
end function

function get_time(x) result(ret)
  integer :: ret
  real*8 :: x
  x=TIMCO
  ret=0
end function

function set_rearth(x) result(ret)
  integer :: ret
  real*8 :: x
  REARTH=x
  LENDEG = REARTH * PI / 180.
  ret=0
end function
function get_rearth(x) result(ret)
  integer :: ret
  real*8 :: x
  x=REARTH
  ret=0
end function

function set_wdip(x) result(ret)
  integer :: ret
  real*8 :: x
  real :: ALTMP, DEGCNV
  WDIP=x
  WDIP = DEGCNV(WDIP)!
  ALTMP = WDIP / 360.
  WDIP = PI2 * (ALTMP - NINT(ALTMP))
  ret=0
end function
function get_wdip(x) result(ret)
  integer :: ret
  real*8 :: x
  real :: DEGCNV
  x=DEGCNV(WDIP/PI2*360.)
  ret=0
end function

! note real vs double
function get_exc_value(i,x) result(ret)
  integer :: ret,i
  real*8 :: x
  x=EXCFLD(i)
  ret=0  
end function
! note real vs double
function set_exc_value(i,x) result(ret)
  integer :: ret,i
  real*8 :: x
  EXCFLD(i)=x
  ret=0  
end function

end module
