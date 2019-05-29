module pop_interface

! !DESCRIPTION:
! This is the interface to the main driver for the Parallel Ocean Program (POP).
!

! !USES:
  use POP_KindsMod
  use POP_ErrorMod
  use POP_FieldMod
  use POP_GridHorzMod
  use POP_CommMod
  use POP_InitMod
  use POP_HaloMod

  use kinds_mod, only: int_kind, r8
  use blocks
  use communicate, only: my_task, master_task, MPI_COMM_OCN
  use gather_scatter
  use grid
  use constants, only: grav
  use domain, only: distrb_clinic, nprocs_clinic, nprocs_tropic, clinic_distribution_type, &
         tropic_distribution_type, ew_boundary_type, ns_boundary_type
  use forcing_fields, only: SMF, SMFT, lsmft_avail, STF
  use forcing_shf, only: set_shf, SHF_QSW, shf_filename, shf_data_type, shf_interp_freq, shf_interp_type
  use forcing_ws, only: set_ws, ws_filename, ws_data_type, ws_data_next, ws_data_update, ws_interp_freq, ws_interp_type, &
      ws_interp_next, ws_interp_last, ws_interp_inc
  use forcing_sfwf, only: sfwf_filename, sfwf_data_type, sfwf_interp_freq, sfwf_interp_type, fwf_imposed
  use forcing_tools, only: never
  use initial, only: init_ts_option, init_ts_file, init_ts_file_fmt
  use io_types, only: nml_filename
  use operators, only: wcalc
  use prognostic, only: TRACER, PSURF, UVEL, VVEL, UBTROP, VBTROP, RHO, curtime
  use restart, only: restart_freq_opt, restart_freq, restart_outfile
  use tavg, only: tavg_freq_opt, tavg_freq, tavg_outfile
  use movie, only: movie_freq_opt, movie_freq, movie_outfile
  use timers, only: timer_print_all, get_timer, timer_start, timer_stop
  use time_management
!   use time_management, only: get_time_flag_id, check_time_flag, &
!       nsteps_run, stdout, exit_pop, override_time_flag
  use step_mod, only: step
  use diagnostics, only: check_KE
  use output, only: output_driver
  use exit_mod, only: sigAbort, sigExit
! use io, only:


  implicit none

!interface
!    subroutine start_timer ( ) bind (c)
!      use iso_c_binding
!    end subroutine start_timer
!    subroutine stop_timer ( time ) bind (c)
!      use iso_c_binding
!      real (c_float) :: time
!    end subroutine stop_timer
!end interface
  



!-----------------------------------------------------------------------
!
! module private variables
!
!-----------------------------------------------------------------------

  integer (POP_i4) :: &
     errorCode ! error flag

  integer (int_kind) :: &
     timer_step,        &! timer number for step
     fstop_now           ! flag id for stop_now flag

  logical :: &
     fupdate_coriolis,  &! flag for update coriolis force
     fupdate_wind_stress ! flag for update wind stress

  logical :: initialized = .false.

  real (r8), dimension (nx_block, ny_block, max_blocks_clinic) :: &
     tau_x, tau_y   ! wind stress on t grid in N/m**2

  real (r8), dimension (nx_block, ny_block, km, max_blocks_clinic) :: &
     WVEL           ! 3D vertical velocity field, derived from UVEL and VVEL using wcalc

  real (r8), dimension (:,:), allocatable :: & 
     WORK_G            ! temporary 2D array for gathered data


contains


!-----------------------------------------------------------------------
!
! pre-initialize the model run
!
!-----------------------------------------------------------------------
function initialize_code() result(ret)
  integer :: ret

  errorCode = POP_Success

  call POP_Initialize0(errorCode)

  if (errorCode /= POP_Success) then
    ret=-1
  else
    ret=0
  endif
end function


!-----------------------------------------------------------------------
!
! fully initialize the model run
!
!-----------------------------------------------------------------------
function commit_parameters() result(ret)
  integer :: ret

  !call POP_CommInitMessageEnvironment

  ret=0
  errorCode = POP_Success

  call POP_Initialize(errorCode)

  fstop_now = get_time_flag_id('stop_now')
  fupdate_coriolis = .false.
  fupdate_wind_stress = .false.

  !allocate space for gather of global field
  if (my_task == master_task) then
    allocate(WORK_G(nx_global,ny_global))
  endif

  !-----------------------------------------------------------------------
  !
  ! start up the main timer
  !
  !-----------------------------------------------------------------------

  call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)
  call timer_start(timer_total)

  call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)

  !-----------------------------------------------------------------------
  !
  ! set windstress to that it can be read before the first step is run
  !
  !-----------------------------------------------------------------------
  if (lsmft_avail) then
    call set_ws(SMF, SMFT=SMFT)
  else
    call set_ws(SMF)
  endif

  !-----------------------------------------------------------------------
  !
  ! set surface heat flux so that it can be read before the first step
  ! this also sets SHF_QSW, although not passed as an argument
  !
  !-----------------------------------------------------------------------
  call set_shf(STF)


  initialized = .true.

  ! ensure the forcings can be read
  ret = prepare_parameters()

  ! ensure full grid info is present on master_task
  call initialize_global_grid

  if (errorCode /= POP_Success) then
    ret=-1
  endif
end function




!-----------------------------------------------------------------------
!
! advance the model in time until the model has run for 'tend' seconds
!
!-----------------------------------------------------------------------
function evolve_model(tend) result(ret)
  integer :: ret
  real(8), intent(in) :: tend

  ! stepsize_next is the time in seconds for next timestep 
  ! it is set and updated by POP's time management, this way
  ! (instead of using dtt) half steps are correctly supported
  do while ( tsecond < tend - stepsize_next*0.5 .and. &
             errorCode == POP_Success )

    call timer_start(timer_step)
    call step(errorCode)
    call timer_stop(timer_step)

    !***
    !*** exit if energy is blowing
    !***
    if (check_KE(100.0_r8)) then
      call override_time_flag(fstop_now,value=.true.)
      call output_driver(errorCode)
      call exit_POP(sigAbort,'ERROR: k.e. > 100 ')
    endif

    !-----------------------------------------------------------------------
    !
    ! write restart dumps and archiving
    !
    !-----------------------------------------------------------------------
    call output_driver(errorCode)

  enddo

  if (errorCode == POP_Success) then
    ret=0
  else
    ret=-1
  endif
end function


!-----------------------------------------------------------------------
!
! print timing information and clean up various environments if
! they have been used
!
!-----------------------------------------------------------------------
function cleanup_code() result(ret)
  integer :: ret

  !override time flag fstop_now and produce final output
  call override_time_flag(fstop_now,value=.true.)
  call output_driver(errorCode)

  !produce final output
  call timer_stop(timer_total)
  call timer_print_all(stats=.true.)

  call POP_ErrorPrint(errorCode, printTask=POP_masterTask)

  !call exit_POP(sigExit,'Successful completion of POP run')
  if (my_task == master_task) then
    write(6,*) 'Successful completion of POP run'  
  endif

  ret=0
end function




!-----------------------------------------------------------------------
!
! Returns the elapsed time in the model in seconds
!
!-----------------------------------------------------------------------
function get_model_time(tout) result(ret)
  integer :: ret
  real(8), intent(out) :: tout

  tout = tsecond 

  ret=0
end function


!-----------------------------------------------------------------------
!
! Returns the default number of seconds in a full time step
!
!-----------------------------------------------------------------------
function get_timestep(dtout) result(ret)
  integer :: ret
  real(8), intent(out) :: dtout

  dtout = dtt

  ret=0
end function

!-----------------------------------------------------------------------
!
! Returns number of seconds in the next step, can be less than default for half steps
!
!-----------------------------------------------------------------------
function get_timestep_next(dtout) result(ret)
  integer :: ret
  real(8), intent(out) :: dtout

  dtout = stepsize_next

  ret=0
end function





!-----------------------------------------------------------------------
!
! Getter and setter for wind stress per grid point
! Be sure to call prepare_wind_stress() before calling this
!
!-----------------------------------------------------------------------
function get_node_wind_stress(g_i, g_j, tau_x_, tau_y_, n) result(ret)
  integer :: ret,n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: tau_x_, tau_y_

  call get_gridded_variable_vector(g_i, g_j, tau_x, tau_x_, n)
  call get_gridded_variable_vector(g_i, g_j, tau_y, tau_y_, n)

  ret=0
end function

!----------------------------------------------------------------------
!
! This function converts SMF or SMFT to wind stress 
! The result is stored in the internal variables tau_x and tau_y
! which can in turn be accessed by get_node_wind_stress()
!
!----------------------------------------------------------------------
function prepare_wind_stress() result(ret)
  integer :: ret

  tau_x =  SMF(:,:,1,:)/momentum_factor * cos(ANGLE) - SMF(:,:,2,:)/momentum_factor * sin(ANGLE)
  tau_y =  ( SMF(:,:,2,:)/momentum_factor + tau_x * sin(ANGLE) ) / cos (ANGLE)

  if (errorCode == POP_Success) then
    ret=0
  else
    ret=-1
  endif
end function
function set_node_wind_stress (g_i, g_j, tau_x_, tau_y_, n) result(ret)
  integer :: ret, n
  real*8, dimension(n), intent(in) :: tau_x_, tau_y_
  integer, dimension(n), intent(in) :: g_i, g_j ! global i and j indexes for the gridpoint

  call set_gridded_variable_vector(g_i, g_j, tau_x, tau_x_, n)
  call set_gridded_variable_vector(g_i, g_j, tau_y, tau_y_, n)

  fupdate_wind_stress = .true.

  ret=0
end function

function set_wind_stress() result(ret)
  integer :: ret

  fupdate_wind_stress = .false.

  !if we want to overwrite the wind stress, then we have to make sure that SMF is
  !not overwritten by some later interpolation routine, therefore, we have to
  !enforce the following settings
  ws_data_type = 'none'
  ws_data_next = never
  ws_data_update = never
  ws_interp_freq = 'never'
  ws_interp_next = never
  ws_interp_last = never
  ws_interp_inc  = c0

  ! rotate_wind_stress applies tau_x and tau_y to the U-grid
  ! variable SMF (surface momentum fluxes (wind stress))
  ! wind stress in N/m^2 is converted to vel flux (cm^2/s^2)
  ! the method expects tau_x and tau_y to be supplied for the U grid

  ! set_gridded_variable only sets grid points in the physical domain
  ! still need to update halo values
  call POP_HaloUpdate(tau_x(:,:,:), POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)

  call POP_HaloUpdate(tau_y(:,:,:), POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)

  call rotate_wind_stress(tau_x, tau_y)

  if (errorCode == POP_Success) then
    ret=0
  else
    ret=-1
  endif
end function


! rotate_wind_stress()
!   rotates true zonal/meridional wind stress into local
!   coordinates, converts to dyne/cm**2
!
! adapted from forcing_coupled.F90, but modified to:
! 1. not apply the land mask RCALC
! 2. to have tau_x and tau_y be supplied on the U-grid instead of the T-grid
subroutine rotate_wind_stress(WORK1, WORK2)

  real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) ::   &
      WORK1, WORK2        ! contains tau_x and tau_y from coupler on U-grid

  integer :: errorCode

  !-----------------------------------------------------------------------
  !
  !  rotate and convert
  !
  !-----------------------------------------------------------------------

   SMF(:,:,1,:) = (WORK1(:,:,:)*cos(ANGLE(:,:,:)) +  &
                   WORK2(:,:,:)*sin(ANGLE(:,:,:)))*  &
                   momentum_factor
   SMF(:,:,2,:) = (WORK2(:,:,:)*cos(ANGLE(:,:,:)) -  &
                   WORK1(:,:,:)*sin(ANGLE(:,:,:)))*  &
                   momentum_factor

  !-----------------------------------------------------------------------
  !
  !  perform halo updates following the vector rotation
  !
  !-----------------------------------------------------------------------

  call POP_HaloUpdate(SMF(:,:,1,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_POP_r8)

  call POP_HaloUpdate(SMF(:,:,2,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_POP_r8)

end subroutine rotate_wind_stress





!-----------------------------------------------------------------------
!
! Generic getter and setters for gridded variables
!
!-----------------------------------------------------------------------
subroutine get_gridded_variable(g_i, g_j, grid, value)
  integer, intent(in) :: g_i, g_j
  real*8, intent(in), dimension(:,:,:) :: grid
  real*8, intent(out) :: value

  integer :: i,j,iblock
  type (block) :: this_block

  value = 0.0 !important for MPI_SUM used later

  ! search rather than iterate over all points
  do iblock=1, nblocks_clinic
    this_block = get_block(blocks_clinic(iblock),iblock)
  
    ! i_glob and j_glob denote the global index of gridpoint within a block
    ! the ghost cells are also given indices but they may have very different values on cyclic borders etc
    ! for now, only look at the grid points that are not ghost cells
    !
    ! nx_block and ny_block are the dimensions of the block
    ! nghost+1 is the first index within the block that is part of the physical domain
    ! nx_block-nghost is the last index within the block that is part of the physical domain
    if (this_block%i_glob(nghost+1) <= g_i .and. this_block%i_glob(nx_block-nghost) >= g_i .and. &
        this_block%j_glob(nghost+1) <= g_j .and. this_block%j_glob(ny_block-nghost) >= g_j ) then 

      !find index in this block from global index (1+ because indexes start at 1 in Fortran)
      i = 1+nghost + g_i - this_block%i_glob(1+nghost) 
      j = 1+nghost + g_j - this_block%j_glob(1+nghost)
      
      !debugging: check if this went ok
      !if (this_block%i_glob(i) /= g_i .or. this_block%j_glob(j) /= g_j) then
      !  write (*,*) 'Error: g_i =', g_i, ', g_j=', g_j, ' and found i=', this_block%i_glob(i), ' j=', this_block%j_glob(j) 
      !endif

      value = grid(i,j,iblock)

    endif

  enddo

end subroutine get_gridded_variable

subroutine set_gridded_variable(g_i, g_j, grid, value)
  integer, intent(in) :: g_i, g_j
  real*8, intent(inout), dimension(:,:,:) :: grid
  real*8, intent(in) :: value

  integer :: i,j,iblock
  type (block) :: this_block

  ! search rather than iterate over all points
  do iblock=1, nblocks_clinic
    this_block = get_block(blocks_clinic(iblock),iblock)
  
    if (this_block%i_glob(nghost+1) <= g_i .and. this_block%i_glob(nx_block-nghost) >= g_i .and. &
        this_block%j_glob(nghost+1) <= g_j .and. this_block%j_glob(ny_block-nghost) >= g_j ) then 

      i = 1+nghost + g_i - this_block%i_glob(1+nghost) 
      j = 1+nghost + g_j - this_block%j_glob(1+nghost)

      !debugging: check if this went ok
      !if (this_block%i_glob(i) /= g_i .or. this_block%j_glob(j) /= g_j) then
      !  write (*,*) 'Error: g_i =', g_i, ', g_j=', g_j, ' and found i=', this_block%i_glob(i), ' j=', this_block%j_glob(j) 
      !endif

      grid(i,j,iblock) = value
      
    endif

  enddo

end subroutine set_gridded_variable


subroutine get_gridded_variable_3D(g_i, g_j, k, grid, value)
  integer, intent(in) :: g_i, g_j, k
  real*8, intent(in), dimension(:,:,:,:) :: grid
  real*8, intent(out) :: value

  integer :: i,j,iblock
  type (block) :: this_block

  do iblock=1, nblocks_clinic
    this_block = get_block(blocks_clinic(iblock),iblock)
    if (this_block%i_glob(nghost+1) <= g_i .and. this_block%i_glob(nx_block-nghost) >= g_i .and. &
        this_block%j_glob(nghost+1) <= g_j .and. this_block%j_glob(ny_block-nghost) >= g_j ) then 

      i = 1+nghost + g_i - this_block%i_glob(1+nghost) 
      j = 1+nghost + g_j - this_block%j_glob(1+nghost)
      
      value = grid(i,j,k,iblock)

    endif
  enddo
end subroutine get_gridded_variable_3D
subroutine set_gridded_variable_3D(g_i, g_j, k, grid, value)
  integer, intent(in) :: g_i, g_j, k
  real*8, intent(inout), dimension(:,:,:,:) :: grid
  real*8, intent(in) :: value

  integer :: i,j,iblock
  type (block) :: this_block

  do iblock=1, nblocks_clinic
    this_block = get_block(blocks_clinic(iblock),iblock)
    if (this_block%i_glob(nghost+1) <= g_i .and. this_block%i_glob(nx_block-nghost) >= g_i .and. &
        this_block%j_glob(nghost+1) <= g_j .and. this_block%j_glob(ny_block-nghost) >= g_j ) then 

      i = 1+nghost + g_i - this_block%i_glob(1+nghost) 
      j = 1+nghost + g_j - this_block%j_glob(1+nghost)

      grid(i,j,k,iblock) = value
    endif
  enddo
  
end subroutine set_gridded_variable_3D



subroutine get_gather(g_i, g_j, grid, value_, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: value_
  real*8, dimension(:,:,:), intent(inout) :: grid

  integer :: ii
  value_ = 0.0 !first zero the array

  call gather_global(WORK_G, grid, master_task, distrb_clinic)
  if (my_task == master_task) then
    do ii=1,n
      value_(ii) = WORK_G(g_i(ii),g_j(ii))
    enddo
  endif
end subroutine get_gather

subroutine get_gather_3D(g_i, g_j, k, grid, value_, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j, k
  real*8, dimension(n), intent(out) :: value_
  real*8, dimension(:,:,:,:), intent(inout) :: grid

  integer :: ii, ki

  do ki=1,km
    if (ANY(k(:) == ki)) then
      call gather_global(WORK_G, grid(:,:,ki,:), master_task, distrb_clinic)
      if (my_task == master_task) then
        do ii=1,n
          if (k(ii) == ki) then
            value_(ii) = WORK_G(g_i(ii),g_j(ii))
          endif
        enddo
      endif
    endif
  enddo

end subroutine get_gather_3D


subroutine get_gridded_variable_vector(g_i, g_j, grid, value, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, intent(in), dimension(:,:,:) :: grid
  real*8, dimension(n), intent(out) :: value

  integer :: ii,ierr

  include 'mpif.h'

  value = 0.0 !important for MPI_SUM used later

  do ii=1,n
    call get_gridded_variable(g_i(ii), g_j(ii), grid, value(ii))
  enddo

  !only the output of the node with MPI rank 0 is significant in AMUSE, reduce to single value per element
  if(my_task == master_task) then
    call MPI_REDUCE(MPI_IN_PLACE, value, n, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  else
    call MPI_REDUCE(value, 0, n, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  endif

end subroutine get_gridded_variable_vector

subroutine set_gridded_variable_vector(g_i, g_j, grid, value, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, intent(inout), dimension(:,:,:) :: grid
  real*8, dimension(n), intent(in) :: value

  integer :: ii

  do ii=1,n
    call set_gridded_variable(g_i(ii), g_j(ii), grid, value(ii))
  enddo

end subroutine set_gridded_variable_vector

subroutine get_gridded_variable_vector_3D(g_i, g_j, k, grid, value, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j, k
  real*8, intent(in), dimension(:,:,:,:) :: grid
  real*8, dimension(n), intent(out) :: value

  integer :: ii,ierr

  include 'mpif.h'

  value = 0.0 !important for MPI_SUM used later

  do ii=1,n
    call get_gridded_variable_3D(g_i(ii), g_j(ii), k(ii), grid, value(ii))
  enddo

  !only the output of the node with MPI rank 0 is significant in AMUSE, reduce to single value per element
  if(my_task == master_task) then
    call MPI_REDUCE(MPI_IN_PLACE, value, n, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  else
    call MPI_REDUCE(value, 0, n, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  endif
end subroutine get_gridded_variable_vector_3D

subroutine set_gridded_variable_vector_3D(g_i, g_j, k, grid, value, n)
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j, k
  real*8, intent(inout), dimension(:,:,:,:) :: grid
  real*8, dimension(n), intent(in) :: value
  integer :: ii

  do ii=1,n
    call set_gridded_variable_3D(g_i(ii), g_j(ii), k(ii), grid, value(ii))
  enddo

end subroutine set_gridded_variable_vector_3D




!-----------------------------------------------------------------------
!
! Getter and setter for Coriolis Force
!
! these are currently implemented using the FCOR on the U grid
! set_coriolis_f() can be used to update the FCORT on the T grid based
! on the modified FCOR for the U grid
!
!-----------------------------------------------------------------------
function get_node_coriolis_f(g_i, g_j, corif_, n) result(ret)
  integer :: ret, n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: corif_

  if (initialized .eqv. .false.) then
    call exit_POP(sigAbort, 'Error: get_node_coriolis_f() called before initialize_code()')
  endif

  call get_gridded_variable_vector(g_i, g_j, FCOR, corif_, n)

  ret=0
end function
function set_node_coriolis_f(g_i, g_j, corif_, n) result(ret)
  integer :: ret, n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(in) :: corif_

  call set_gridded_variable_vector(g_i, g_j, FCOR, corif_, n)

  fupdate_coriolis = .true.

  ret=0
end function

! set_coriolis_f updates FCORT (T-grid) based on FCOR (U-grid)
function set_coriolis_f() result(ret)
  integer :: ret
  integer :: iblock

  fupdate_coriolis = .false.

  !first update ghost cells for FCOR on U-grid
  call POP_HaloUpdate(FCOR, POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)

  !convert to T-grid
  do iblock=1, nblocks_clinic
    call ugrid_to_tgrid(FCORT, FCOR, iblock)
  enddo

  !update halo cells on T-grid
  call POP_HaloUpdate(FCORT, POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)

  if (errorCode == POP_Success) then
    ret=0
  else
    ret=-1
  endif
end function


!-----------------------------------------------------------------------
!
! Recommit forcings will call all functions that had forcings set to
! temporary locations and will ensure these are properly converted to
! internal variables used by POP
!
!-----------------------------------------------------------------------
function recommit_forcings() result(ret)
  integer :: ret
  ret = 0

  if (fupdate_wind_stress) then
    ret = set_wind_stress()
  endif

  if (ret /= 0) return
  
  if (fupdate_coriolis) then
    ret = set_coriolis_f()
  endif
end function


!-----------------------------------------------------------------------
!
!  Getters for node and element surface state
!
!-----------------------------------------------------------------------
function get_node_barotropic_vel(g_i, g_j, uvel_, vvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: uvel_, vvel_

  if (n < nx_global*ny_global) then
    call get_gridded_variable_vector(g_i, g_j, UBTROP(:,:,curtime,:), uvel_, n)
    call get_gridded_variable_vector(g_i, g_j, VBTROP(:,:,curtime,:), vvel_, n)
  else
    call get_gather(g_i, g_j, UBTROP(:,:,curtime,:), uvel_, n)
    call get_gather(g_i, g_j, VBTROP(:,:,curtime,:), vvel_, n)
  endif

  ret=0
end function
function get_node_surface_state(g_i, g_j, ssh_, uvel_, vvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: ssh_, uvel_, vvel_

  if (n < nx_global*ny_global) then
    call get_gridded_variable_vector(g_i, g_j, PSURF(:,:,curtime,:), ssh_, n)
    ssh_ = ssh_ / grav

    call get_gridded_variable_vector(g_i, g_j, UVEL(:,:,1,curtime,:), uvel_, n)
    call get_gridded_variable_vector(g_i, g_j, VVEL(:,:,1,curtime,:), vvel_, n)
  else
    call get_gather(g_i, g_j, PSURF(:,:,curtime,:), ssh_, n)
    ssh_ = ssh_ / grav
    call get_gather(g_i, g_j, UVEL(:,:,1,curtime,:), uvel_, n)
    call get_gather(g_i, g_j, VVEL(:,:,1,curtime,:), vvel_, n)
  endif

  ret=0
end function
function get_element_surface_state(g_i, g_j, temp_, salt_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: temp_, salt_

  if (my_task == master_task) then
    write (*,*) 'get_element_surface_state() called with n=', n
  endif

  if (n < nx_global*ny_global) then
    call get_gridded_variable_vector(g_i, g_j, TRACER(:,:,1,1,curtime,:), temp_, n)
    call get_gridded_variable_vector(g_i, g_j, TRACER(:,:,1,2,curtime,:), salt_, n)
  else
    call get_gather(g_i, g_j, TRACER(:,:,1,1,curtime,:), temp_, n)
    call get_gather(g_i, g_j, TRACER(:,:,1,2,curtime,:), salt_, n)
  endif

  ret=0
end function

function get_element_surface_heat_flux(g_i, g_j, shf_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: shf_

  real*8, dimension(nx_block, ny_block, max_blocks_clinic) :: WORK1
  integer :: iblock

  do iblock=1, nblocks_clinic
    where (KMT(:,:,iblock) > 0)
      WORK1(:,:,iblock) = (STF(:,:,1,iblock)+SHF_QSW(:,:,iblock))/hflux_factor ! W/m^2
    elsewhere
      WORK1(:,:,iblock) = c0
    end where
  end do

  call get_gridded_variable_vector(g_i, g_j, WORK1, shf_, n)

  ret=0
end function


!-----------------------------------------------------------------------
!
! Getters for grid info
!
! Nodes are on the U grid
! Elements on the T grid
!
!-----------------------------------------------------------------------
function get_node_position(g_i, g_j, lat_, lon_, n) result(ret)
  integer :: ret,n,ii
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: lat_, lon_

!  real :: time
  integer :: ierr

!  if (my_task == master_task) then
!    write (*,*) 'get_node_position() called n=', n
!  endif

!  if (n < nx_global*ny_global) then
!    call get_gridded_variable_vector(g_i, g_j, ULAT, lat_, n)
!    call get_gridded_variable_vector(g_i, g_j, ULON, lon_, n)
!  else

!  time = 0.0
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call start_timer
!   call get_gather(g_i, g_j, ULAT, lat_, n)
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call stop_timer(time)

!  if (my_task == master_task) then
!    write(*,*) 'get_gather took: ', time, 'ms .'
!  endif

!  time = 0.0
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call start_timer
!  call get_gridded_variable_vector(g_i, g_j, ULAT, lat_, n)
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call stop_timer(time)

!  if (my_task == master_task) then
!    write(*,*) 'get_gridded_variable_vector took: ', time, ' ms.'
!  endif
!   call get_gather(g_i, g_j, ULON, lon_, n)

  if (my_task == master_task) then
    do ii=1,n
      lon_(ii) = ULON_G(g_i(ii),g_j(ii))
      lat_(ii) = ULAT_G(g_i(ii),g_j(ii))
    enddo
  endif

!  endif

  ret=0
end function

function get_element_position(g_i, g_j, lat_, lon_, n) result(ret)
  integer :: ret,n
  integer :: ii
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: lat_, lon_

  !if (n < nx_global*ny_global) then
  !  call get_gridded_variable_vector(g_i, g_j, TLAT, lat_, n)
  !  call get_gridded_variable_vector(g_i, g_j, TLON, lon_, n)
  !else
  ! call get_gather(g_i, g_j, TLAT, lat_, n)
  ! call get_gather(g_i, g_j, TLON, lon_, n)
  !endif
  if (my_task == master_task) then
    do ii=1,n
      lon_(ii) = TLON_G(g_i(ii),g_j(ii))
      lat_(ii) = TLAT_G(g_i(ii),g_j(ii))
    enddo
  endif

  ret=0
end function

function get_node_vposition(g_i, g_j, k, depth_, n) result(ret)
  integer :: ret,n
  integer, dimension(n), intent(in) :: g_i, g_j, k
  real*8, dimension(n), intent(out) :: depth_

  ret=0
  if (partial_bottom_cells) then
    call get_gridded_variable_vector_3D(g_i, g_j, k, DZU, depth_, n)
  else
    ret = get_zt(k, depth_, n)
  endif

end function

function get_element_vposition(g_i, g_j, k, depth_, n) result(ret)
  integer :: ret,n
  integer, dimension(n), intent(in) :: g_i, g_j, k
  real*8, dimension(n), intent(out) :: depth_

  ret=0
  if (partial_bottom_cells) then
    call get_gridded_variable_vector_3D(g_i, g_j, k, DZT, depth_, n)
  else
    ret = get_zt(k, depth_, n)
  endif

end function

function get_number_of_nodes(num) result(ret)
  integer :: ret
  integer, intent(out) :: num
  num = nx_global * ny_global
  ret=0
end function

function get_domain_size(x_, y_) result(ret)
  integer :: ret
  integer, intent(out) :: x_, y_
  x_ = nx_global
  y_ = ny_global
  ret=0
end function

function get_number_of_vertical_levels(km_) result(ret)
  integer :: ret
  integer, intent(out) :: km_
  km_ = km
  ret=0
end function

function get_node_depth(g_i, g_j, depth_, n) result(ret)
  integer :: ret,n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: depth_

  if (n < nx_global*ny_global) then
    call get_gridded_variable_vector(g_i, g_j, HU, depth_, n)
  else
    call get_gather(g_i, g_j, HU, depth_, n)
  endif

  ret=0
end function

function get_element_depth(g_i, g_j, depth_, n) result(ret)
  integer :: ret,n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: depth_

  if (n < nx_global*ny_global) then
    call get_gridded_variable_vector(g_i, g_j, HT, depth_, n)
  else
    call get_gather(g_i, g_j, HT, depth_, n)
  endif

  ret=0
end function

function get_zt(k, zt_, n) result(ret)
  integer :: ret,n,i
  integer, dimension(n), intent(in) :: k
  real*8, dimension(n), intent(out) :: zt_

  zt_(1:n)=zt(k(1:n))

  ret=0
end function





!-----------------------------------------------------------------------
!
! Getters and setters for individual 3D fields
! this code was written optimistically and probably needs adjustment later
!
!-----------------------------------------------------------------------
function get_element3d_temperature(i, j, k, temp_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: temp_

!  real*4 :: time
  integer :: ierr

  if (my_task == master_task) then
    write (*,*) 'get_element3d_temperature() called with n=', n
  endif

!  time = 0.0
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call start_timer
!  call get_gather_3D(i, j, k, TRACER(:,:,:,1,curtime,:), temp_, n)

!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call stop_timer(time)

!  if (my_task == master_task) then
!    write(*,*) 'get_gather_3D took: ', time, 'ms .'
!  endif

!  time = 0.0
!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call start_timer
  call get_gridded_variable_vector_3D(i, j, k, TRACER(:,:,:,1,curtime,:), temp_, n)

!  call MPI_Barrier(MPI_COMM_OCN, ierr)
!  call stop_timer(time)

!  if (my_task == master_task) then
!    write(*,*) 'get_gridded_variable_vector_3D took: ', time, ' ms.'
!  endif

  ret=0
end function
function set_element3d_temperature(i, j, k, temp_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(in) :: temp_

  call set_gridded_variable_vector_3D(i, j, k, TRACER(:,:,:,1,curtime,:), temp_, n)

  ret=0
end function
function get_element3d_salinity(i, j, k, salt_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: salt_

  call get_gridded_variable_vector_3D(i, j, k, TRACER(:,:,:,2,curtime,:), salt_, n)
!  call get_gather_3D(i, j, k, TRACER(:,:,:,2,curtime,:), salt_, n)

  ret=0
end function
function set_element3d_salinity(i, j, k, salt_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(in) :: salt_

  call set_gridded_variable_vector_3D(i, j, k, TRACER(:,:,:,2,curtime,:), salt_, n)

  ret=0
end function

function get_node3d_velocity_xvel(i, j, k, uvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: uvel_

  call get_gridded_variable_vector_3D(i, j, k, UVEL(:,:,:,curtime,:), uvel_, n)
!  call get_gather_3D(i, j, k, UVEL(:,:,:,curtime,:), uvel_, n)

  ret=0
end function
function set_node3d_velocity_xvel(i, j, k, uvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(in) :: uvel_

  call set_gridded_variable_vector_3D(i, j, k, UVEL(:,:,:,curtime,:), uvel_, n)

  ret=0
end function
function get_node3d_velocity_yvel(i, j, k, vvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: vvel_

  call get_gridded_variable_vector_3D(i, j, k, VVEL(:,:,:,curtime,:), vvel_, n)
!  call get_gather_3D(i, j, k, VVEL(:,:,:,curtime,:), vvel_, n)

  ret=0
end function
function set_node3d_velocity_yvel(i, j, k, vvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(in) :: vvel_

  call set_gridded_variable_vector_3D(i, j, k, VVEL(:,:,:,curtime,:), vvel_, n)

  ret=0
end function

function get_node3d_velocity_zvel(i, j, k, wvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: wvel_

  call get_gridded_variable_vector_3D(i, j, k, WVEL(:,:,:,:), wvel_, n)

  ret=0
end function

function get_element3d_density(i, j, k, rho_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(out) :: rho_

  call get_gridded_variable_vector_3D(i, j, k, RHO(:,:,:,curtime,:), rho_, n)
!  call get_gather_3D(i, j, k, RHO(:,:,:,curtime,:), rho_, n)

  ret=0
end function
function set_element3d_density(i, j, k, rho_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: i, j, k
  real*8, dimension(n), intent(in) :: rho_

  call set_gridded_variable_vector_3D(i, j, k, RHO(:,:,:,curtime,:), rho_, n)

  ret=0
end function










!-----------------------------------------------------------------------
!
! State transition functions
!
! Currently these both only call one internal routine, but this is likely
! to be extended in the future
!
!-----------------------------------------------------------------------
function recommit_parameters() result (ret)
  integer :: ret
  ret=0

  ret = recommit_forcings()

end function
function prepare_parameters() result (ret)
  integer :: ret
  ret=0

  if (initialized .eqv. .false.) then
    call exit_POP(sigAbort, 'Error: prepare_parameters() called before initialize_code()')
  endif

  ret = prepare_wind_stress()

  ret = prepare_vertical_velocity()

end function


function prepare_vertical_velocity() result (ret)
  integer :: ret

   real (r8), allocatable, dimension(:,:,:) :: &
      WVELT   ! vertical velocity at T-points

   type (block) :: &
     this_block  ! block info for current block

   integer(int_kind) :: iblock, k

   allocate (WVELT(nx_block,ny_block,km))

   WVELT = c0

   do iblock = 1, nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      call wcalc(WVELT, UVEL(:,:,:,curtime,iblock),   &
                        VVEL(:,:,:,curtime,iblock), this_block)

      do k = 1, km
         call tgrid_to_ugrid(WVEL(:,:,k,iblock), WVELT(:,:,k), iblock)
      enddo

   enddo

   deallocate (WVELT)

  ret = 0
end function






!-----------------------------------------------------------------------
!
! Getters and settings for parameters
!
!-----------------------------------------------------------------------
function get_horiz_grid_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = horiz_grid_opt

  ret=0
end function
function set_horiz_grid_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'file' .OR. option == 'internal') then
    horiz_grid_opt = option
  else 
    ret=-1
  endif
end function
function get_horiz_grid_file(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = horiz_grid_file

  ret=0
end function
function set_horiz_grid_file(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  horiz_grid_opt = 'file'
  horiz_grid_file = option
end function


function get_vert_grid_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = vert_grid_opt

  ret=0
end function
function set_vert_grid_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'file' .OR. option == 'internal') then
    vert_grid_opt = option
  else 
    ret=-1
  endif
end function
function get_vert_grid_file(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = vert_grid_file

  ret=0
end function
function set_vert_grid_file(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  vert_grid_opt = 'file'
  vert_grid_file = option
end function

function get_topography_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = topography_opt

  ret=0
end function
function set_topography_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'file' .OR. option == 'internal') then
    topography_opt = option
  else 
    ret=-1
  endif
end function
function get_topography_file(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = topography_file

  ret=0
end function
function set_topography_file(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  topography_opt = 'file'
  topography_file = option
end function
function get_bottom_cell_file(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = bottom_cell_file

  ret=0
end function
function set_bottom_cell_file(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  partial_bottom_cells = .true.
  bottom_cell_file = option
end function


function get_ts_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = init_ts_option

  ret=0
end function
function set_ts_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'restart' .OR. option == 'internal') then
    init_ts_option = option
  else 
    ret=-1
  endif
end function
function get_ts_file(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = init_ts_file

  ret=0
end function
function set_ts_file(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  init_ts_option = 'restart'
  init_ts_file = option

  ! prevent that this setting can be overwritten by a pointer file
  luse_pointer_files = .false.
end function
function get_ts_file_format(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = init_ts_file_fmt

  ret=0
end function
function set_ts_file_format(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'bin' .OR. option == 'nc') then
    init_ts_file_fmt = option
  else 
    ret=-1
  endif
end function

function set_nprocs(nprocs) result (ret)
  integer :: ret
  integer, intent(in) :: nprocs
  ret=0

  nprocs_clinic = nprocs
  nprocs_tropic = nprocs
end function



function get_distribution(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = clinic_distribution_type

  ret=0
end function
function set_distribution(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'cartesian' .OR. option == 'predefined') then
    clinic_distribution_type = option
    tropic_distribution_type = option
  else
    ret=1
  endif
end function
function get_distribution_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  ret=0

  filename = distribution_file
end function
function set_distribution_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  ret=0

  distribution_file = filename
  clinic_distribution_type = 'predefined'
  tropic_distribution_type = 'predefined'
end function

function get_ew_boundary_type(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = ew_boundary_type

  ret=0
end function
function set_ew_boundary_type(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'closed' .OR. option == 'cyclic') then
    ew_boundary_type = option
  else
    ret=1
  endif
end function
function get_ns_boundary_type(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = ns_boundary_type

  ret=0
end function
function set_ns_boundary_type(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'closed' .OR. option == 'cyclic' .OR. option == 'tripole') then
    ns_boundary_type = option
  else
    ret=1
  endif
end function



function get_restart_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = restart_freq_opt

  ret=0
end function
function set_restart_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'never' .OR. option == 'nyear' .OR. option == 'nmonth' & 
      .OR. option == 'nday' .OR. option == 'nhour' .OR. option == 'nstep') then
    restart_freq_opt = option
  else
    ret=1
  endif
end function
function get_restart_freq_option(option) result (ret)
  integer :: ret
  integer, intent(out) :: option
  
  option = restart_freq

  ret=0
end function
function set_restart_freq_option(option) result (ret)
  integer :: ret
  integer, intent(in) :: option
  ret=0

  restart_freq = option
end function
function get_restart_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  
  filename = restart_outfile

  ret=0
end function
function set_restart_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  ret=0

  restart_outfile = filename
end function



function get_tavg_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = tavg_freq_opt

  ret=0
end function
function set_tavg_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'never' .OR. option == 'nyear' .OR. option == 'nmonth' & 
      .OR. option == 'nday' .OR. option == 'nhour' .OR. option == 'nstep') then
    tavg_freq_opt = option
  else
    ret=1
  endif
end function
function get_tavg_freq_option(option) result (ret)
  integer :: ret
  integer, intent(out) :: option
  
  option = tavg_freq

  ret=0
end function
function set_tavg_freq_option(option) result (ret)
  integer :: ret
  integer, intent(in) :: option
  ret=0

  tavg_freq = option
end function
function get_tavg_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  
  filename = tavg_outfile

  ret=0
end function
function set_tavg_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  ret=0

  tavg_outfile = filename
end function


function get_movie_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = movie_freq_opt

  ret=0
end function
function set_movie_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'never' .OR. option == 'nyear' .OR. option == 'nmonth' & 
      .OR. option == 'nday' .OR. option == 'nhour' .OR. option == 'nstep') then
    movie_freq_opt = option
  else
    ret=1
  endif
end function
function get_movie_freq_option(option) result (ret)
  integer :: ret
  integer, intent(out) :: option
  
  option = movie_freq

  ret=0
end function
function set_movie_freq_option(option) result (ret)
  integer :: ret
  integer, intent(in) :: option
  ret=0

  movie_freq = option
end function
function get_movie_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  
  filename = movie_outfile

  ret=0
end function
function set_movie_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  ret=0

  movie_outfile = filename
end function


function get_runid(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = runid

  ret=0
end function
function set_runid(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  runid = option
end function
function get_dt_option(option) result (ret)
  integer :: ret
  character (char_len), intent(out) :: option
  
  option = dt_option

  ret=0
end function
function set_dt_option(option) result (ret)
  integer :: ret
  character (char_len), intent(in) :: option
  ret=0

  if (option == 'auto_dt' .OR. option == 'steps_per_year' .OR. option == 'steps_per_day' & 
      .OR. option == 'seconds' .OR. option == 'hours') then
    dt_option = option
  else
    ret=1
  endif
end function
function get_dt_count(option) result (ret)
  integer :: ret
  integer, intent(out) :: option
  
  option = dt_count

  ret=0
end function
function set_dt_count(option) result (ret)
  integer :: ret
  integer, intent(in) :: option
  ret=0

  dt_count = option
end function

function change_directory(pathname) result (ret)
  integer :: ret
  character (char_len), intent(in) :: pathname

  call CHDIR(pathname)
  ret=0
end function


function get_shf_filename(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  filename = shf_filename
  ret=0
end function
function get_shf_data_type(type_) result (ret)
  integer :: ret
  character (char_len), intent(out) :: type_
  type_ = shf_data_type
  ret=0
end function
function set_shf_monthly_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  shf_data_type   = 'monthly'
  shf_interp_freq = 'every-timestep'
  shf_interp_type = 'linear'
  shf_filename    = filename
  ret=0
end function

function get_sfwf_filename(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  filename = sfwf_filename
  ret=0
end function
function get_sfwf_data_type(type_) result (ret)
  integer :: ret
  character (char_len), intent(out) :: type_
  type_ = sfwf_data_type
  ret=0
end function
function set_sfwf_monthly_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  sfwf_data_type   = 'monthly'
  sfwf_interp_freq = 'every-timestep'
  sfwf_interp_type = 'linear'
  sfwf_filename    = filename
  ret=0
end function

!imposed freshwater flux, only signifcant when ladjust_precip = .true.
!value is in Sverdrups
function get_fwf_imposed(fwf_) result (ret)
  integer :: ret
  real*8, intent(out) :: fwf_
  fwf_ = fwf_imposed
  ret=0
end function
function set_fwf_imposed(fwf_) result (ret)
  integer :: ret
  real*8, intent(in) :: fwf_
  fwf_imposed = fwf_
  ret=0
end function

function get_ws_filename(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  filename = ws_filename
  ret=0
end function
function get_ws_data_type(type_) result (ret)
  integer :: ret
  character (char_len), intent(out) :: type_
  type_ = ws_data_type
  ret=0
end function
function set_ws_monthly_file(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  ws_data_type    = 'monthly'
  ws_interp_freq  = 'every-timestep'
  ws_interp_type  = 'linear'
  ws_filename = filename
  ret=0
end function


function get_namelist_filename(filename) result (ret)
  integer :: ret
  character (char_len), intent(out) :: filename
  filename = nml_filename
  ret=0
end function
function set_namelist_filename(filename) result (ret)
  integer :: ret
  character (char_len), intent(in) :: filename
  nml_filename = filename
  ret=0
end function




! this routine was created because, for coupling purposes, it is important
! that the lon,lat positions of the complete grid can be requested, including
! the land-only blocks not present in the rest of the POP simulation
subroutine initialize_global_grid

    select case (horiz_grid_opt)
    case ('internal')
      call horiz_grid_internal(.true.)
    case ('file')
      call broadcast_scalar(horiz_grid_file, master_task)
      call read_horiz_grid(horiz_grid_file,.true.)
    case default
      call exit_POP(sigAbort,'ERROR: unknown horizontal grid option')
    end select

    !read_horiz_grid allocates ULAT_G and ULON_G on all nodes
    !deallocate on all but the master
    if (my_task /= master_task) then
        deallocate(ULAT_G, ULON_G)
    endif

    !compute the lat,lon of T points on master
    if (my_task == master_task) then
        allocate (TLAT_G(nx_global,ny_global), &
              TLON_G(nx_global,ny_global))

        call calc_tpoints_global
    endif
end subroutine initialize_global_grid

! copied and modified from grid.F90 to work for global grid
subroutine calc_tpoints_global
! !DESCRIPTION:
!  Calculates lat/lon coordinates of T points from U points
!  using a simple average of four neighbors in Cartesian 3d space.
!
   integer (POP_i4) :: i,j,n

   real (POP_r8) ::                   &
      xc,yc,zc,xs,ys,zs,xw,yw,zw, &! Cartesian coordinates for
      xsw,ysw,zsw,tx,ty,tz,da      !    nbr points

!-----------------------------------------------------------------------
!
!  TLAT_G, TLON_G are southwest 4-point averages of ULAT_G,ULON_G
!  for general grids, must drop into 3-d Cartesian space to prevent
!  problems near the pole
!
!-----------------------------------------------------------------------

    do j=2,ny_global
    do i=2,nx_global

         !***
         !*** convert neighbor U-cell coordinates to 3-d Cartesian coordinates 
         !*** to prevent problems with averaging near the pole
         !***

         zsw = cos(ULAT_G(i-1,j-1))
         xsw = cos(ULON_G(i-1,j-1))*zsw
         ysw = sin(ULON_G(i-1,j-1))*zsw
         zsw = sin(ULAT_G(i-1,j-1))

         zs  = cos(ULAT_G(i  ,j-1))
         xs  = cos(ULON_G(i  ,j-1))*zs
         ys  = sin(ULON_G(i  ,j-1))*zs
         zs  = sin(ULAT_G(i  ,j-1))

         zw  = cos(ULAT_G(i-1,j  ))
         xw  = cos(ULON_G(i-1,j  ))*zw
         yw  = sin(ULON_G(i-1,j  ))*zw
         zw  = sin(ULAT_G(i-1,j  ))

         zc  = cos(ULAT_G(i  ,j  ))
         xc  = cos(ULON_G(i  ,j  ))*zc
         yc  = sin(ULON_G(i  ,j  ))*zc
         zc  = sin(ULAT_G(i  ,j  ))

         !***
         !*** straight 4-point average to T-cell Cartesian coords
         !***

         tx = p25*(xc + xs + xw + xsw)
         ty = p25*(yc + ys + yw + ysw)
         tz = p25*(zc + zs + zw + zsw)

         !***
         !*** convert to lat/lon in radians
         !***

         da = sqrt(tx**2 + ty**2 + tz**2)

         TLAT_G(i,j) = asin(tz/da)

         if (tx /= c0 .or. ty /= c0) then
            TLON_G(i,j) = atan2(ty,tx)
         else
            TLON_G(i,j) = c0
         endif

    end do
    end do

    !***
    !*** for bottom row of domain where sw 4pt average is not valid,
    !*** extrapolate from interior
    !*** NOTE: THIS ASSUMES A CLOSED SOUTH BOUNDARY - WILL NOT
    !***       WORK CORRECTLY FOR CYCLIC OPTION
    !***

    do i=1,nx_global
       TLON_G(i,1) = TLON_G(i,2)
       TLAT_G(i,1) = c2*TLAT_G(i,2) - TLAT_G(i,3)
    end do

    where (TLON_G(:,:) > pi2) TLON_G(:,:) = TLON_G(:,:) - pi2
    where (TLON_G(:,:) < c0 ) TLON_G(:,:) = TLON_G(:,:) + pi2

    where (TLON_G(:,:) > pi2) TLON_G(:,:) = TLON_G(:,:) - pi2
    where (TLON_G(:,:) < c0 ) TLON_G(:,:) = TLON_G(:,:) + pi2

end subroutine calc_tpoints_global

end module pop_interface

