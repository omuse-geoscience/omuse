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
  use constants, only: grav
  use domain, only: distrb_clinic
  use forcing_coupled, only: rotate_wind_stress
  use forcing_fields, only: SMF, SMFT, lsmft_avail
  use prognostic, only: TRACER, PSURF, UVEL, VVEL, RHO, curtime
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

  real (r8), dimension (nx_block, ny_block, max_blocks_clinic) :: &
     tau_x, tau_y    ! wind stress on t grid in N/m**2

  real (r8), dimension (:,:), allocatable :: & 
     WORK_G            ! temporary 2D array for gathered data



contains

!-----------------------------------------------------------------------
!
! initialize the model run
!
!-----------------------------------------------------------------------
function initialize_code() result(ret)
  integer :: ret

  !call POP_CommInitMessageEnvironment

  errorCode = POP_Success

  call POP_Initialize(errorCode)

  fstop_now = get_time_flag_id('stop_now')
  fupdate_coriolis = .false.
  fupdate_wind_stress = .false.

  ! allocate space for global gathers
  !if (my_task == master_task) then
    allocate(WORK_G(nx_global,ny_global))
  !endif

  !-----------------------------------------------------------------------
  !
  ! start up the main timer
  !
  !-----------------------------------------------------------------------

  call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)
  call timer_start(timer_total)

  call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)

  ret=0
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
function get_node_wind_stress(g_i, g_j, tau_x_, tau_y_) result(ret)
  integer, intent(in) :: g_i, g_j
  real*8, intent(out) :: tau_x_, tau_y_
  integer :: ret

  call get_gridded_variable(g_i, g_j, tau_x, tau_x_)
  call get_gridded_variable(g_i, g_j, tau_y, tau_y_)

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

  ! local variables
  integer :: iblock

  ! one problem may be that wind stress is zeroed above land points by
  ! rotate_wind_stress()

    if (lsmft_avail) then !use SMFT

      ! do nothing to prepare SMFT

    else  ! convert SMF into SMFT

      do iblock=1,nblocks_clinic
        call ugrid_to_tgrid(SMFT(:,:,1,iblock), SMF(:,:,1,iblock), iblock)
        call ugrid_to_tgrid(SMFT(:,:,2,iblock), SMF(:,:,2,iblock), iblock)
      enddo

      call POP_HaloUpdate(SMFT(:,:,1,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)

      call POP_HaloUpdate(SMFT(:,:,2,:),POP_haloClinic,  &
                      POP_gridHorzLocCenter,          &
                      POP_fieldKindVector, errorCode, &
                      fillValue = 0.0_r8)
    endif

    ! SMFT has been prepared now compute tau_x and tau_y from SMFT

    ! the following computation is simply the inverted computation
    ! of what happens within rotate_wind_stress():
    !    SMFT(:,:,1,:) = (WORK1(:,:,:)*cos(ANGLET(:,:,:)) +  &
    !                     WORK2(:,:,:)*sin(ANGLET(:,:,:)))*  &
    !                     RCALCT(:,:,:)*momentum_factor
    !    SMFT(:,:,2,:) = (WORK2(:,:,:)*cos(ANGLET(:,:,:)) -  &
    !                     WORK1(:,:,:)*sin(ANGLET(:,:,:)))*  &
    !                     RCALCT(:,:,:)*momentum_factor
    ! where RCALCT is a 0 or 1 mask for land or ocean points
    ! and WORK1 and WORK2 are tau_x and tau_y
    tau_x =  SMFT(:,:,1,:)/momentum_factor * cos(ANGLET) - SMFT(:,:,2,:)/momentum_factor * sin(ANGLET)
    tau_y = ( (SMFT(:,:,1,:) / momentum_factor) - tau_x * cos(ANGLET) ) / sin(ANGLET)

  if (errorCode == POP_Success) then
    ret=0
  else
    ret=-1
  endif
end function
function set_node_wind_stress (g_i, g_j, tau_x_, tau_y_) result(ret)
  real*8, intent(in) :: tau_x_, tau_y_
  integer, intent(in) :: g_i, g_j ! global i and j indexes for the gridpoint
  integer :: ret

  call set_gridded_variable(g_i, g_j, tau_x, tau_x_)
  call set_gridded_variable(g_i, g_j, tau_y, tau_y_)

  fupdate_wind_stress = .true.

  ret=0
end function

function set_wind_stress() result(ret)
  integer :: ret

  fupdate_wind_stress = .false.

  ! rotate_wind_stress applies tau_x and tau_y to the U-grid
  ! variable SMF (surface momentum fluxes (wind stress))
  ! wind stress in N/m^2 is converted to vel flux (cm^2/s^2)
  ! the method expects tau_x and tau_y to be supplied for the T grid

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


!-----------------------------------------------------------------------
!
! Generic getter and setters for gridded variables
!
!-----------------------------------------------------------------------
subroutine get_gridded_variable(g_i, g_j, grid, value)
  integer, intent(in) :: g_i, g_j
  real*8, intent(in), dimension(nx_block, ny_block, max_blocks_clinic) :: grid
  real*8, intent(out) :: value

  integer :: i,j,iblock,ierr
  type (block) :: this_block

  include 'mpif.h'

  value = 0.0 !important for MPI_SUM used later

  ! search rather than iterate over all points
  do iblock=1, nblocks_clinic
    this_block = get_block(blocks_clinic(iblock),iblock)
  
    !if (my_task == master_task) then  !debug info
      !write(*,*) 'Master is looking for index', g_i, ',', g_j
      !write(*,*) 'Start of block i index is:', this_block%i_glob(1)
      !write(*,*) 'Start of block j index is:', this_block%j_glob(1)
      !write(*,*) 'Start of domain in block i:', this_block%i_glob(1+nghost)
      !write(*,*) 'Start of domain in block j:', this_block%i_glob(1+nghost)
      !write(*,*) 'End of block i index is:', this_block%i_glob(nx_block)
      !write(*,*) 'End of block j index is:', this_block%j_glob(nx_block)
      !write(*,*) 'End of domain in block i:', this_block%i_glob(nx_block-nghost)
      !write(*,*) 'End of domain in block j:', this_block%i_glob(ny_block-nghost)
    !endif  

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
      if (this_block%i_glob(i) /= g_i .or. this_block%j_glob(j) /= g_j) then
        write (*,*) 'Error: g_i =', g_i, ', g_j=', g_j, ' and found i=', this_block%i_glob(i), ' j=', this_block%j_glob(j) 
      endif

      value = grid(i,j,iblock)

    endif

  enddo

  !only the output of the node with MPI rank 0 is significant in AMUSE, reduce to single value
  if(my_task == master_task) then
    call MPI_REDUCE(MPI_IN_PLACE, value, 1, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  else
    call MPI_REDUCE(value, 0, 1, MPI_DBL, MPI_SUM, master_task, MPI_COMM_OCN, ierr)
  endif

!  if (ierr /= MPI_SUCCESS) then
!   ret=-1
!  else
!   ret=0
!  endif
end subroutine get_gridded_variable

subroutine set_gridded_variable(g_i, g_j, grid, value)
  integer, intent(in) :: g_i, g_j
  real*8, intent(out), dimension(nx_block, ny_block, max_blocks_clinic) :: grid
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
      if (this_block%i_glob(i) /= g_i .or. this_block%j_glob(j) /= g_j) then
        write (*,*) 'Error: g_i =', g_i, ', g_j=', g_j, ' and found i=', this_block%i_glob(i), ' j=', this_block%j_glob(j) 
      endif

      grid(i,j,iblock) = value
      
    endif

  enddo

end subroutine set_gridded_variable


subroutine get_gridded_variable_3D(g_i, g_j, k, grid, value)
  integer, intent(in) :: g_i, g_j, k
  real*8, intent(in), dimension(nx_block, ny_block, km, max_blocks_clinic) :: grid
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
  real*8, intent(out), dimension(nx_block, ny_block, km, max_blocks_clinic) :: grid
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


!-----------------------------------------------------------------------
!
! Getter and setter for Coriolis Force
!
! these are currently implemented using the FCOR on the U grid
! set_coriolis_f() can be used to update the FCORT on the T grid based
! on the modified FCOR for the U grid
!
!-----------------------------------------------------------------------
function get_node_coriolis_f(g_i, g_j, corif_) result(ret)
  integer :: ret
  integer, intent(in) :: g_i, g_j
  real*8, intent(out) :: corif_

  call get_gridded_variable(g_i, g_j, FCOR, corif_)

  ret=0
end function
function set_node_coriolis_f(g_i, g_j, corif_) result(ret)
  integer :: ret
  integer, intent(in) :: g_i, g_j
  real*8, intent(out) :: corif_

  call set_gridded_variable(g_i, g_j, FCOR, corif_)

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
function get_node_surface_state(g_i, g_j, ssh_, uvel_, vvel_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: ssh_, uvel_, vvel_
  integer :: ii

  !local variable:
  real*8 :: psurf_

  if (n == 1) then
    call get_gridded_variable(g_i(1), g_j(1), PSURF(:,:,curtime,:), psurf_)
    ssh_(1) = psurf_ / grav

    call get_gridded_variable(g_i(1), g_j(1), UVEL(:,:,1,curtime,:), uvel_(1))
    call get_gridded_variable(g_i(1), g_j(1), VVEL(:,:,1,curtime,:), vvel_(1))
  else
    call gather_global(WORK_G, PSURF(:,:,curtime,:), master_task, distrb_clinic)
    if (my_task == master_task) then
      do ii=1,n
        ssh_(ii) = WORK_G(g_i(ii),g_j(ii)) / grav
      enddo
    endif
    call gather_global(WORK_G, UVEL(:,:,1,curtime,:), master_task, distrb_clinic)
    if (my_task == master_task) then
      do ii=1,n
        uvel_(ii) = WORK_G(g_i(ii),g_j(ii))
      enddo
    endif
    call gather_global(WORK_G, VVEL(:,:,1,curtime,:), master_task, distrb_clinic)
    if (my_task == master_task) then
      do ii=1,n
        vvel_(ii) = WORK_G(g_i(ii),g_j(ii))
      enddo
    endif
 
  endif

  ret=0
end function
function get_element_surface_state(g_i, g_j, temp_, salt_, n) result (ret)
  integer :: ret
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: g_i, g_j
  real*8, dimension(n), intent(out) :: temp_, salt_
  integer :: ii

  if (n == 1) then
    !first 1 is depth level, second is tracer index (temp or salt)
    call get_gridded_variable(g_i(1), g_j(1), TRACER(:,:,1,1,curtime,:), temp_(1))
    call get_gridded_variable(g_i(1), g_j(1), TRACER(:,:,1,2,curtime,:), salt_(1))
  else
    call gather_global(WORK_G, TRACER(:,:,1,1,curtime,:), master_task, distrb_clinic)
    if (my_task == master_task) then
      do ii=1,n
        temp_(ii) = WORK_G(g_i(ii),g_j(ii))
      enddo
    endif
    call gather_global(WORK_G, TRACER(:,:,1,2,curtime,:), master_task, distrb_clinic)
    if (my_task == master_task) then
      do ii=1,n
        salt_(ii) = WORK_G(g_i(ii),g_j(ii))
      enddo
    endif
  endif

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
function get_node_position(g_i, g_j, lat_, lon_) result(ret)
  integer :: ret
  integer, intent(in) :: g_i, g_j
  real*8, intent(out) :: lat_, lon_

  call get_gridded_variable(g_i, g_j, ULAT, lat_)
  call get_gridded_variable(g_i, g_j, ULON, lon_)

  ret=0
end function

function get_element_position(g_i, g_j, lat_, lon_) result(ret)
  integer :: ret
  integer, intent(in) :: g_i, g_j
  real*8, intent(out) :: lat_, lon_

  call get_gridded_variable(g_i, g_j, TLAT, lat_)
  call get_gridded_variable(g_i, g_j, TLON, lon_)

  ret=0
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







!-----------------------------------------------------------------------
!
! Getters and setters for 3D fields
!
!-----------------------------------------------------------------------
function get_temperature(i, j, k, temp_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(out) :: temp_

  call get_gridded_variable_3D(i, j, k, TRACER(:,:,:,1,curtime,:), temp_)

  ret=0
end function
function set_temperature(i, j, k, temp_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(in) :: temp_

  call set_gridded_variable_3D(i, j, k, TRACER(:,:,:,1,curtime,:), temp_)

  ret=0
end function
function get_salinity(i, j, k, salt_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(out) :: salt_

  call get_gridded_variable_3D(i, j, k, TRACER(:,:,:,2,curtime,:), salt_)

  ret=0
end function
function set_salinity(i, j, k, salt_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(in) :: salt_

  call set_gridded_variable_3D(i, j, k, TRACER(:,:,:,2,curtime,:), salt_)

  ret=0
end function

function get_xvel(i, j, k, uvel_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(out) :: uvel_

  call get_gridded_variable_3D(i, j, k, UVEL(:,:,:,curtime,:), uvel_)

  ret=0
end function
function set_xvel(i, j, k, uvel_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(in) :: uvel_

  call set_gridded_variable_3D(i, j, k, UVEL(:,:,:,curtime,:), uvel_)

  ret=0
end function
function get_yvel(i, j, k, vvel_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(out) :: vvel_

  call get_gridded_variable_3D(i, j, k, VVEL(:,:,:,curtime,:), vvel_)

  ret=0
end function
function set_yvel(i, j, k, vvel_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(in) :: vvel_

  call set_gridded_variable_3D(i, j, k, VVEL(:,:,:,curtime,:), vvel_)

  ret=0
end function

function get_density(i, j, k, rho_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(out) :: rho_

  call get_gridded_variable_3D(i, j, k, RHO(:,:,:,curtime,:), rho_)

  ret=0
end function
function set_density(i, j, k, rho_) result (ret)
  integer :: ret
  integer, intent(in) :: i, j, k
  real*8, intent(in) :: rho_

  call set_gridded_variable_3D(i, j, k, RHO(:,:,:,curtime,:), rho_)

  ret=0
end function







end module pop_interface
