module MOM_dart_driver

use MOM_coms, only : set_PElist, broadcast, root_PE
use MOM_domains, only : global_field, pass_var
use MOM_ensemble_manager, only : get_ensemble_id, get_ensemble_size
use MOM_ensemble_manager, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_error_handler, only : is_root_pe, MOM_error, FATAL
use MOM_file_parser, only : get_param, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_time_manager, only : get_time, get_date, time_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

use mpp_mod, only: mpp_declare_pelist ! TODO: replace with wrap

use filter_mod, only: online_filter_init, online_filter_assim

implicit none; private
public :: dart_CS, dart_init, dart_assim_step, dart_end

character(len=40)  :: mdl = "MOM_dart_driver"
#include <MOM_memory.h>

type :: dart_CS
  integer :: nk = 1
  integer :: interval_sec
  integer(kind=8) :: nxg, nyg, nz
  integer(kind=8) :: model_size
  character(len=200) :: obs_dir
  real(kind=8), dimension(:), allocatable :: X
end type

contains

subroutine dart_init(G, GV, param_file, CS)
  type(ocean_grid_type), pointer :: G  !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Pointer to DA vrtical grid
  type(param_file_type), intent(in) :: param_file
  type(dart_CS), pointer :: CS

  integer :: ens_info(6)
  integer :: npes_pm
  integer :: ensemble_size !< Size of the ensemble
  integer :: this_id = 0 !< id of the current ensemble member
  integer :: dart_comm
  integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
  integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members

  allocate(CS)
  call get_param(param_file, mdl, "DART_ASSIM_INTERVAL", CS%interval_sec, &
    "Time between DART analysis steps.", units="s", default=86400)
  call get_param(param_file, mdl, "DART_OBS_DIR", CS%obs_dir, &
    "Directory with obs_seq files", default="obs")

  ! similar to MOM_oda_driver.F90
  ens_info = get_ensemble_size()
  ensemble_size = ens_info(1)
  npes_pm = ens_info(3)
  this_id = get_ensemble_id()
  allocate(ensemble_pelist(ensemble_size,npes_pm))
  allocate(filter_pelist(ensemble_size*npes_pm))
  call get_ensemble_pelist(ensemble_pelist, 'ocean')
  call get_ensemble_filter_pelist(filter_pelist, 'ocean')

  call set_PElist(filter_pelist)
  call mpp_declare_pelist(ensemble_pelist(:, 1), 'dart_roots', dart_comm)
  call set_PElist(ensemble_pelist(this_id, :))

  CS%nxg = G%ieg - G%isg + 1
  CS%nyg = G%jeg - G%jsg + 1
  CS%nz = GV%ke

  if (is_root_PE()) then
    call online_filter_init(dart_comm, CS%model_size)
    if (CS%model_size .ne. 2_8 * CS%nxg * CS%nyg * CS%nz) then
      ! Factor of 2 accounts for 2 variables (T and S). If size disagrees
      ! DART could be trying to assimilate more variables (change in input.nml).
      write(*, *) "Model size: ", CS%model_size
      write(*, *) "x, y, z", CS%nxg, CS%nyg, CS%nz
      call MOM_error(FATAL, "dart_init: Size of DART state is inconsistent with MOM grid")
    endif
    allocate(CS%X(CS%model_size))
  endif

end subroutine dart_init

subroutine dart_assim_step(Time, G, GV, US, tv, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ocean_grid_type), pointer :: G  !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Pointer to DA vertical grid
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(dart_CS), pointer :: CS

  integer :: nxy
  real, dimension(:, :), allocatable :: global_2d
  integer :: year, month, day, hr, min, sec
  character(len=300) :: obs_file

  call get_time(Time, sec, day)
  if (mod(sec, CS%interval_sec) .ne. 0) return
  call get_date(Time, year, month, day, hr, min, sec)
  write(obs_file, '(a,"/obs_seq.",i4.4,3i2.2,".out")') trim(CS%obs_dir), year, month, day, hr

  ! Total horizontal grid size
  nxy = CS%nxg * CS%nyg
  allocate(global_2d(CS%nxg, CS%nyg))

  call gather_t(G, tv%S, US%S_to_ppt, 0, CS, global_2d) ! offset 0 (salinity is first)
  call gather_t(G, tv%T, US%C_to_degC, nxy*CS%nz, CS, global_2d) ! temperature second

  ! TODO: assim step here
  if (is_root_PE()) then
    call online_filter_assim(CS%X, year, month, day, hr, min, sec, trim(obs_file))
  endif

  call scatter_t(G, tv%S, US%S_to_ppt, 0, CS, global_2d, floor_zero=.true.) ! offset 0 (salinity is first)
  call scatter_t(G, tv%T, US%C_to_degC, nxy*CS%nz, CS, global_2d) ! temperature second

  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  deallocate(global_2d)

end subroutine dart_assim_step

subroutine gather_t(G, a, scale, off, CS, global_2d)
  type(ocean_grid_type), pointer :: G  !< ocean grid structure
  type(dart_CS), pointer :: CS
  real, dimension(SZI_(G),SZJ_(G), CS%nk), intent(in) :: a
  real :: scale
  integer, intent(in) :: off
  real, intent(inout), dimension(:, :) :: global_2d
  integer :: k, nxy
  ! Total horizontal grid size
  nxy = CS%nxg * CS%nyg
  ! similar to MOM_oda_driver
  do k=1, CS%nz
    global_2d(:, :) = 0.0
    call global_field(G%Domain%mpp_domain,  a(:, :, k), global_2d)
    if (is_root_PE()) then
      CS%X(off+(k-1)*nxy+1:off+k*nxy) = scale * reshape(global_2d, (/nxy/))
    endif
  enddo
end subroutine gather_t

subroutine scatter_t(G, a, scale, off, CS, global_2d, floor_zero)
  type(ocean_grid_type), pointer :: G  !< ocean grid structure
  type(dart_CS), pointer :: CS
  real, dimension(SZI_(G),SZJ_(G), CS%nk), intent(out) :: a
  real :: scale
  integer, intent(in) :: off
  real, intent(inout), dimension(:, :) :: global_2d
  logical, intent(in), optional :: floor_zero
  integer :: i, j, k, nxy
  ! Total horizontal grid size
  nxy = CS%nxg * CS%nyg
  ! similar to MOM_oda_driver
  do k=1, CS%nz
    if (is_root_PE()) then
      global_2d = reshape(CS%X(off+(k-1)*nxy+1:off+k*nxy), (/CS%nxg, CS%nyg/))
    endif
    call broadcast(global_2d, nxy, root_PE())
    do j=G%jsc, G%jec; do i=G%isc, G%iec;
      a(i, j, k) = global_2d(i+G%idg_offset-G%isg+1, j+G%jdg_offset-G%jsg+1) / scale
      if (present(floor_zero)) then
        if (floor_zero) a(i, j, k) = max(a(i, j, k), 0.)
      endif
    enddo; enddo
  enddo
end subroutine scatter_t

subroutine dart_end(CS)
  type(dart_CS), pointer :: CS
  if (associated(CS)) then
    if (allocated(CS%X)) deallocate(CS%X)
    deallocate(CS)
  endif
  !write(*,*) "test"
end subroutine dart_end

end module MOM_dart_driver