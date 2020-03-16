module OzoneLunaMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates ozone-induced stress for use with LUNA model.
  !
  ! will be updated
  !
  ! Developed by Stefanie Falk.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use clm_varcon  , only : spval
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use OzoneBaseMod, only : ozone_base_type
  use abortutils  , only : endrun

  implicit none
  save
  private

  ! !PUBLIC TYPES:
  type, extends(ozone_base_type), public :: ozone_luna_type
     private
     ! Private data members
     ! May need some more later on
     real(r8), pointer :: o3uptakesha_patch(:) ! ozone dose, shaded leaves (mmol O3/m^2)
     real(r8), pointer :: o3uptakesun_patch(:) ! ozone dose, sunlit leaves (mmol O3/m^2)
     real(r8), pointer, public :: o3coefvcmaxsha_patch(:)  ! ozone coefficient for max. carboxylation rate, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefvcmaxsun_patch(:)  ! ozone coefficient for max. carboxylation rate, sunlit leaves (0 - 1)
     real(r8), pointer, public :: o3coefjmaxsha_patch(:)  ! ozone coefficient for max. electron transport rate, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefjmaxsun_patch(:)  ! ozone coefficient for max. electron transport rate, sunlit leaves (0 - 1)
  

   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: Restart
     procedure, public :: CalcOzoneStress
     
     ! Private routines
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
  end type ozone_luna_type

  interface ozone_luna_type
     module procedure constructor
  end interface ozone_luna_type

  ! !PRIVATE TYPES:
  
  ! TODO(wjs, 2014-09-29) This parameter will eventually become a spatially-varying
  ! value, obtained from ATM
  real(r8), parameter :: forc_ozone = 100._r8 * 1.e-9_r8  ! ozone partial pressure [mol/mol]

  ! TODO(wjs, 2014-09-29) The following parameters should eventually be moved to the
  ! params file. Parameters differentiated on veg type should be put on the params file
  ! with a pft dimension.

  ! o3:h2o resistance ratio defined by Sitch et al. 2007
  ! This should actually be defined directly via the diffusivities in air DH2O : DO3
  real(r8), parameter :: ko3 = 1.67_r8

  ! A maximum damage may have to be set for the damage factor.
  real(r8), parameter :: o3_damage_limit = 0.4_r8

  ! Data is only available for broadleaf species as of now 2020-03
  ! o3 intercepts and slopes for JmaxO3/Jmax0
  real(r8), parameter :: needleleafJmaxInt   = 1._r8  ! units = unitless 
  real(r8), parameter :: needleleafJmaxSlope = 0._r8      ! units = per mmol m^-2
  real(r8), parameter :: broadleafJmaxInt    = 1._r8  ! units = unitless  
  real(r8), parameter :: broadleafJmaxSlope  = -0.0037_r8      ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyJmaxInt     = 1._r8  ! units = unitless
  real(r8), parameter :: nonwoodyJmaxSlope   = 0._r8 ! units = per mmol m^-2
  
  ! o3 intercepts and slopes for VcmaxO3/Vcmax0
  real(r8), parameter :: needleleafVcmaxInt   = 1._r8  ! units = unitless 
  real(r8), parameter :: needleleafVcmaxSlope = 0._r8      ! units = per mmol m^-2
  real(r8), parameter :: broadleafVcmaxInt    = 1._r8  ! units = unitless  
  real(r8), parameter :: broadleafVcmaxSlope  = -0.0093_r8      ! units = per mmol m^-2
  real(r8), parameter :: nonwoodyVcmaxInt     = 1._r8  ! units = unitless
  real(r8), parameter :: nonwoodyVcmaxSlope   = 0._r8 ! units = per mmol m^-2

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains
  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================

  !-----------------------------------------------------------------------
  function constructor() result(ozone_luna)
    !
    ! !DESCRIPTION:
    ! Return an instance of ozone_luna_type
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ozone_luna_type) :: ozone_luna  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'constructor'
    !-----------------------------------------------------------------------

    ! DO NOTHING (simply return a variable of the appropriate type)

    ! Eventually this should call the Init routine (or replace the Init routine
    ! entirely). But I think it would be confusing to do that until we switch everything
    ! to use a constructor rather than the init routine.
    
  end function constructor

!-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize ozone data structure
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !-----------------------------------------------------------------------

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    
  end subroutine Init

!-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for ozone data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    call this%InitAllocateBase(bounds)
    
    allocate(this%o3uptakesha_patch(begp:endp)) ; this%o3uptakesha_patch(:) = nan
    allocate(this%o3uptakesun_patch(begp:endp)) ; this%o3uptakesun_patch(:) = nan
    allocate(this%o3coefvcmaxsha_patch(begp:endp))  ; this%o3coefvcmaxsha_patch(:) = nan
    allocate(this%o3coefvcmaxsun_patch(begp:endp))  ; this%o3coefvcmaxsun_patch(:) = nan
    allocate(this%o3coefjmaxsha_patch(begp:endp))  ; this%o3coefjmaxsha_patch(:) = nan
    allocate(this%o3coefjmaxsun_patch(begp:endp))  ; this%o3coefjmaxsun_patch(:) = nan
    !allocate(this%tlai_old_patch(begp:endp))    ; this%tlai_old_patch(:) = nan

  end subroutine InitAllocate

 !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize ozone history variables
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp
    endp = bounds%endp

    this%o3uptakesun_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3UPTAKESUN', units='mmol/m^2', &
         avgflag='A', long_name='total ozone flux into sunlit leaves', &
         ptr_patch=this%o3uptakesun_patch)

    this%o3uptakesha_patch(begp:endp) = spval
    call hist_addfld1d (fname='O3UPTAKESHA', units='mmol/m^2', &
         avgflag='A', long_name='total ozone flux into shaded leaves', &
         ptr_patch=this%o3uptakesha_patch)

  end subroutine InitHistory
 
!-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for ozone
    !
    ! !ARGUMENTS:
    class(ozone_luna_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    call this%InitColdBase(bounds)

    this%o3uptakesha_patch(begp:endp) = 0._r8
    this%o3uptakesun_patch(begp:endp) = 0._r8
    this%o3coefvcmaxsha_patch(begp:endp) = 0._r8
    this%o3coefvcmaxsun_patch(begp:endp) = 0._r8
    this%o3coefjmaxsha_patch(begp:endp) = 0._r8
    this%o3coefjmaxsun_patch(begp:endp) = 0._r8
    
    !this%tlai_old_patch(begp:endp) = 0._r8

  end subroutine InitCold

!-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart of ozone variables.
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_inqvdlen, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(ozone_luna_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------
    !call restartvar(ncid=ncid, flag=flag, varname='o3_tlaiold', xtype=ncd_double, &
    !     dim1name='pft', &
    !     long_name='one-sided leaf area index, from previous timestep, for ozone calculations', units='', &
    !     readvar=readvar, interpinic_flag='interp', data=this%tlai_old_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone uptake for shaded leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3uptakesun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone uptake for sunlit leaves', units='mmol m^-3', &
         readvar=readvar, interpinic_flag='interp', data=this%o3uptakesun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefvcmaxsun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. carboxylation rate for sunlit leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefvcmaxsun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefjmaxsun', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. electron transport rate for sunlit leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefjmaxsun_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefvcmaxsha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. carboxylation rate for shaded leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefvcmaxsha_patch)

    call restartvar(ncid=ncid, flag=flag, varname='o3coefjmaxsha', xtype=ncd_double, &
         dim1name='pft', &
         long_name='ozone coefficient for max. electron transport rate for shaded leaves', units='unitless', &
         readvar=readvar, interpinic_flag='interp', data=this%o3coefjmaxsha_patch)

  end subroutine Restart

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine CalcOzoneStress(this, bounds, num_exposedvegp, filter_exposedvegp, &
       forc_pbot, forc_th, rssun, rssha, rb, ram, tlai)

    class(ozone_luna_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    integer  , intent(in) :: num_exposedvegp           ! number of points in filter_exposedvegp
    integer  , intent(in) :: filter_exposedvegp(:)     ! patch filter for non-snow-covered veg
    real(r8) , intent(in) :: forc_pbot( bounds%begc: ) ! atmospheric pressure (Pa)
    real(r8) , intent(in) :: forc_th( bounds%begc: )   ! atmospheric potential temperature (K)
    real(r8) , intent(in) :: rssun( bounds%begp: )     ! leaf stomatal resistance, sunlit leaves (s/m)
    real(r8) , intent(in) :: rssha( bounds%begp: )     ! leaf stomatal resistance, shaded leaves (s/m)
    real(r8) , intent(in) :: rb( bounds%begp: )        ! boundary layer resistance (s/m)
    real(r8) , intent(in) :: ram( bounds%begp: )       ! aerodynamical resistance (s/m)
    real(r8) , intent(in) :: tlai( bounds%begp: )      ! one-sided leaf area index, no burying by snow

    ! Enforce expected array sizes (mainly so that a debug-mode threaded test with
    ! ozone-off can pick up problems with the call to this routine)
    SHR_ASSERT_ALL((ubound(forc_pbot) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(forc_th) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rssun) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rssha) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(rb) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ram) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(tlai) == (/bounds%endp/)), errMsg(sourcefile, __LINE__))

    ! Explicitly set outputs to 1. This isn't really needed, because they should still be
    ! at 1 from cold-start initialization, but do this for clarity here.

    this%o3coefvsha_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefvsun_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefgsha_patch(bounds%begp:bounds%endp) = 1._r8
    this%o3coefgsun_patch(bounds%begp:bounds%endp) = 1._r8

  end subroutine CalcOzoneStress
  

end module OzoneLunaMod
