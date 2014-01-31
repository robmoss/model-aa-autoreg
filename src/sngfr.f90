!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

!> Calculate the afferent arteriole response over a range of conditions.
!! The results can be used as lookup tables for the afferent arteriole in
!! larger models, allowing SNGFR to be determined with no computational cost.
module sngfr

  use model_base, only: dp
  use model_glom, only: kPa_of_mmHg

  implicit none

  private

  public sngfr_table

  real(kind=dp), parameter, private :: aa_interp = 0.2_dp

contains

  ! Write a lookup table of net afferent resistance for a range of myogenic
  ! response slopes, TGF activation levels, and arterial pressures.
  subroutine write_table(unit)

    use model_base

    integer, intent(in) :: unit

    integer, parameter :: nsegs = 100
    real(kind=dp), parameter :: slopes(2) = (/ 1.0_dp, 2.0_dp /)
    real(kind=dp), parameter :: alphas(2) = (/ 1.0_dp, 10.0_dp /)
    type(aa_params) :: p
    type(aa_state) :: st
    real(kind=dp) :: ptgf

    integer :: slope, alpha, tgf, pa

    open(unit, file="sngfr_lookup_base.ssv")
    write (unit, "(4A6, A12)") "slope", "alpha", "ptgf", "rpp", "res"
    do slope = 1, 2
       p = aa_params_init(shape=1, myo_slope=slopes(slope))
       do alpha = 2, 2
          p%alpha = alphas(alpha)
          do tgf = 2, 22
             ptgf = real(tgf, dp) / 200_dp
             do pa = 40, 180
                ! Convert the pressure from mmHg to kPa
                p%Pa = kPa_of_mmHg(real(pa, dp))
                st = aa_state_of(p, aa_interp, nsegs, st, ptgf=ptgf)
                write (unit, "(2F6.2, F6.3, F6.1, ES20.12)") &
                     slopes(slope), alphas(alpha), ptgf, real(pa, dp), &
                     sum(st%res)
                call flush(unit)
             end do
          end do
       end do
    end do
    close(unit)
  end subroutine write_table


  ! Write a lookup table of net afferent resistance for a range of myogenic
  ! response slopes, TGF activation levels, and arterial pressures.
  subroutine write_table_glom(unit)

    use model_glom

    integer, intent(in) :: unit
    integer, parameter :: nsegs = 100, shape = 1
    real(kind=dp), parameter :: myo_slope = 2.5_dp, alpha = 10_dp
    ! Kf chosen such that SNGFR ~= 30 and 36 nL/min, respectively.
    real(kind=dp), parameter :: kfs(2) = (/ 1.6_dp, 1.9_dp /)
    type(glom_params) :: p
    type(glom_state) :: st

    integer :: kf, tgf, pa
    real(kind=dp) :: ptgf

    p = glom_params_init(shape=shape, myo_slope=myo_slope)
    p%aa%k = p%aa%k * 5_dp
    p%aa%alpha = alpha

    open(unit, file="sngfr_lookup_glom.ssv")
    write (unit, "(3A6, 8A20)") "neph", "ptgf", "rpp", "sngfr", "aa_res", &
         "pg0", "pg1", "aff_bf", "aff_pf", "eff_bf", "eff_pf"

    do kf = 1, 2
       p%kf = kfs(kf)

       do tgf = 2, 22
          ptgf = real(tgf, dp) / 200_dp

          do pa = 40, 130
             ! Convert the pressure from mmHg to kPa
             p%aa%Pa = kPa_of_mmHg(real(pa, dp))

             st = glom_state_of(p, aa_interp, nsegs, ptgf, st)

             write (unit, "(I6, F6.3, I6, 8ES20.12)") &
                  kf, ptgf, pa, st%sngfr, st%aff_res, &
                  st%pg0, st%pg1, st%aff_bf, st%aff_pf, st%eff_bf, st%eff_pf
             call flush(unit)
          end do

          ! Reset the normalised radius when finished iterating over pa.
          st%aa%rad = 1.0_dp
       end do
    end do

    close(unit)
  end subroutine write_table_glom

  subroutine sngfr_table(opts)
    use cmdline, only: options
    type(options), intent(in) :: opts
    integer, parameter :: unit = 42

    if (opts%glomerulus) then
       call write_table_glom(unit)
    else
       call write_table(unit)
    end if
  end subroutine sngfr_table

end module sngfr
