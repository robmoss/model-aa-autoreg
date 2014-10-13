!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

!> An extension of the afferent arteriolar model of renal autoregulation that
!! accounts for filtration along the glomerular capillary bed.
module model_glom

  use model_base, only: aa_params, aa_params_init, aa_state, aa_state_init

  implicit none

  !> The kind parameter for double-precision floating point numbers.
  integer, parameter :: dp = selected_real_kind(15, 307)

  !> Scaling factor to convert from kPa to mmHg.
  real(kind=dp), parameter :: kPa_to_mmHg = 7.50061683d0

  !> Scaling factor to convert from kPa . s / um^3 to mmHg . min / nLkPa.
  real(kind=dp), parameter :: res_kPa_to_mmHg = 1.2501028d5

  !> Model parameters specific to glomerular filtration.
  type glom_params
     !> Afferent arteriole parameters.
     type(aa_params) :: aa
     !> Glomerular properties: resistance of the glomerular capillary bed and
     !! efferent arteriole, glomerular filtration coefficient, and the
     !! hydrostatic pressure in Bowman's capsule.
     real(kind=dp) :: res_gc, res_ea, kf, p_bc
     !> Plasma properties: hematocrit and protein concentration (g/dl).
     real(kind=dp) :: hct, c_p
  end type glom_params

  !> State of the glomerular portion of the model, used by the iterative
  !! solver to determine the model SNGFR.
  type glom_sngfr
     !> Single-nephron glomerular filtration rate.
     real(kind=dp) :: sngfr
     !> Hydrostatic pressure at the beginning and end of the glomerulus.
     real(kind=dp) :: p_gc0, p_gc1
     !> Afferent blood and plasma flow.
     real(kind=dp) :: bf_aff, pf_aff
     !> Efferent blood and plasma flow.
     real(kind=dp) :: bf_eff, pf_eff
  end type glom_sngfr

  !> Model state.
  type glom_state
     !> The number of afferent arteriole segments.
     integer :: n
     !> The state of the afferent arteriole model.
     type(aa_state) :: aa
     !> Afferent blood and plasma flow.
     real(kind=dp) :: aff_bf, aff_pf
     !> Efferent blood and plasma flow.
     real(kind=dp) :: eff_bf, eff_pf
     !> Hydrostatic pressure at the beginning and end of the glomerulus.
     real(kind=dp) :: pg0, pg1
     !> Single-nephron glomerular filtration rate.
     real(kind=dp) :: sngfr
     !> Net resistance of the afferent arteriole.
     real(kind=dp) :: aff_res
  end type glom_state

contains

  !> Initialise the mode state.
  elemental function glom_state_init(n, realloc) result(s)
    integer, intent(in) :: n
    logical, intent(in), optional :: realloc
    type(glom_state) :: s

    if (present(realloc)) then
       s%aa = aa_state_init(n, realloc)
    else
       s%aa = aa_state_init(n)
    end if

    s%n = n
  end function glom_state_init


  !> Return the model parameters for "shape 1" (default) or "shape 2".
  elemental function glom_params_init(shape, myo_slope) result(p)
    integer, intent(in), optional :: shape
    real(kind=dp), intent(in), optional :: myo_slope
    type(glom_params) :: p
    integer :: sh
    real(kind=dp) :: slope

    real(kind=dp), parameter :: myo_op = 8.6_dp, myo_range = 7.4_dp

    if (present(shape)) then
       sh = shape
    else
       sh = 1
    end if

    if (present(myo_slope)) then
       slope = myo_slope
    else
       slope = 1_dp
    end if

    p%aa = aa_params_init(shape=shape, myo_slope=slope)

    ! Make the TGF response decay beyond the terminal 10% of the arteriole.
    p%aa%alpha = 10.0_dp

    p%res_gc = 0_dp
    p%res_ea = 0.4_dp
    ! Set kf = 1.6d0 --> SNGFR ~ 30 nL/min
    !     kf = 1.9d0 --> SNGFR ~ 36 nL/min
    p%kf = 1.6d0
    p%p_bc = 15.0d0

    p%hct = 0.45d0 ! Blood hematocrit.
    p%c_p = 5.5d0 ! Plasma protein concentration (g/dL).

  end function glom_params_init

  !> Convert from kPa to mmHg.
  elemental function mmHg_of_kPa(p_kPa) result(p_mmHg)
    real(kind=dp), intent(in) :: p_kPa
    real(kind=dp) :: p_mmHg

    p_mmHg = p_kPa * kPa_to_mmHg
  end function mmHg_of_kPa

  !> Convert from mmHg to kPa.
  elemental function kPa_of_mmHg(p_mmHg) result(p_kPa)
    real(kind=dp), intent(in) :: p_mmHg
    real(kind=dp) :: p_kPa

    p_kPa = p_mmHg / kPa_to_mmHg
  end function kPa_of_mmHg

  !> Convert from kPa . s / um^3 to mmHg . min / nLkPa,
  elemental function res_mmHg_of_kPa(res_kPa) result(res_mmHg)
    real(kind=dp), intent(in) :: res_kPa
    real(kind=dp) :: res_mmHg

    res_mmHg = res_kPa * res_kPa_to_mmHg
  end function res_mmHg_of_kPa

  !> Calculate the colloid oncotic pressure exerted by the plasma protein.
  !! @param[in] Cp The plasma protein concentration (g/dL).
  elemental function p_oncotic(Cp)
    real(kind=dp), intent(in) :: Cp
    real(kind=dp), parameter :: a = 2.1d0, b = 0.16d0, c = 0.009d0
    real(kind=dp) p_oncotic

    p_oncotic = a * Cp + b * Cp**2.0d0 + c * Cp**3.0d0
  end function p_oncotic

  !> Calculate the SNGFR (nL/min), efferent plasma flow (nL/min) and efferent
  !! hydrostatic pressure (mmHg).
  !! @param[in] p     The model parameters.
  !! @param[in] pg0   The initial hydrostatic pressure (mmHg).
  !! @param[in] c_p   The afferent protein concentration (g/dL).
  !! @param[in] q_aff The afferent plasma flow (nL/min).
  !! @param[in] b_aff The afferent blood flow (nL/min).
  elemental function sngfr_of(p, pg0, c_p, q_aff, b_aff)
    type(glom_params), intent(in) :: p
    real(kind=dp), intent(in) :: pg0, c_p, q_aff, b_aff
    type(glom_sngfr) :: sngfr_of

    ! The number of integration steps along the glomerular capillary bed.
    integer, parameter :: NS = 50
    ! The size of each integration step, assuming that the glomerular
    ! capillary bed has unit length.
    real(kind=dp), parameter :: dx = 1.0d0 / NS

    ! The hydrostatic pressure, plasma flow and blood flow along the
    ! glomerular capillary bed, and integration-step variables.
    real(kind=dp) :: PG_x(0:NS), QG_x(0:NS), BG_x(0:NS), p_onc, dQ, dP

    integer :: i

    ! The initial hydrostatic pressure and plasma flow are given.
    PG_x(0) = pg0
    QG_x(0) = q_aff
    BG_x(0) = b_aff

    do i = 1, NS
       p_onc = p_oncotic(c_p * q_aff / QG_x(i - 1))
       ! The decrease in plasma flow due to some fraction being filtered.
       dQ = - p%kf * (PG_x(i - 1) - p_onc - p%p_bc)
       ! The decrease in hydrostatic pressure due to blood flow.
       dP = - p%res_gc * BG_x(i - 1)
       QG_x(i) = QG_x(i - 1) + dx * dQ
       BG_x(i) = BG_x(i - 1) + dx * dQ
       PG_x(i) = PG_x(i - 1) + dx * dP
    end do

    sngfr_of%sngfr = q_aff - QG_x(NS)
    sngfr_of%p_gc0 = pg0
    sngfr_of%p_gc1 = PG_x(NS)
    sngfr_of%bf_aff = b_aff
    sngfr_of%pf_aff = q_aff
    sngfr_of%bf_eff = BG_x(NS)
    sngfr_of%pf_eff = QG_x(NS)

  end function sngfr_of


  !> Calculate the SNGFR (nL/min) for a nephron by solving for the hydrostatic
  !! pressure at the terminal end of the afferent arteriole.
  !! @param[in] p      The model parameters.
  !! @param[in] res_aa The net afferent arteriole resistance (mmHg/(nl/min)).
  function solve_sngfr(p, res_aa) result(soln)
    type(glom_params), intent(in) :: p
    real(kind=dp), intent(in) :: res_aa
    type(glom_sngfr) :: soln

    real(kind=dp), parameter :: toln = 1.0d-3
    real(kind=dp), parameter :: err_interp = 0.01d0
    real(kind=dp) :: l_pg, u_pg, l_err, u_err, m_pg, m_err, pa_mmHg, pv_mmHg
    integer :: i

    ! Convert from kPa to mmHg.
    pa_mmHg = mmHg_of_kPa(p%aa%Pa)
    pv_mmHg = mmHg_of_kPa(p%aa%Pv)

    ! For the initial bounds on PG(0), inspect intervals with widths of
    ! RPP * 5%, from 4% of RPP to 99% of RPP.
    do i = 1, 19
       l_pg = (real(i, dp) / 20_dp - 0.01_dp) * pa_mmHg
       u_pg = (real(i + 1, dp) / 20_dp - 0.01_dp) * pa_mmHg
       l_err = error_of_sngfr(sngfr_of_guess(l_pg))
       u_err = error_of_sngfr(sngfr_of_guess(u_pg))
       if (l_err * u_err < 0) then
          exit
       end if
    end do

    if (l_err * u_err > 0) then
       write (*,*) "    BISECTION ERROR", l_err, u_err
       call exit(1)
    else
       do while (abs(u_pg - l_pg) > 0.00001_dp)
          m_pg = 0.5_dp * (u_pg + l_pg)
          soln = sngfr_of_guess(m_pg)
          m_err = error_of_sngfr(soln)
          if (m_err * u_err > 0) then
             u_pg = m_pg
          else
             l_pg = m_pg
          end if
       end do

       return
    end if

  contains

      !> Given a guess for PG(0), calculate the SNGFR.
      elemental function sngfr_of_guess(pg0) result(soln)
        real(kind=dp), intent(in) :: pg0
        type(glom_sngfr) :: soln

        real(kind=dp) :: dPaff, Baff, Qaff

        dPaff = pa_mmHg - pg0
        Baff = dPaff / res_aa
        Qaff = Baff * (1.0d0 - p%hct)

        soln = sngfr_of(p, pg0, p%c_p, Qaff, Baff)
      end function sngfr_of_guess

      !> Given a calculated SNGFR, determine the conservation error.
      elemental function error_of_sngfr(soln) result(err)
        type(glom_sngfr), intent(in) :: soln
        real(kind=dp) :: err

        real(kind=dp) :: Baff, Beff, dPaff

        dPaff = pa_mmHg - soln%p_gc0
        Baff = dPaff / res_aa
        Beff = Baff - soln%sngfr

        err = soln%p_gc1 - Beff * p%res_ea - pv_mmHg

      end function error_of_sngfr

  end function solve_sngfr

  !> Calculate the steady-state profiles of radius, resistance, pressure and
  !! combined TGF and myogenic activation, along the afferent arteriole.
  function glom_state_of(p, interp_frac, segments, tgf_act, s_prev) result(s)
    type(glom_params), intent(in) :: p
    real(kind=dp), intent(in) :: interp_frac, tgf_act
    integer, intent(in) :: segments
    type(glom_state), intent(in), optional :: s_prev
    type(glom_state) :: s
    real(kind=dp) :: x, delta
    real(kind=dp), parameter :: toln = 9d-7 ! Tolerance for final solution.
    integer :: step, i

    s = glom_state_init(segments)

    ! Initialise the (normalised) inner radius to the reference value (1.0).
    if (present(s_prev) .and. s_prev%n == segments) then
       s%aa%rad = s_prev%aa%rad
    else
       s%aa%rad = 1.0_dp
    end if

    delta = update_radius(s%aa%rad)
    step = 1

    do while (delta > toln)
       delta = update_radius(s%aa%rad)
       step = step + 1
    end do

  contains

    !> Given the current (normalised) radius of each segment, update the
    !! net resistance, hydrostatic pressure, and the combined activation of
    !! the myogenic and TGF responses.
    !! Adjust the normalised radius of each segment to account for the change
    !! in activation, and return the maximal absolute change in (normalised)
    !! segment radius.
    function update_radius(radius) result(delta)

      use model_base, only: aa_seg_inner_rho, aa_myogenic_act, &
           aa_seg_norm_radius

      real(kind=dp), intent(inout) :: radius(segments)
      real(kind=dp) :: new_radius(segments), delta, x, res_aff, pa_glom
      integer :: i
      logical :: sum_mask(segments)
      type(glom_sngfr) :: soln

      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         ! The absolute radius of the segment.
         s%aa%abs_rad(i) = radius(i) * aa_seg_inner_rho(p%aa, x)
         ! The hemodynamic resistance of the segment.
         s%aa%res(i) = p%aa%k / (s%aa%abs_rad(i) ** 4_dp * segments)
      end do

      ! Calculate the net resistance of the afferent arteriole.
      res_aff = sum(s%aa%res)
      s%aff_res = res_mmHg_of_kPa(res_aff)

      ! Given this net resistance, determine what the hydrostatic pressure
      ! must be at the end of the afferent arteriole in order to satisfy
      ! conservation of energy.
      soln = solve_sngfr(p, s%aff_res)

      ! Update the model state.
      pa_glom = kPa_of_mmHg(soln%p_gc0)
      s%sngfr = soln%sngfr
      s%pg0 = soln%p_gc0
      s%pg1 = soln%p_gc1
      s%aff_bf = soln%bf_aff
      s%aff_pf = soln%pf_aff
      s%eff_bf = soln%bf_eff
      s%eff_pf = soln%pf_eff

      ! Calculate the hydrostatic pressure at the mid-point of each segment.
      sum_mask(:) = .false.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         if (i > 1) then
            sum_mask(i - 1) = .true.
         end if
         ! The hydrostatic pressure at the mid-point of the segment.
         s%aa%pa(i) = p%aa%Pa - (p%aa%Pa - pa_glom) &
              * (sum(s%aa%res, sum_mask) + s%aa%res(i) * 0.5_dp) / res_aff
      end do

      ! Calculate the combined activation of the myogenic and TGF responses.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         s%aa%myo_act(i) = aa_myogenic_act(p%aa, s%aa%pa(i))
         if (x >= (1_dp - 1_dp / p%aa%alpha)) then
            s%aa%tgf_act(i) = tgf_act * (1_dp - p%aa%alpha * (1_dp - x))
         else
            s%aa%tgf_act(i) = 0_dp
         end if
         s%aa%net_act(i) = s%aa%myo_act(i) + s%aa%tgf_act(i)
      end do

      ! Finally, update the (normalised) radius of each segment.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         new_radius(i) = aa_seg_norm_radius(p%aa, s%aa%net_act(i), &
              s%aa%pa(i), x)
      end do

      delta = maxval(abs(radius - new_radius))
      radius = (1_dp - interp_frac) * abs(radius) &
           + interp_frac * abs(new_radius)

    end function update_radius

  end function glom_state_of

end module model_glom
