!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

!> An implementation of the afferent arteriolar model of renal autoregulation
!! published by Feldberg et al., AJP Renal 269(4 Pt 2): F581-93, Oct 1995.
!! http://www.ncbi.nlm.nih.gov/pubmed/7485545
module model_base

  implicit none

  !> The kind parameter for double-precision floating point numbers.
  integer, parameter :: dp = selected_real_kind(15, 307)

  !> Model parameters, named as per Table 1 of Feldberg et al.
  type aa_params
     !> Constants that govern the passive stress due to elasticity.
     real(kind=dp) :: C1, C2, alpha1, alpha2
     !> Constants the govern the active stress due to contraction.
     real(kind=dp) :: sigma_max, epsilon_max, epsilon_w
     !> Relationship between transmural pressure and psi.
     real(kind=dp) :: Pmax, Pmin
     !> Afferent and venous hydrostatic pressures.
     real(kind=dp) :: Pa, Pv
     !> Efferent arteriolar resistance.
     real(kind=dp) :: Re
     !> Afferent arteriole resistance scaling factor.
     real(kind=dp) :: k
     !> Afferent arteriole inner radius when there is no transmural pressure.
     real(kind=dp) :: radius_i_00, radius_i_09, radius_i_10
     !> Afferent arteriole outer radius when there is no transmural pressure.
     real(kind=dp) :: radius_o_00, radius_o_09, radius_o_10
     !> Whether to use a linear radius profile (Eq 15) or a non-linear radius
     !! profile that narrows distinctly at the glomerular end (Eq 16).
     logical :: radius_nl
     !> Parameters for the linear TGF response at x = 1.
     real(kind=dp) :: a, b
     !> Lower and upper bounds for the TGF activation.
     real(kind=dp) :: tgf_min, tgf_max
     !> Determines how rapidly the TGF activation decreases with distance from
     !! the glomerulus.
     real(kind=dp) :: alpha
     integer :: shape
  end type aa_params

  !> State of the spatially-distributed model.
  type aa_state
     !> The number of afferent arteriole segments.
     integer :: n
     !> The normalised and absolute radii, resistance, transmural pressure,
     !! net activation, myogenic activation and TGF activation for every
     !! segment of the afferent arteriole.
     real(kind=dp), dimension(:), allocatable :: rad, abs_rad, res, pa, &
          net_act, myo_act, tgf_act
     !> The net afferent blood flow.
     real(kind=dp) :: aff_bf
  end type aa_state

contains

  !> Initialise the model state.
  pure elemental function aa_state_init(n, realloc) result(s)
    integer, intent(in) :: n
    logical, intent(in), optional :: realloc
    type(aa_state) :: s

    if (allocated(s%rad) .and. present(realloc) .and. realloc) then
       deallocate(s%rad)
       deallocate(s%abs_rad)
       deallocate(s%res)
       deallocate(s%pa)
       deallocate(s%net_act)
       deallocate(s%myo_act)
       deallocate(s%tgf_act)
    end if

    if (.not. allocated(s%rad)) then
       s%n = n
       allocate(s%rad(n))
       allocate(s%abs_rad(n))
       allocate(s%res(n))
       allocate(s%pa(n))
       allocate(s%net_act(n))
       allocate(s%myo_act(n))
       allocate(s%tgf_act(n))
       ! Set the initial (normalised) radius of each segment to 1.
       s%rad = 1.0_dp
    end if
  end function aa_state_init


  !> Return the model parameters for "shape 1" (default) or "shape 2".
  pure function aa_params_init(shape, myo_slope) result(p)
    integer, intent(in), optional :: shape
    real(kind=dp), intent(in), optional :: myo_slope
    type(aa_params) :: p
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

    p%shape = sh

    p%C1 = 2.7_dp
    p%alpha1 = 2.0_dp
    p%C2 = 0.25_dp
    p%alpha2 = 12.0_dp

    p%sigma_max = 100.0_dp
    p%epsilon_max = 0.25_dp
    p%epsilon_w = 0.9_dp

    p%Pmin = myo_op - myo_range / slope
    p%Pmax = myo_op + myo_range / slope

    p%Pa = 12.5_dp
    p%Pv = 0.5_dp

    p%Re = 1.63d-6

    p%a = 1.45d-4
    p%b = - 0.884_dp
    p%tgf_min = 0.01_dp
    p%tgf_max = 0.11_dp
    p%alpha = 1.0_dp

    if (sh == 1) then
       p%k = 4.18d-4
       p%radius_nl = .false.
       p%radius_i_00 = 10.0_dp
       p%radius_i_09 = 9.1_dp
       p%radius_i_10 = 9.0_dp
       p%radius_o_00 = 13.0_dp
       p%radius_o_09 = 11.8_dp
       p%radius_o_10 = 11.7_dp
    else if (sh == 2) then
       p%k = 2.49d-4
       p%radius_nl = .true.
       p%radius_i_00 = 10.0_dp
       p%radius_i_09 = 9.0_dp
       p%radius_i_10 = 8.0_dp
       p%radius_o_00 = 13.0_dp
       p%radius_o_09 = 11.7_dp
       p%radius_o_10 = 11.5_dp
    else
       ! Indirectly signal an error by setting the vessel radii to zero.
       p%k = 0_dp
       p%radius_nl = .false.
       p%radius_i_00 = 0_dp
       p%radius_i_09 = 0_dp
       p%radius_i_10 = 0_dp
       p%radius_o_00 = 0_dp
       p%radius_o_09 = 0_dp
       p%radius_o_10 = 0_dp
    end if

    ! Adjust the parameters (as published) so that the model results match the
    ! figures from the paper.
    p%a = p%a * 1d3
    p%k = p%k * 10_dp
  end function aa_params_init

  !> Calculate the passive (elastic) strain (Eq 4).
  pure function aa_seg_passive_strain(p, epsilon) result(strain)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: epsilon
    real(kind=dp) :: strain

    strain = p%C1 * (exp(p%alpha1 * epsilon) - 1.0d0) &
         + p%C2 * (exp(p%alpha2 * epsilon) - 1.0d0)
  end function aa_seg_passive_strain

  !> Calculate the active strain (Eq 5).
  pure function aa_seg_active_strain(p, epsilon, psi) result(strain)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: epsilon, psi
    real(kind=dp) :: strain

    if (abs(epsilon - p%epsilon_max) < p%epsilon_w) then
       strain = psi * p%sigma_max * &
            (1.0_dp - abs((epsilon - p%epsilon_max) / p%epsilon_w))
    else
       strain = 0.0_dp
    end if
  end function aa_seg_active_strain

  !> Calculate the (baseline) inner radius of the arteriole segment.
  pure function aa_seg_inner_rho(p, x) result(radius)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: radius

    if (p%radius_nl) then
       if (x < 0.9) then
          radius = p%radius_i_00 &
               + (x / 0.9_dp) * (p%radius_i_09 - p%radius_i_00)
       else
          radius = p%radius_i_09 &
               + ((x - 0.9_dp) / 0.1_dp) * (p%radius_i_10 - p%radius_i_09)
       end if
    else
       radius = p%radius_i_00 + x * (p%radius_i_10 - p%radius_i_00)
    end if
  end function aa_seg_inner_rho

  !> Calculate the (baseline) outer radius of the arteriole segment.
  pure function aa_seg_outer_rho(p, x) result(radius)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: radius

    if (p%radius_nl) then
       if (x < 0.9) then
          radius = p%radius_o_00 &
               + (x / 0.9_dp) * (p%radius_o_09 - p%radius_o_00)
       else
          radius = p%radius_o_09 &
               + ((x - 0.9_dp) / 0.1_dp) * (p%radius_o_10 - p%radius_o_09)
       end if
    else
       radius = p%radius_o_00 + x * (p%radius_o_10 - p%radius_o_00)
    end if
  end function aa_seg_outer_rho

  !> Calculate the net transmural pressure across the arterial wall (Eq 13).
  pure function aa_seg_transmural(p, psi, radius, aa_x) result(pressure)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: psi, radius
    real(kind=dp), intent(in), optional :: aa_x
    real(kind=dp) :: pressure

    integer, parameter :: steps = 100
    integer :: step
    real(kind=dp) :: x, z, z_min, z_max, dz
    real(kind=dp) :: epsilon, denom
    real(kind=dp) :: p_el(0:steps), p_act(0:steps)

    if (present(aa_x)) then
       x = aa_x
    else
       x = 0.0_dp
    end if

    z_min = 1.0_dp
    z_max = aa_seg_outer_rho(p, x) / aa_seg_inner_rho(p, x)
    dz = (z_max - z_min) / real(steps, dp)

    do step = 0, steps
       z = z_min + (z_max - z_min) * real(step, dp) / real(steps, dp)
       epsilon = (1_dp / z) * sqrt(radius**2_dp - 1_dp + z**2_dp) - 1_dp
       denom = z * (1_dp + epsilon)

       p_el(step) = dz * aa_seg_passive_strain(p, epsilon) / denom
       p_act(step) = dz * aa_seg_active_strain(p, epsilon, 1.0_dp) / denom
    end do

    ! Use the trapezoidal approximation of the integral.
    p_el(0) = 0.5_dp * p_el(0)
    p_el(steps) = 0.5_dp * p_el(steps)
    p_act(0) = 0.5_dp * p_act(0)
    p_act(steps) = 0.5_dp * p_act(steps)

    pressure = sum(p_el) + psi * sum(p_act)

  end function aa_seg_transmural

  !> Calculate the activation of the myogenic response (Eq 14).
  pure function aa_myogenic_act(p, pressure) result(psi)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: pressure
    real(kind=dp) :: psi

    if (pressure < p%Pmin) then
       psi = 0_dp
    else if (pressure > p%Pmax) then
       psi = 1_dp
    else
       psi = (pressure - p%Pmin) / (p%Pmax - p%Pmin)
    end if
  end function aa_myogenic_act

  !> Calculate the activation of the TGF response (Eqs 18, 19).
  pure function aa_tgf_act(p, pressure_glom, aa_x) result(psi)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: pressure_glom
    real(kind=dp), intent(in), optional :: aa_x
    real(kind=dp) :: psi, x

    if (present(aa_x)) then
       x = aa_x
    else
       x = 0.0_dp
    end if

    ! Calculate the TGF activation at the distal end of the arteriole (x=1).
    psi = p%a * pressure_glom + p%b
    ! Restrict the TGF activation at x=1 to the defined bounds.
    psi = max(p%tgf_min, min(p%tgf_max, psi))

    ! Account for the distance from the glomerulus.
    if (x >= (1_dp - 1_dp / p%alpha)) then
       psi = psi * (1_dp - p%alpha * (1_dp - x))
    else
       psi = 0_dp
    end if
  end function aa_tgf_act

  !> Given the current activation (psi) and transmural pressure (pressure),
  !! determine the normalised radius of the AA segment (satisfying Eq 13).
  function aa_seg_norm_radius(p, psi, pressure, x) result(soln)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: psi, pressure
    real(kind=dp), intent(in), optional :: x

    real(kind=dp), parameter :: bisect_toln = 1d-10
    real(kind=dp) :: lbound, ubound, midp, error, soln
    real(kind=dp) :: lerr, uerr
    integer :: i

    lbound = 0.001_dp
    ubound = 1.8_dp
    if (present(x)) then
       lerr = aa_seg_transmural(p, psi, lbound, x) - pressure
       uerr = aa_seg_transmural(p, psi, ubound, x) - pressure
    else
       lerr = aa_seg_transmural(p, psi, lbound) - pressure
       uerr = aa_seg_transmural(p, psi, ubound) - pressure
    end if

    if (lerr * uerr > 0_dp) then
       ! They have identical signs.
       write (*, *) "Bisection error: both bounds have the same error sign."
       if (present(x)) then
          write (*,*) psi, pressure, x
       else
          write (*,*) psi, pressure
       end if

       call exit(1)
    end if

    do while (abs(ubound - lbound) > bisect_toln)
       midp = 0.5_dp * (ubound + lbound)
       if (present(x)) then
          error = aa_seg_transmural(p, psi, midp, x) - pressure
       else
          error = aa_seg_transmural(p, psi, midp) - pressure
       end if

       if (error < 0) then
          lbound = midp
       else
          ubound = midp
       end if

    end do

    soln = 0.5_dp * (ubound + lbound)

  end function aa_seg_norm_radius

  !> Calculate the steady-state profiles of radius, resistance, pressure and
  !! combined TGF and myogenic activation, along the afferent arteriole.
  function aa_state_of(p, interp_frac, segments, prev_state, ptgf) result(s)
    type(aa_params), intent(in) :: p
    real(kind=dp), intent(in) :: interp_frac
    integer, intent(in) :: segments
    type(aa_state), intent(in), optional :: prev_state
    real(kind=dp), intent(in), optional :: ptgf
    type(aa_state) :: s
    real(kind=dp) :: x, delta, old_r0
    real(kind=dp), parameter :: toln = 1d-8 ! Tolerance for final solution.
    integer :: step, i

    s = aa_state_init(segments)

    ! Initialise the (normalised) inner radius to the reference value (1.0).
    if (present(prev_state) .and. prev_state%n == segments) then
       s%rad = prev_state%rad
    else
       s%rad = 1.0_dp
    end if

    old_r0 = s%rad(1)
    delta = update_radius(s%rad)
    step = 1

    do while (delta > toln)
       old_r0 = s%rad(1)
       delta = update_radius(s%rad)
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
      real(kind=dp), intent(inout) :: radius(segments)
      real(kind=dp) :: new_radius(segments), delta, x, res_aff, pa_glom
      integer :: i
      logical :: sum_mask(segments)

      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         ! The absolute radius of the segment.
         s%abs_rad(i) = radius(i) * aa_seg_inner_rho(p, x)
         ! The hemodynamic resistance of the segment.
         s%res(i) = p%k / (s%abs_rad(i) ** 4_dp * segments)
      end do

      ! Calculate the net resistance of the afferent arteriole.
      res_aff = sum(s%res)

      ! Calculate the hydrostatic pressure at the mid-point of each segment.
      sum_mask(:) = .false.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         if (i > 1) then
            sum_mask(i - 1) = .true.
         end if
         ! The hydrostatic pressure at the mid-point of the segment.
         s%pa(i) = p%Pa - (p%Pa - p%Pv) &
              * (sum(s%res, sum_mask) + s%res(i) * 0.5_dp) / (res_aff + p%Re)
      end do

      ! Calculate the hydrostatic pressure in the glomerulus.
      pa_glom = p%Pv + (p%Pa - p%Pv) * p%Re / (res_aff + p%Re)

      ! Calculate the afferent blood flow (nl/s).
      s%aff_bf = (p%Pa - pa_glom) / sum(s%res) * 1d-6

      ! Calculate the combined activation of the myogenic and TGF responses.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         s%myo_act(i) = aa_myogenic_act(p, s%pa(i))
         if (present(ptgf)) then
            if (x >= (1_dp - 1_dp / p%alpha)) then
               s%tgf_act(i) = ptgf * (1_dp - p%alpha * (1_dp - x))
            else
               s%tgf_act(i) = 0_dp
            end if
         else
            s%tgf_act(i) = aa_tgf_act(p, pa_glom, x)
         end if
         s%net_act(i) = s%myo_act(i) + s%tgf_act(i)
      end do

      ! Finally, update the (normalised) radius of each segment.
      do i = 1, segments
         x = real(i, dp) / real(segments, dp)
         new_radius(i) = aa_seg_norm_radius(p, s%net_act(i), s%pa(i), x)
      end do

      delta = maxval(abs(radius - new_radius))
      radius = (1_dp - interp_frac) * abs(radius) &
           + interp_frac * abs(new_radius)

    end function update_radius

  end function aa_state_of

end module model_base
