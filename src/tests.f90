!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

!> A set of test cases for the afferent arteriolar model.
module tests

  use model_base, only: dp

  implicit none

  private

  public run_tests

  real(kind=dp), parameter, private :: aa_interp = 0.2_dp

contains

  !> Reproduce the data from Figure 1.
  subroutine test_stress_strain(p, unit)

    use model_base

    type(aa_params), intent(in) :: p
    integer, intent(in) :: unit

    integer :: pstep, estep
    real(kind=dp) :: psi, epsilon, strain_p, strain_a

    open(unit, file="stress_strain.ssv")

    do pstep = 1, 5
       psi = real(pstep, dp) / 5_dp

       do estep = -100, 80
          epsilon = real(estep, dp) / 100_dp
          strain_p = aa_seg_passive_strain(p, epsilon)
          strain_a = aa_seg_active_strain(p, epsilon, psi)

          write (unit, "(F3.1, F7.2, F7.1, F7.1)") &
               psi, epsilon, strain_p + strain_a, strain_a
       end do
    end do
    close(unit)
  end subroutine test_stress_strain


  ! Reproduce the data from Figure 2.
  subroutine test_psi_pressure(p, unit)

    use model_base

    type(aa_params), intent(in) :: p
    integer, intent(in) :: unit

    integer :: pstep, rstep
    real(kind=dp) :: psi, radius, pressure

    open(unit, file="pressure.ssv")

    do pstep = 1, 6
       psi = real(pstep, dp) / 5_dp

       ! When radius = 0, epsilon = -1 and the denominator becomes 0.
       do rstep = 1, 180
          radius = real(rstep, dp) / 100_dp
          pressure = aa_seg_transmural(p, psi, radius)

          write (unit, "(F3.1, F6.2, F12.1)") psi, radius, pressure
       end do
    end do
    close(unit)
  end subroutine test_psi_pressure


  ! Reproduce the data from Figure 4.
  subroutine test_pressure_radius(p0, unit)

    use model_base

    type(aa_params), intent(in) :: p0
    integer, intent(in) :: unit

    type(aa_params) :: p
    integer :: pstep
    real(kind=dp) :: pressure, psi, radius

    p = p0

    open(unit, file="radius.ssv")
    do pstep = 1, 300
       pressure = real(pstep, dp) / 10_dp
       psi = aa_myogenic_act(p, pressure)
       radius = aa_seg_norm_radius(p, psi, pressure)
       write (unit, "(F6.2, F7.3, F5.1)") pressure, radius, 1.0_dp
       call flush(unit)
    end do

    p%Pmin = 4.9_dp
    p%pMax = 12.3_dp
    do pstep = 1, 300
       pressure = real(pstep, dp) / 10_dp
       psi = aa_myogenic_act(p, pressure)
       radius = aa_seg_norm_radius(p, psi, pressure)
       write (unit, "(F6.2, F7.3, F5.1)") pressure, radius, 2.0_dp
       call flush(unit)
    end do

    p%Pmin = -6.2_dp
    p%pMax = 23.4_dp
    do pstep = 1, 300
       pressure = real(pstep, dp) / 10_dp
       psi = aa_myogenic_act(p, pressure)
       radius = aa_seg_norm_radius(p, psi, pressure)
       write (unit, "(F6.2, F7.3, F5.1)") pressure, radius, 0.5_dp
       call flush(unit)
    end do

    close(unit)
  end subroutine test_pressure_radius


  ! Reproduce the data from Figure 6.
  subroutine test_afferent_bf(unit, verbose)

    use model_base

    integer, intent(in) :: unit
    logical, intent(in) :: verbose

    integer, parameter :: nslope = 3
    real(kind=dp), parameter :: slope(nslope) = (/ 0.5_dp, 1.0_dp, 2.0_dp /)
    integer, parameter :: ns = 7, nsegs(ns) = (/ 100, 50, 25, 10, 5, 4, 3 /)
    type(aa_params) :: p
    type(aa_state) :: states(size(nsegs))
    integer :: i, j, s

    open(unit, file="aa_bf.ssv")
    write (unit, "(A4, 3A8)") "pa", "Q", "n", "slope"

    do s = 1, nslope
       states = aa_state_init(nsegs)
       p = aa_params_init(shape=1, myo_slope = slope(s))

       do i = 50, 220, 5
          p%Pa = real(i, dp) / 10_dp

          if (verbose) write (*, "(F4.1, A4)", advance="no") p%Pa, "kPa"

          do j = 1, size(states)
             states(j) = aa_state_of(p, aa_interp, states(j)%n, states(j))
             write (unit, "(F4.1, F8.3, I8, F8.3)") p%Pa, states(j)%aff_bf, &
                  states(j)%n, slope(s)
             if (verbose) write (*, "(4F7.3)", advance="no") states(j)%aff_bf
          end do
          if (verbose) write (*, "(A5)") "nl/s"
       end do

    end do

    close(unit)
  end subroutine test_afferent_bf


  subroutine write_aa_state(p, unit, verbose, ptgf)

    use model_base

    type(aa_params), intent(in) :: p
    integer, intent(in) :: unit
    logical, intent(in) :: verbose
    real(kind=dp), intent(in), optional :: ptgf

    integer :: i, n
    real(kind=dp) :: x
    type(aa_state) :: s

    n = 100
    if (present(ptgf)) then
       s = aa_state_of(p, aa_interp, n, ptgf=ptgf)
    else
       s = aa_state_of(p, aa_interp, n)
    end if

    if (verbose) then
       write (*, "(F4.1, A4, F8.3, A5)") p%Pa, "kPa", s%aff_bf, "nL/s"
    end if

    do i = 1, n
       x = real(i, dp) / real(n, dp)
       write (unit, "(F8.4, 7ES15.6, I2, F5.1)") x, s%rad(i), &
            s%abs_rad(i), s%res(i), s%pa(i), s%net_act(i), &
            s%myo_act(i), s%tgf_act(i), p%shape, p%Pa
    end do

  end subroutine write_aa_state


  subroutine write_aa_state_glom(p, unit, verbose, ptgf)

    use model_glom

    type(glom_params), intent(in) :: p
    integer, intent(in) :: unit
    logical, intent(in) :: verbose
    real(kind=dp), intent(in) :: ptgf

    integer :: i, n
    real(kind=dp) :: x
    type(glom_state) :: s

    n = 100
    s = glom_state_of(p, aa_interp, n, ptgf)

    if (verbose) then
       write (*, "(F4.1, A4, F8.3, A5)") p%aa%Pa, "kPa", s%aff_bf, "nL/s"
    end if

    do i = 1, n
       x = real(i, dp) / real(n, dp)
       write (unit, "(F8.4, 7ES15.6, I2, F5.1, F10.4)") x, s%aa%rad(i), &
            s%aa%abs_rad(i), s%aa%res(i), s%aa%pa(i), s%aa%net_act(i), &
            s%aa%myo_act(i), s%aa%tgf_act(i), p%aa%shape, p%aa%Pa, s%sngfr
    end do

  end subroutine write_aa_state_glom


  ! Reproduce the data from Figures 7 and 8.
  subroutine test_aa_profiles(unit, verbose)

    use model_base

    integer, intent(in) :: unit
    logical, intent(in) :: verbose

    type(aa_params) :: p
    real(kind=dp), parameter :: ptgf = 0.06_dp
    integer :: i

    open(unit, file="profiles.ssv")

    p = aa_params_init(shape=1, myo_slope=1.0_dp)

    ! Iterate from 50 to 150 mmHg in 25 mmHg increments.
    do i = 2, 6
       p%Pa = real(10 * i, dp) / 3.0_dp
       call write_aa_state(p, unit, verbose, ptgf=ptgf)
    end do

    close(unit)
  end subroutine test_aa_profiles


  subroutine test_aa_profiles_glom(unit, verbose)

    use model_glom

    integer, intent(in) :: unit
    logical, intent(in) :: verbose

    type(glom_params) :: p
    real(kind=dp), parameter :: ptgf = 0.06_dp
    integer :: i

    open(unit, file="profiles_glom.ssv")

    p = glom_params_init(shape=1, myo_slope=2.0_dp)

    ! Iterate from 50 to 150 mmHg in 25 mmHg increments.
    do i = 2, 6
       p%aa%Pa = real(10 * i, dp) / 3.0_dp
       call write_aa_state_glom(p, unit, verbose, ptgf=ptgf)
    end do

    close(unit)
  end subroutine test_aa_profiles_glom


  !> Run all of the tests.
  subroutine run_tests(opts, unit)

    use cmdline, only: options
    use model_base, only: aa_params, aa_params_init

    type(options) :: opts
    integer, intent(in), optional :: unit
    integer :: u
    type(aa_params) :: p

    if (present(unit)) then
       u = unit
    else
       u = 42
    end if

    p = aa_params_init()

    if (opts%verbose) write (*,"(A)") "Running tests ..."

    if (opts%verbose) write (*,"(A)", advance='no') "stress_strain ..."
    call test_stress_strain(p, u)
    if (opts%verbose) write (*,"(A)") " done"

    if (opts%verbose) write (*,"(A)", advance='no') "psi_pressure ..."
    call test_psi_pressure(p, u)
    if (opts%verbose) write (*,"(A)") " done"

    if (opts%verbose) write (*,"(A)", advance='no') "pressure_radius ..."
    call test_pressure_radius(p, u)
    if (opts%verbose) write (*,"(A)") " done"

    if (opts%blood_flow) then
       if (opts%verbose) write (*,"(A)", advance='no') "afferent_bf ..."
       call test_afferent_bf(u, opts%verbose)
       if (opts%verbose) write (*,"(A)") " done"
    end if

    if (opts%verbose) write (*,"(A)", advance='no') "aa_profiles ..."
    call test_aa_profiles(u, opts%verbose)
    call test_aa_profiles_glom(u, opts%verbose)
    if (opts%verbose) write (*,"(A)") " done"

    if (opts%verbose) write (*,"(A)") "Run 'plot.R' to plot the test results."

  end subroutine run_tests

end module tests
