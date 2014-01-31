!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

!> Parser for command-line option.
module cmdline

  private

  integer, parameter :: ARG_MAX = 32

  type options
     character(len=ARG_MAX) :: command
     logical :: verbose
     logical :: glomerulus
     logical :: blood_flow
  end type options

  public :: ARG_MAX, options, parse_cmdline, print_help

contains

  function parse_cmdline() result(opts)
    type(options) :: opts
    integer :: i, status
    character(len=ARG_MAX) :: arg

    ! Set the default options.
    opts%command = ""
    opts%verbose = .false.
    opts%glomerulus = .false.
    opts%blood_flow = .false.

    do i = 1, command_argument_count()
       call get_command_argument(i, arg, status=status)
       if (status < 0) then
          ! The argument was truncated to fit into the buffer.
          write (*,"(A)") "ERROR: truncated command-line argument."
          call exit(1)
       elseif (status > 0) then
          ! The argument retrieval failed.
          write (*,"(A)") "ERROR: unable to retrieve command-line argument."
          call exit(1)
       end if
       call update_opts(opts, arg)
    end do

  end function parse_cmdline

  subroutine print_help()
    write (*,"(A)") ""
    write (*,"(A)") "USAGE"
    write (*,"(A)") "      run_model COMMAND [OPTIONS]"
    write (*,"(A)") ""
    write (*,"(A)") "COMMANDS"
    write (*,"(A)") "      test    Run model test cases."
    write (*,"(A)") "      sngfr   Generate a lookup table for model SNGFR."
    write (*,"(A)") "              This can take a very long time."
    write (*,"(A)") ""
    write (*,"(A)") "OPTIONS"
    write (*,"(A)") "      --help  Display help about commands and options."
    write (*,"(A)") ""
    write (*,"(A)") "      --verbose"
    write (*,"(A)") "              Write additional output."
    write (*,"(A)") ""
    write (*,"(A)") "  Testing Options"
    write (*,"(A)") "      --blood-flow"
    write (*,"(A)") "              Calculate afferent blood flow."
    write (*,"(A)") "              This can take several minutes to complete."
    write (*,"(A)") "                    Default: false."
    write (*,"(A)") ""
    write (*,"(A)") "  Lookup Table Options"
    write (*,"(A)") "      --glom  Use the explicit glomerulus model."
    write (*,"(A)") "                    Default: false."
    write (*,"(A)") ""
  end subroutine print_help

  subroutine update_opts(opts, arg)
    type(options), intent(inout) :: opts
    character(len=ARG_MAX), intent(in) :: arg

    select case (trim(arg))
       case ("test", "sngfr")
          opts%command = arg
       case ("--help")
          call print_help()
          call exit(0)
       case ("--verbose")
          opts%verbose = .true.
       case ("--glom")
          opts%glomerulus = .true.
       case ("--blood-flow")
          opts%blood_flow = .true.
       case default
          write (*,"(3A)") "ERROR: invalid argument '", trim(arg), "'."
          call print_help()
          call exit(0)
    end select

  end subroutine update_opts

end module cmdline
