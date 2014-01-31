!
! Copyright 2014, Robert Moss.
! Distributed under the BSD 3-Clause license (see LICENSE).
!

program main

  use cmdline, only: ARG_MAX, options, parse_cmdline, print_help
  use tests, only: run_tests
  use sngfr, only: sngfr_table

  type(options) :: opts

  opts = parse_cmdline()

  select case (trim(opts%command))
     case ("test")
        call run_tests(opts)
     case ("sngfr")
        call sngfr_table(opts)
     case default
        call print_help()
        call exit(0)
  end select

end program main
