module config_param

  ! Ludovico Fossa 2026
  ! performs zeroth-order interpolation of a three-dimensional cartesian field

  ! domain in
  integer, parameter :: nxt_in = 256, nyt_in = 256, nzt_in = 512
  ! domain out
  integer, parameter :: nxt_out = 512, nyt_out = 512, nzt_out = 1024
  ! number of components and bytes per scalar
  integer, parameter :: n_comp = 3, n_bytes = 8
  ! I/O directory
  character(len=120) :: datadir = '_'
  character(len=120) :: prefix = '_', suffix = '_'
  integer, parameter :: iStep = 0

end module config_param