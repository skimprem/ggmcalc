program ggmcalc

  use ggmcalc_mod
  use color_mod
  use help_mod
  use nrtype

  implicit none

  character(len=500) :: arg, number_to_string
  character(len=:), allocatable :: coef_file, ellipsoid_config, input_file
  integer :: stdout = 6, i
  character*(256) productionname
  real(longdp) :: N_Precision, zeta_Precision
  logical :: ellipsoidal_true, T_B_IT_true, writing
  integer(i4b) :: Iteration

  coef_file = ''
  ellipsoid_config = ''
  input_file = ''
  ellipsoidal_true = .true.
  T_B_IT_true = .false.
  Iteration = 1000
  N_Precision = 1e-3_longdp
  Zeta_Precision = 1e-3_longdp
  productionname = 'default'

  i = 0
  do while(i < command_argument_count())
    i = i + 1
    call get_command_argument(i, arg)
    select case(trim(adjustl(arg)))
    case('--coeffs')
      i = i + 1
      call get_command_argument(i, arg)
      coef_file = trim(adjustl(arg))
    case('--ellip_conf')
      i = i + 1
      call get_command_argument(i, arg)
      ellipsoid_config = trim(adjustl(arg))
    case('--input_file')
      i = i + 1
      call get_command_argument(i, arg)
      input_file = trim(adjustl(arg))
    case('--ortho_height')
      ellipsoidal_true = .false.
    case('--topo')
        T_B_IT_true = .true.
    case('--iters')
      i = i + 1
      call get_command_argument(i, arg)
      read(arg, *) Iteration 
    case('--n_prec')
      i = i + 1
      call get_command_argument(i, arg)
      read(arg, *) N_Precision
    case('--z_prec')
      i = i + 1
      call get_command_argument(i, arg)
      read(arg, *) Zeta_Precision
    case('--output')
      i = i + 1
      call get_command_argument(i, arg)
      productionname = trim(adjustl(arg))
    case('--help')
      write(stdout, '(a)') 'Help launch!'
      write(stdout, '(a)') print_help()
      stop
    end select
  end do

  if (len(trim(coef_file)) == 0) then
    write(stdout, *) 'Coeff file not set!'
    write(stdout, *) print_help()
    stop
  end if

  if (len(trim(ellipsoid_config)) == 0) then
    write(stdout, *) 'Ellipsoidal configuration file not set!'
    write(stdout, *) print_help()
    stop
  end if
  
  if (len(trim(input_file)) == 0) then
    write(stdout, *) 'Input file not set!'
    write(stdout, *) print_help()
    stop
  end if

  write(stdout, '(2a)') 'coeff file: ', color(coef_file, c_green)
  write(stdout, '(2a)') 'ellipsoid config: ', color(ellipsoid_config, c_green)
  write(stdout, '(2a)') 'input file: ', color(input_file, c_green)

  call ggmcalc_main(coef_file, ellipsoid_config, input_file, ellipsoidal_true, T_B_IT_true, iteration, N_Precision, zeta_Precision, productionname)

end program ggmcalc