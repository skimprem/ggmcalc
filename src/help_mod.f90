module help_mod

      implicit none

    contains

      function print_help() result(out)
      character(len=:), allocatable :: out

      out = 'Usage: ggmcalc [OPTION]... [FILE]...'
      out = out // char(10) // ''
      out = out // char(10) // 'Reads EGM coefficients and point coordinates from specified FILEs,'
      out = out // char(10) // 'then calculates height anomalies, geoid undulations, and other values.'
      out = out // char(10) // ''
      out = out // char(10) // 'Options:'
      out = out // char(10) // '  --coeffs        Specify the EGM coefficients file in ICGEM format.'
      out = out // char(10) // '  --ellip_conf    Specify the ellipsoid configuration file with the following format:'
      out = out // char(10) // ''
      out = out // char(10) // '                  3986004.418d8'
      out = out // char(10) // '                  6378137.0d0'
      out = out // char(10) // '                  6356752.3142d0'
      out = out // char(10) // '                  298.257223563'
      out = out // char(10) // '                  7292115.0d-11'
      out = out // char(10) // '                  62636851.7146d0'
      out = out // char(10) // '                  WGS_84'
      out = out // char(10) // ''
      out = out // char(10) // '  --input_file    Specify the input file(s) with coordinates and other required fields.'
      out = out // char(10) // ''
      out = out // char(10) // '  --ortho_height  [optional] To set the input height in orthometric system.'
      out = out // char(10) // '                  Default is ELLIPSOIDAL.'
      out = out // char(10) // ''
      out = out // char(10) // '  --topo          [optional] To set topography, bathymetry and Ice Thickness as'
      out = out // char(10) // '                  separately. Default is FALSE'
      out = out // char(10) // ''
      out = out // char(10) // '  --iters         [optional] To set iteration number. Default is 0.'
      out = out // char(10) // ''
      out = out // char(10) // '  --n_prec        [optional] To set geoid undulations precision. Default is 0.001.'
      out = out // char(10) // ''
      out = out // char(10) // '  --z_prec        [optional] To set height anomaly precision. Default is 0.001.'
      out = out // char(10) // ''
      out = out // char(10) // '  --output        Specify the suffix name for the output CSV file with results.'
      out = out // char(10) // ''
      out = out // char(10) // '  --help          Display this help message and exit.'

      end function print_help

end module help_mod
