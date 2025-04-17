!-------------------------------------------------------------------------------------------------
!	GGMCalc, A program to compute some components of gravity field using Global Geopotential Models
!
!	Copyright (C) 2010-2018 Siamak Moazezi
!	
!	GGMCalc is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 3 of the License, or
!	(at your option) any later version.
!
!	GGMCalc is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with GGMCalc.  If not, see <http://www.gnu.org/licenses/>.
!
!	Contact info: http://www.sourceforge.net/projects/xgravity
!	              Siamak Moazezi <s.moazezi@srbiau.ac.ir>
!-------------------------------------------------------------------------------------------------

module ggmcalc_mod
contains

   subroutine ggmcalc_main(coef_file, ellipsoid_config, input_file, ellipsoidal_true, T_B_IT_true, iteration, N_Precision, zeta_Precision, productionname)

     use nrtype
     use inout_mod
     use Legendre_mod
     use C_2n_mod
     use gamma_mod
     use W_mod
     use U_mod
     use undulation_mod
     use height_anomaly_mod
     use gravity_disturbance_mod
     use gravity_anomaly_mod
     use ProgressBar_mod
     use duration
     implicit none
     integer(i4b) :: nmin=0, nmax, nmax_p, n, m, ii, jj, nn
     real(longdp) :: GM_ellipsoid, a, b, f_reciprocal, omega, U_0, m_ellipsoid, gamma_0
     real(longdp) :: GM, R, c, s, sigma_c, sigma_s
     real(longdp), allocatable, dimension(:) :: C_bar_nm
     real(longdp), allocatable, dimension(:,:) :: p_bar_nm, dp_bar_nm
     integer(i4b) :: Iteration
     real(longdp) :: N_Precision, zeta_Precision, Latitude, Longitude, h, depth, ice, rho, N_O, N_E
     real(longdp), allocatable, dimension(:) :: Latitude_p, Latitude_p_geodetic, Longitude_p, h_p, Undulation_p, Undulation_Simple_Assumed_Density_p, Undulation_Assumed_Density_p, Gravity_Disturbance_p, Classical_Gravity_Anomaly_p, Classical_Gravity_Anomaly_2_p, Molodensky_Gravity_Anomaly_p
     real(longdp), allocatable, dimension(:) :: depth_p, ice_p
     character*(1) :: temp
     character(len=50) :: str1, str2, product_type, modelname, errors, norm, tide_system, coefficient_kind, ellipsoidname
     character(len=256) :: buffer, label, productionname, outputname, outputlogname
     integer(i4b) :: pos, Number_p
     logical :: ellipsoidal_true, T_B_IT_true, writing

     integer :: days1, hours1, minutes1, days2, hours2, minutes2, i, stdout = 6
     real(longdp) :: seconds1, seconds2
!      character(len=:), allocatable :: coef_file, ellipsoid_config, input_file
     character(len=512), intent(in) :: coef_file, ellipsoid_config, input_file
     character(len=500) :: number_to_string

    ! nproc = omp_get_num_procs()
    ! write(*,*) ' Number of available threads ',nproc
    ! write(*,*) " Insert number of threads "
    ! read(*,*) nproc
    ! write(*,*)' Threads used ', nproc
    ! call omp_set_num_threads(nproc)

    !!----------------------------------------------------------------------------------------------------
    !!---------------------------- Reading of Reference Ellipsoid Parameters -----------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)', advance='no') 'Read Ellipsoid Parameters ... '

     open(12, file=ellipsoid_config)
    !!----------------------------------------------------------------------------------------------------
     read(12, *) GM_ellipsoid
     read(12, *) a
     read(12, *) b
     read(12, *) f_reciprocal
     read(12, *) omega
     read(12, *) U_0
     read(12, *) ellipsoidname
    !!----------------------------------------------------------------------------------------------------
     close(12)
    !!----------------------------------------------------------------------------------------------------
     m_ellipsoid = (omega**2.0_longdp * a**2.0_longdp *b) / GM_ellipsoid
     gamma_0 = GM_ellipsoid / a**2.0_longdp
    !!----------------------------------------------------------------------------------------------------
    !!------------------------- End of Reading of Reference Ellipsoid Parameters -------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)') 'done!'
     write(stdout, *) 'Ellipsoid Parameters:'
     write(stdout, *) '----------------------------------------------------------------------'
     write(stdout, *) 'Reference Ellipsoid Name: ', trim(adjustl(ellipsoidname))
     write(number_to_string, '(ES50.6)') GM_ellipsoid
     write(stdout, *) 'Ellipsoid Gravity Constant: ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.3)') a
     write(stdout, *) 'Semimajor Axis of the Ellipsoid: ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.3)') b
     write(stdout, *) 'Semiminor Axis of the Ellipsoid: ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.9)') f_reciprocal
     write(stdout, *) 'Reciprocal Flattening: ', trim(adjustl(number_to_string))
     write(number_to_string, '(ES50.7)') omega
     write(stdout, *) 'Angular Velocity of the Earth: ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.3)') U_0
     write(stdout, *) 'Normal Potential at the Ellipsoid: ', trim(adjustl(number_to_string))
     write(number_to_string, '(ES50.5)') m_ellipsoid
     write(stdout, *) 'm = ω**2 a**2 b/(GM): ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.10)') gamma_0
     write(stdout, *) 'Normal Gravity at the Equator: ', trim(adjustl(number_to_string))
     write(stdout, *) '----------------------------------------------------------------------'

    !!----------------------------------------------------------------------------------------------------
    !!----------------------------- Reading of Geopotential Model Parameters -----------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)', advance='no') 'Read Geopotential Model Parameters ... '

     product_type = 'Unknown'
     modelname = 'Unknown'
     errors = 'Unknown'
     norm = 'Unknown'
     tide_system = 'Unknown'
     open(12, file=coef_file)

     do
      read(12, '(A)') buffer
      pos = scan(buffer, ' ')
      if (pos /= 0) then
              label = buffer(1:pos)
              buffer = adjustl(buffer(pos+1:))
              select case (label)
              case ('product_type')
                      read(buffer, *) str2
                      product_type = trim(str2)
              case ('modelname')
                      read(buffer, *) str2
                      modelname = trim(str2)
              case ('earth_gravity_constant')
                      read(buffer, *) str2
                      read(str2, *) GM
              case ('radius')
                      read(buffer, *) str2
                      read(str2, *) R
              case ('max_degree')
                      read(buffer, *) str2
                      read(str2, *) nmax
              case ('errors')
                      read(buffer, *) str2
                      errors = trim(str2)
              case ('norm')
                      read(buffer, *) str2
                      norm = trim(str2)
              case ('tide_system')
                      read(buffer, *) str2
                      tide_system = trim(str2)
              !key    L    M         C                  S           sigma C    sigma S
              case ('end_of_head')
                      exit
              end select
      end if
     end do

    !!----------------------------------------------------------------------------------------------------
     allocate(C_bar_nm(nmax*(nmax+2)+1))
     do ii = 0, nmax
      do jj = 0, ii !nmax
              if (jj == 0) then
                      C_bar_nm((nmax - jj) * (nmax - jj + 1) + nmax - ii + 1) = 0.0_longdp
              else
                      nn = (nmax - jj) * (nmax - jj + 1) + 2*(nmax - ii) + 1
                      C_bar_nm(nn) = 0.0_longdp
                      C_bar_nm(nn+1) = 0.0_longdp
              end if
      end do
     end do
    !!----------------------------------------------------------------------------------------------------
 100     continue
     if (index(errors, 'no') /= 0) then
      read(12, *, end=200) coefficient_kind, n, m, c, s
     else
      read(12, *, end=200) coefficient_kind, n, m, c, s, sigma_c, sigma_s
     end if

     if (n < nmin .or. n > nmax) goto 100
     if (index(coefficient_kind, 'gfc') == 0) goto 100

     if (m == 0) then
      C_bar_nm((nmax - m) * (nmax - m + 1) + nmax - n + 1) = c
     else
      nn = (nmax - m) * (nmax - m + 1) + 2*(nmax - n) + 1
      C_bar_nm(nn) = c
      C_bar_nm(nn+1) = s
     end if
     goto 100
 200     continue
    !!----------------------------------------------------------------------------------------------------
     close(12)
    !!----------------------------------------------------------------------------------------------------
    !!-------------------------- End of Reading of Geopotential Model Parameters -------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)') 'done! '

     write(stdout, *) 'Geopotential Model Parameters:'
     write(stdout, *) '----------------------------------------------------------------------'
     write(stdout, *) 'Product Type: ', trim(adjustl(product_type))
     write(stdout, *) 'Model Name: ', trim(adjustl(modelname))
     write(number_to_string, '(ES50.10)') GM
     write(stdout, *) 'Earth Gravity Constant: ', trim(adjustl(number_to_string))
     write(number_to_string, '(F50.3)') R
     write(stdout, *) 'Radius of the Earth: ', trim(adjustl(number_to_string))
     write(number_to_string, *) nmax
     write(stdout, *) 'Maximum Degree of Model: ', trim(adjustl(number_to_string))
     write(stdout, *) 'Errors Type of Model: ', trim(adjustl(errors))
     write(stdout, *) 'Normalization Type: ', trim(adjustl(norm))
     write(stdout, *) 'Tide System: ', trim(adjustl(tide_system))
     write(stdout, *) '----------------------------------------------------------------------'

     if (ellipsoidal_true) then
      write(stdout, '(a)') 'Ellipsoidal height expected!'
     else
      write(stdout, '(a)') 'Orthometric height expected!'
     end if

     if (T_B_IT_true) then
      write(stdout, '(a)') 'Topography, Bathymetry, and Ice Thickness is Available Separately!'

      allocate(depth_p(nmax_p))
      allocate(ice_p(nmax_p))

      do ii = 1, nmax_p
              depth_p(ii) = 0.0_longdp
              ice_p(ii) = 0.0_longdp
      end do

     end if

     if (Iteration == 0) then
       write(number_to_string, *) Iteration
       write(stdout, *) 'Number of Iteration, s: ', trim(adjustl(number_to_string))
       write(number_to_string, '(F50.4)') N_Precision
       write(stdout, *) 'Undulation Precision, m: ', trim(adjustl(number_to_string))
       write(number_to_string, '(F50.4)') zeta_Precision
       write(stdout, *)'Height Anomaly Precision, m: ', trim(adjustl(number_to_string))
     else
       N_Precision = 0.0_longdp
       Zeta_Precision = 0.0_longdp
     end if

    !!----------------------------------------------------------------------------------------------------
    !!---------------------------------- Reading of Points Coordinates -----------------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)', advance='no') 'Read Data Points ... '

     open(12, file=input_file)

     nmax_p = 0
 301     continue
     nmax_p = nmax_p + 1
     read(12, *, end=401) temp
     goto 301
 401     continue

     allocate(latitude_p(nmax_p))
     allocate(latitude_p_geodetic(nmax_p))
     allocate(longitude_p(nmax_p))
     allocate(h_p(nmax_p))
     allocate(Undulation_p(nmax_p))
     allocate(Undulation_Simple_Assumed_Density_p(nmax_p))
     allocate(Undulation_Assumed_Density_p(nmax_p))
     allocate(Height_Anomaly_p(nmax_p))
     allocate(Gravity_Disturbance_p(nmax_p))
     allocate(Classical_Gravity_Anomaly_p(nmax_p))
     allocate(Classical_Gravity_Anomaly_2_p(nmax_p))
     allocate(Molodensky_Gravity_Anomaly_p(nmax_p))
    !!----------------------------------------------------------------------------------------------------
     do ii = 1, nmax_p
      latitude_p(ii) = 0._longdp
      latitude_p_geodetic(ii) = 0._longdp
      longitude_p(ii) = 0._longdp
      h_p(ii) = 0.0_longdp
      Undulation_p(ii) = 0.0_longdp
      Undulation_Simple_Assumed_Density_p(ii) = 0.0_longdp
      Undulation_Assumed_Density_p(ii) = 0.0_longdp
      Height_Anomaly_p(ii) = 0.0_longdp
      Gravity_Disturbance_p(ii) = 0.0_longdp
      Classical_Gravity_Anomaly_p(ii) = 0.0_longdp
      Classical_Gravity_Anomaly_2_p(ii) = 0.0_longdp
      Molodensky_Gravity_Anomaly_p(ii) = 0.0_longdp
     end do
    !!----------------------------------------------------------------------------------------------------

     close(12)
     open(12, file=input_file)

     ii = 0
 300     continue
     if (T_B_IT_true) then
      !------------============= Reading Topography, Bathymetry, and Ice Thickness =============------------
      read(12, *, end=400) latitude, longitude, h, depth, ice
      !------------============= Reading Topography, Bathymetry, and Ice Thickness =============------------
     else
      read(12, *, end=400) latitude, longitude, h
     end if

     ii = ii + 1

     Latitude_p_geodetic(ii) = Latitude
     Latitude_p(ii) = Latitude
     Longitude_p(ii) = Longitude
     h_p(ii) = h

     if (T_B_IT_true) then
      depth_p(ii) = depth
      ice_p(ii) = ice
     end if

     goto 300
 400     continue
     nmax_p = ii
     ii = 0
    !!write(stdout, *) nmax_p
    !!pause
    !!----------------------------------------------------------------------------------------------------
     close(12)
    !!----------------------------------------------------------------------------------------------------
    !!------------------------------- End of Reading of Points Coordinates -------------------------------
    !!----------------------------------------------------------------------------------------------------

     write(stdout, '(a)') 'done! '

     write(number_to_string, *) nmax_p
     write(stdout, *) 'Number of Data Points: ', trim(adjustl(number_to_string))

     allocate(p_bar_nm(nmax+1, nmax+1))
     do ii = 0, nmax
      do jj = 0, nmax
              p_bar_nm(ii+1, jj+1) = 0.0_longdp
      end do
     end do

     allocate(dp_bar_nm(nmax+1, nmax+1))
     do ii = 0, nmax
      do jj = 0, nmax
              dp_bar_nm(ii+1, jj+1) = 0.0_longdp
      end do
     end do

     write(stdout, *) 'Production Name: ', trim(adjustl(productionname))
     write(stdout, '(a)') 'Computing ... '
     write(stdout, *) ''

     Number_p = 0
     call progress(Number_p, nmax_p) ! generate the progress bar.

     writing = .false.

     call exec_time(days1, hours1, minutes1, seconds1)

    !!$OMP PARALLEL num_threads(200)
    !!$OMP DO private(ii)
     do ii = 1, nmax_p
    ! if (writing .eqv. .false.) then
    ! 	writing = .true.
              call progress(Number_p, nmax_p) ! generate the progress bar.
    ! 	writing = .false.
    ! end if

      call Legendre(p_bar_nm, dp_bar_nm, 10, h_p(ii), Latitude_p(ii), longitude_p(ii), a, b)

      if (h_p(ii) >= 0.0_longdp) then
              rho = 2670.0_longdp
      else
              rho = 1030.0_longdp
      end if

      Undulation_p(ii) = undulation(latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm, Iteration, N_Precision)

      if (ellipsoidal_true) then
              N_O = 0.0_longdp
              if (T_B_IT_true) then
                N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - Undulation_p(ii))**2.0_longdp - 2670.0_longdp * (depth_p(ii) + Undulation_p(ii))**2.0_longdp + 1030.0_longdp * (depth_p(ii) + Undulation_p(ii))**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
                N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - N_E)**2.0_longdp - 2670.0_longdp * (depth_p(ii) + N_E)**2.0_longdp + 1030.0_longdp * (depth_p(ii) + N_E)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
              else
                N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - Undulation_p(ii))**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
                N_E = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - N_E)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
              end if
      else
              if (T_B_IT_true) then
                N_O = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * h_p(ii)**2.0_longdp - 2670.0_longdp * depth_p(ii)**2.0_longdp + 1030.0_longdp * depth_p(ii)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
              else
                N_O = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * h_p(ii)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
              end if
              N_E = 0.0_longdp
      end if

      Undulation_Simple_Assumed_Density_p(ii) = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * rho * (h_p(ii) - N_E)**2.0_longdp) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
      Height_Anomaly_p(ii) = height_anomaly(h_p(ii) + N_O, latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm, Iteration, zeta_Precision)
      Gravity_Disturbance_p(ii) = gravity_disturbance(h_p(ii) + N_O, latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp
      Classical_Gravity_Anomaly_p(ii) = ga_classic_first_approximation(Undulation_Simple_Assumed_Density_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100
      Classical_Gravity_Anomaly_2_p(ii) = ga_classic_second_approximation(h_p(ii) + N_O, Undulation_Simple_Assumed_Density_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp
      Molodensky_Gravity_Anomaly_p(ii) = gravity_anomaly_molodensky(h_p(ii) + N_O, Height_Anomaly_p(ii), latitude_p(ii), longitude_p(ii), GM, R, GM_ellipsoid, a, b, f_reciprocal, omega, C_bar_nm, nmax, p_bar_nm, dp_bar_nm) * 100000.0_longdp

      if (T_B_IT_true) then
        Undulation_Assumed_Density_p(ii) = Undulation_p(ii) - ((2.0_longdp*PI_D * 6.67428E-11 * (2670.0_longdp * (h_p(ii) - N_E)**2.0_longdp - 2670.0_longdp * (depth_p(ii) + N_E)**2.0_longdp + 1030.0_longdp * (depth_p(ii) + N_E)**2.0_longdp - 2670.0_longdp * ice_p(ii)**2.0_longdp + 900.0_longdp * ice_p(ii)**2.0_longdp)) / gamma_h(latitude_p(ii), 0.0_longdp, GM_ellipsoid, a, b, f_reciprocal, omega))
      end if

      Number_p = Number_p + 1
     end do
    !!$OMP END DO NOWAIT
    !!$OMP END PARALLEL

     call exec_time(days1, hours1, minutes1, seconds1)
     call exec_time(days2, hours2, minutes2, seconds2)

     call progress_close

     write(stdout, '(a)', advance='no') 'Saving Results ... '

    write(outputname, *) 'output_', trim(modelname), '_', trim(productionname), '.dat'
    write(outputlogname, *) 'output_', trim(modelname), '_', trim(productionname), '.log'
     open(12, file=trim(adjustl(outputname)))

     if (T_B_IT_true) then
        write(12,*) 'Latitude(degree) Longitude(degree) Height(m) Water_Depth(m) Ice_Thickness(m) Undulation(m) Undulation_with_Simple_Assumed_Density(m) Undulation_with_Assumed_Density(m) Heigh_Anomaly(m) Gravity_Disturbance(mGal) Classical_Gravity_Anomaly(mGal) Classical_Gravity_Anomaly_2(mGal) Molodensky_Gravity_Anomaly(mGal)'
     else
	write(12,*) 'Latitude(degree) Longitude(degree) Height(m) Undulation(m) Undulation_with_Simple_Assumed_Density(m) Height_Anomaly(m) Gravity_Disturbance(mGal) Classical_Gravity_Anomaly(mGal) Classical_Gravity_Anomaly_2(mGal) Molodensky_Gravity_Anomaly(mGal)'
     endif

     do ii = 1, nmax_p
    ! call progress(ii, nmax_p) ! generate the progress bar.

      if (T_B_IT_true) then
        write(12, *) latitude_p_geodetic(ii), longitude_p(ii), h_p(ii), depth_p(ii), ice_p(ii), Undulation_p(ii), Undulation_Simple_Assumed_Density_p(ii), Undulation_Assumed_Density_p(ii), Height_Anomaly_p(ii), Gravity_Disturbance_p(ii), Classical_Gravity_Anomaly_p(ii), Classical_Gravity_Anomaly_2_p(ii), Molodensky_Gravity_Anomaly_p(ii)
      else
        write(12, *) latitude_p_geodetic(ii), longitude_p(ii), h_p(ii), Undulation_p(ii), Undulation_Simple_Assumed_Density_p(ii), Height_Anomaly_p(ii), Gravity_Disturbance_p(ii), Classical_Gravity_Anomaly_p(ii), Classical_Gravity_Anomaly_2_p(ii), Molodensky_Gravity_Anomaly_p(ii)
      end if
     end do

    ! call progress_close

     close(12)

     write(stdout, '(a)') 'done!'

     write(stdout, 4141) days1, hours1, minutes1, seconds1
 4141     format(' Computing time: ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')
     call exec_time(days2, hours2, minutes2, seconds2)
     write(stdout, 5151) days2, hours2, minutes2, seconds2
 5151     format(' Saving time:    ', I3, ' days ', I3, ' hours ', I3, ' minutes ', F7.3, ' seconds')

     open(12, file=trim(adjustl(outputlogname)))

     close(12)

     deallocate(C_bar_nm)
     deallocate(latitude_p)
     deallocate(latitude_p_geodetic)
     deallocate(longitude_p)
     deallocate(h_p)
     deallocate(Undulation_p)
     deallocate(Undulation_Simple_Assumed_Density_p)
     deallocate(Undulation_Assumed_Density_p)
!      deallocate(Height_Anomaly_p)
     deallocate(Classical_Gravity_Anomaly_p)
     deallocate(Classical_Gravity_Anomaly_2_p)
     deallocate(Molodensky_Gravity_Anomaly_p)
     deallocate(p_bar_nm)
     deallocate(dp_bar_nm)

     if (T_B_IT_true) then
      deallocate(depth_p)
      deallocate(ice_p)
     end if

   end subroutine ggmcalc_main
end module ggmcalc_mod