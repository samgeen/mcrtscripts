! Trying to make this compile with f2py

subroutine solve_cooling_ramses(nH,T2,zsolar,dt,gammain,ncell,deltaT2out)
    !=======================================================================
        use cooling_module
        implicit none
      ! BRIDGE FUNCTION WITH SAME INTERFACE AS SOLVE_COOLING 
      ! Input/output variables to this function
      ! nH - hydrogen number density in PHYSICAL units
      ! T2 - temperature / mu in PHYSICAL units
      ! zsolar - Metallicity in solar units (Zphys / 0.02)
      ! dt - cooling timestep in seconds
      ! ncell - number of elements in the vector
      ! deltaT2 - temperature change in K/mu (??)
      integer::ncell
      real(kind=8),intent(in)::dt,gammain
      !real(kind=8),dimension(1:ncell),intent(in)::nH,T2,zsolar
      real(kind=8),dimension(1:ncell)::nH,T2,zsolar,boost,deltaT2
      real(kind=8),dimension(1:ncell),intent(out)::deltaT2out


      ! Set up bridge function stuff and call main function
      ! gamma = gammain
      boost = 1d0
      call set_table(1d0)
      call solve_cooling(nH,T2,zsolar,boost,dt,deltaT2,ncell)
      deltaT2out = deltaT2
 
 end subroutine solve_cooling_ramses