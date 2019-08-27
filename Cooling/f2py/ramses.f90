! HACKED TO ALLOW JUMPING OUT WITHOUT HARD STOPPING PROGRAM

module data
  use pm_commons
  use amr_commons
  use hydro_commons
  use rt_parameters
  use rt_cooling_module
  implicit none

! Leaf cell data
  integer,allocatable,dimension(:)::ind_leaf,ind_cell
  integer,allocatable,dimension(:,:)::ind_hydro
  real(kind=8),allocatable,dimension(:,:,:)::hydros
  integer::nleaf,ncache,i,j,ngrid,iskip,size,ivar,ix,iy
  integer,dimension(8)::dx,dy
! PARTICLE DATA
!f2py   real(kind=8),allocatable,dimension(:,:)::xp       ! Positions
!f2py   real(kind=8),allocatable,dimension(:,:)::vp       ! Velocities
!f2py   real(kind=8),allocatable,dimension(:)  ::mp       ! Masses
!f2py   real(kind=8),allocatable,dimension(:)  ::tp       ! Birth epoch
! UNITS DATA (RAW VALUES ONLY)
!f2py   real(kind=8)::units_density
!f2py   real(kind=8)::units_time
!f2py   real(kind=8)::units_length
! AMR DATA
!f2py   integer::ngrid_current
!f2py   real(kind=8),allocatable,dimension(:,:)::xg       ! Position
!f2py   real(kind=8)::t                                   ! Time variable
!f2py   integer::nlevelmax                                ! Max level
!f2py   real(kind=8)::boxlen                              ! Box length
! HYDRO DATA
!f2py   real(kind=8),allocatable,dimension(:,:)::uold     ! Hydro vars
! EXTRAS
!f2py   integer,allocatable,dimension(:)::ind_leaf        ! Leaf IDs
!f2py   real(kind=8),allocatable,dimension(:,:,:)::hydros ! Hydro values
!f2py   integer::nleaf                                    ! Num leaf IDs
! NEEDED FOR RT COOLING
!f2py   integer::ndim                                     ! Number of dimensions
!f2py   integer::nIons                                    ! Number of ion packets

contains
  subroutine startup
    !f2py threadsafe
    !call set_jump_init
    call read_params
    call init_sim
    size = 2**nlevelmax
    allocate(ind_leaf(1:ngridmax))
    allocate(ind_cell(1:ngridmax))
    allocate(ind_hydro(1:ngridmax,1:twotondim))
    allocate(hydros(1:size,1:size,1:nvar))
    nleaf = 0
    ind_leaf=0
    ind_cell=0
    ind_hydro=0
    dx = 0
    dy = 0
    dx(1) = 0 ; dy(1) = 0
    dx(2) = 1 ; dy(2) = 0
    dx(3) = 0 ; dy(3) = 1
    dx(4) = 1 ; dy(4) = 1
  end subroutine startup

  subroutine init_sim
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(kind=8)::n_step
  integer::ilevel,idim,ivar,info,tot_pt
  real(kind=8)::tt1,tt2,muspt,muspt_this_step,wallsec,dumpsec
  real(kind=4)::real_mem,real_mem_tot
  real(kind=8),save::tstart=0
#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
  ! for calculating total run time
  if (tstart.eq.0.0) then
     tstart = MPI_WTIME()
  end if
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(use_grackle==0)then
     if(cooling.and..not.neq_chem) &
        call set_table(dble(aexp))    ! Initialize cooling look up table
  endif
#else
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again
  if(extinction) call init_radiative ! Geometrical corrections in cooling_fine (VV)  

#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration'

  999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

  end subroutine init_sim

  subroutine rt_cooling(T2in, xionin, Npin, Fpin, p_gasin        &
                           & ,nHin, Zsolarin, dtin, a_expin,ncellin &
                           & ,T2out, dNpdtout, dFpdtout)
! Semi-implicitly solve for new temperature, ionization states, 
! photon density/flux, and gas velocity in a number of cells.
! Parameters: 
! T2     <=> T/mu [K] 
! xion   <=> NION ionization fractions 
! Np     <=> NGROUPS photon number densities [cm-3]
! Fp     <=> NGROUPS * ndim photon number fluxes [cm-2 s-1]
! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
! dNpdt   =>  Op split increment in photon densities during dt
! dFpdt   =>  Op split increment in photon flux magnitudes during dt
! nH      =>  Hydrogen number densities [cm-3]
! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
! Zsolar  =>  Cell metallicities [solar fraction]
! dt      =>  Timestep size             [s]
! a_exp   =>  Cosmic expansion
! nCell   =>  Number of cells (length of all the above vectors)
!
! We use a slightly modified method of Anninos et al. (1997).
!-------------------------------------------------------------------------

  implicit none  
  ! WARNING! HARD-CODED SO NDIM=3,NION=3,NGROUPS=5
  !--------------------------------------------------------
  ! Interface variables
  !--------------------------------------------------------
  integer,parameter::nvectorIn=500
  integer,parameter::ndimIn=3
  integer,parameter::nIonsIn=3
  integer,parameter::nGroupsIn=5
  real(kind=8),intent(in),dimension(nvectorIn):: T2in
  real(kind=8),dimension(1:nIonsIn,1:nvectorIn),intent(in):: xionin
  real(kind=8),dimension(1:nGroupsIn,1:nvectorIn),intent(in):: Npin
  real(kind=8),dimension(1:ndimIn,1:nGroupsIn,1:nvectorIn),intent(in):: Fpin
  real(kind=8),dimension(1:ndimIn,1:nvectorIn),intent(in):: p_gasin
  real(kind=8),dimension(1:nvectorIn),intent(in):: nHin, Zsolarin
  !logical,allocatable(:),intent(in):: c_switchin
  real(kind=8),intent(in)::dtin, a_expin
  integer,intent(in):: ncellin
  !--------------------------------------------------------
  real(kind=8),dimension(1:nvectorIn),intent(out)::T2out
  real(kind=8),dimension(1:nGroupsIn,1:nvectorIn),intent(out):: dNpdtout
  real(kind=8),dimension(1:ndimIn,1:nGroupsIn,1:nvectorIn),intent(out):: dFpdtout
  !--------------------------------------------------------
  ! Cooling function variables
  !--------------------------------------------------------
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:nIons, 1:nvector):: xion
  real(dp),dimension(1:nGroups, 1:nvector):: Np, dNpdt
  real(dp),dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
  real(dp),dimension(1:ndim, 1:nvector):: p_gas
  real(dp),dimension(1:nvector):: nH, Zsolar
  logical,dimension(1:nvector):: c_switch
  real(dp)::dt, a_exp
  integer::ncell !--------------------------------------------------------

  ! Set up input
  if(verbose) write(*,*) "Entering cooling"
  ncell = ncellin
  T2(1:ncell) = T2in(1:ncell)
  xion(1:nIons,1:ncell) = xionin(1:nIons,1:ncell)
  Np(1:ngroups,1:ncell) = Npin(1:ngroups,1:ncell)
  Fp(1:ndim,1:ngroups,1:ncell) = Fpin(1:ndim,1:ngroups,1:ncell)
  dNpdt = 0.0d0
  dFpdt = 0d0
  p_gas(1:ndim,1:ncell) = p_gasin(1:ndim,1:ncell)
  nH(1:ncell) = nHin(1:ncell)
  Zsolar(1:ncell) = Zsolarin(1:ncell)
  c_switch = .true. ! Always cool
  dt = dtin
  a_exp = a_expin

  if(verbose) write(*,*) "Solving cooling"
  call rt_solve_cooling(T2, xion, Np, Fp, p_gas, dNpdt, dFpdt        &
                           ,nH, c_switch, Zsolar, dt, a_exp, nCell)
  if(verbose) write(*,*) "Cooling solved"

  ! Copy outputs
  ! allocate(T2out(ncell))
  ! allocate(dNpdtout(1:nGroups,1:ncell))
  ! allocate(dFpdtout(1:ndim,1:ngroups,1:ncell))
  T2out(1:ncell) = T2(1:ncell)
  dNpdtout(1:nGroups,1:ncell) = dNpdt(1:nGroups,1:ncell)
  dFpdtout(1:ndim,1:ngroups,1:ncell) = dFpdt(1:ndim,1:ngroups,1:ncell)

  end subroutine rt_cooling

  ! subroutine step
  !   implicit none
  !   integer::nx_loc
  !   real(dp),dimension(1:3)::skip_loc
  !   !f2py threadsafe
  !   ! Call step
  !   call set_jump_step
  !   ! Make extra stuff for the python modules
  !   nleaf=0
  !   ngrid=active(nlevelmax)%ngrid
  !   size = 2**nlevelmax
  !   nx_loc=(icoarse_max-icoarse_min+1)
  !   skip_loc=(/0.0d0,0.0d0,0.0d0/)
  !   if(ndim>0)skip_loc(1)=dble(icoarse_min)
  !   if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  !   if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  !   do i=1,ngrid
  !      ind_cell(i)=active(nlevelmax)%igrid(i-1)
  !      !if(son(ind_cell(i))==0)then
  !         nleaf=nleaf+1
  !         ind_leaf(i)=ind_cell(i)
  !         do j=1,twotondim
  !            iskip=ncoarse+(j-1)*ngridmax
  !            ind_hydro(i,j)=ind_leaf(i)+iskip
  !            ix = int((xg(ind_cell(i),1)-0.5-skip_loc(1))*(size/2))+size/4
  !            iy = int((xg(ind_cell(i),2)-0.5-skip_loc(2))*(size/2))+size/4
  !            ix = ix*2+dx(j)
  !            iy = iy*2+dy(j)
  !            !write(*,*)ix,iy
  !            !write(*,*),xg(ind_cell(i),1),xg(ind_cell(i),2)
  !            if ((ix > 0) .and. (ix <= size) .and. &
  !                 (iy > 0) .and. (iy <= size)) then
  !               do ivar=1,nvar
  !                  hydros(ix,iy,ivar) = uold(ind_leaf(i)+iskip,ivar)
  !               end do
  !            endif
  !         end do
  !      !end if
  !   end do
  ! end subroutine step
    
end module data
