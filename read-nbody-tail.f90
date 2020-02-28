!**************************************************************
! read_nbody
! Author: S.Mohammad Hoseini Rad
! smhr313@gmail.com
! Nov 2012, IASBS, Zanjan
! Last modification: February 28, 2019 / 9 Esfand 1398
!**************************************************************
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the
!     Free Software Foundation, Inc.,
!     59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


! ******************************************
program read_nbody
   use kdtree2_module
   implicit none
! ******************************************
   real*8::start ! used for cpu_time
   integer*4::i,j,k,ll ! loop counter variables
   integer*4::IO,code, err, debug
   integer::tscreen ! frequency of output on screen
   integer::tout ! frequency of output on harddisk
   integer,dimension(13)::buff ! used for measuring size of input file
   integer*4::status ! used for measuring size of input file
   integer::loop_index
   integer::INIT_NTOT
   integer::INIT_NTAIL
! ******************************************
! NBODY6 outputs (here inputs). see output.f in Ncode directory of NBODY6
! before removing stars with NAME < 1 and outside tidal rdius.
   integer*4::iiNTOT
   integer,dimension(:),allocatable::iiNAME,iiKSTAR
   real*8,dimension(:),allocatable::iiBODYS,iiRADIUS,iiZLMSTY
   real*8,dimension(:,:),allocatable::iiXS,iiVS
! ******************************************
   integer*4::iNTOT
   integer,dimension(:),allocatable::iNAME,iKSTAR
   real*8,dimension(:),allocatable::iBODYS,iRADIUS,iZLMSTY
   real*8,dimension(:,:),allocatable::iXS,iVS
! ******************************************
! NBODY6 tidal tail outputs (here inputs). see output.f in Ncode directory of NBODY6
! before removing stars with NAME < 1
   integer*4::iNTAIL
   integer,dimension(:),allocatable::iNAMEt
   real*8,dimension(:),allocatable::iBODYSt
   real*8,dimension(:,:),allocatable::iXSt,iVSt
! ******************************************
! NBODY6 outputs (here inputs), after removing stars with NAME < 1 and outside tidal rdius
   integer*4::NTOT,NK,N,MODEL,NRUN
   integer,dimension(:),allocatable::NAME,KSTAR
   real*8,dimension(:),allocatable::AS,BODYS,RADIUS,ZLMSTY
   real*8,dimension(:,:),allocatable::XS,VS
   real*4,dimension(:,:),allocatable::XSS,VSS
   real*4,dimension(:),allocatable::BODYSS,ASS  ! used for converting single precision to double precision ( if code = 2)
! ******************************************
! NBODY6 tidal tail outputs (here inputs), after removing stars with NAME < 1
   integer*4::NTAIL
   integer,dimension(:),allocatable::NAMEt
   real*8,dimension(:),allocatable::ASt,BODYSt
   real*8,dimension(:,:),allocatable::XSt,VSt
   real*4,dimension(:,:),allocatable::XSSt,VSSt
   real*4,dimension(:),allocatable::BODYSSt,ASSt  ! used for converting single precision to double precision ( if code = 2)
! ******************************************
   real*8,dimension(3)::RG,VG
   real*8::rrt,OMEGA,mtot_temp,rt_temp
   integer::ntot_temp
! ******************************************
   real*8::mtot, imtot, T6, mtot0
   integer,dimension(:,:),allocatable::list_neighbor_idx ! neighbor index array
   real,dimension(:,:),allocatable::list_neighbor_dis ! neighbor distance array
   character(len=100)::output_file,input_file
   character(len=100)::output_file_tail
   integer::n_neighbor !number of nearest neighbors in Casterano & Hut method
   real*8,dimension(3)::Xd ! density center coordinates
   real*8,dimension(5)::LR ! Lagrangian radii array
   real*8,dimension(:),allocatable::r,ir ! position vector of stars after and befor tidal radius detemination
! ******************************************
   type(kdtree2), pointer::tree
   real, allocatable:: mydata(:,:)
   type(kdtree2_result), allocatable :: results(:)
! ******************************************
   real*8,dimension(:),allocatable::irho, rho, rho2
   real*8::sixth_neighbor_distance, m_tot_5th, irho_tot, rho2_tot_in_rh
   real*8,dimension(3)::irho_times_dis_tot
! ******************************************
   real*8::mtemp, rtemp, rstep, rc, rt
   integer::check, diss_check, major_output
! ******************************************
   real::Ndiss, Mdiss
! ******************************************
   character(len=40)::dir_command
! ******************************************
   real*8::rx, ry, rz, vx, vy, vz, Ebin, a, ecc, ecc2, Lx, Ly, Lz, fbin0, fbint
   real*8,dimension(:,:),allocatable::del_r, del_v2, nei_BODYS, L2
   real*8,dimension(:,:),allocatable::idel_r, idel_v2, inei_bodys, iL2
   integer*8,dimension(:,:),allocatable::nei_NAME, inei_NAME
   integer::nid, find_binary, Nbin
   real*8::binary_energy_criterion
! *****************************************
   real*8::Mstar, Tstar, Rstar, Vstar
   real,parameter::G = 0.0043009211 ! In pc*km^2/s^2*M_sun
   real,parameter::AU = 206264.806 ! pc to AU
! ******************************************
   integer:: res_INIT_NTOT
   real*8:: res_mtot0
! ******************************************

   call cpu_time(start)

! ******************************************
!** Initialization of internal variables ***
! Do not change these variables
   IO = 0
   err = 0
   status = 0
   n_neighbor = 6 ! number of nearest neighbors in Casterano & Hut (1985) and Aarseth (2001) method
   loop_index = 0
   INIT_NTOT = 0; iNTOT = 0; imtot = 0
   rstep = 0.01 ! bin length in pc for calculating lagrangian radii
   res_INIT_NTOT = 0
   res_mtot0 = 0.

! ************ Code options ****************
! ******************************************
   code = 4         ! 1:NBODY6 custom, 2: NBODY6, 3: NBODY6 new custom ver, 4: NBODY6 new size_scale ver, 5: NBODY6 size_scale ver with separate OUT33
   tscreen = 1      ! time interval of output on screen & harddisk. (NBODY6 custom: Myr, NBODY6: N-body unit)
   tout = 1    ! time interval of output on screen & harddisk (suppressed)
   ! it also prevent neighbor arrays and kdtree2 pointers to be allocated
   ! so if you have many time snapshots, by increasing 'tout'
   ! you can pass the memory overflow problem. (NBODY6 custom: Myr, NBODY6: N-body unit)
   debug = 1        ! 1: debug mode
   diss_check = 1   ! 1: check dissolution of cluster; 0: no check
   Ndiss = 0        ! by reaching to this fraction of initial number of stars, the cluster is considered as a dissolved cluster;
   ! 0: Suppress this option
   Mdiss = 0.001    ! By reaching to this fraction of initial mass, the cluster is considered as a dissolved cluster;
   ! 0: Suppress this option
   major_output = 1
! Please note just use one of Ndiss or Mdiss options
   find_binary = 1  ! 0: Skip to find binary; 1: Find binaries
   binary_energy_criterion = -0.001
   mtot0 = 0.

   ! ********* Input & Output files. *********
   ! *****************************************
   call get_model_name ( input_file, res_INIT_NTOT, res_mtot0 )
   open(1,file=input_file, form='unformatted')
   open(2,file='overview.txt') ! this file includes some cluster characteristics
   open(3,file='density_center.txt') ! this file includes density center coordinate of cluster during evolution
   open(7,file='radii.txt') ! this file includes Lagrangian radii
   open(10,file='binary.txt')
   open(33,file='OUT33', form='unformatted')
! ******************************************
   dir_command = 'mkdir snapshot'
   call system ( dir_command )
   dir_command = 'mkdir tail'
   write(*,*)dir_command
   call system ( dir_command )

   call STAT (input_file, buff,status) ! calculating size of input file
   if (status == 0) then
      write(*,*)
      write(*,*)"Size of file is", buff(8), "Byte."
   endif
! ******************************************
   write(2,*) "    T_NB    T_Myr   N        M      M_ratio    Rh       Rt       Rc     Rc_O_Rh     RC    Rh_O_Rt   fbin0    fbint"
   write(7,*)" T(NBODY) Time(Myr) LR_0.01   LR_0.10   LR_0.30  LR_0.50   LR_0.70     Rt         Rc       Rc/Rh      RC      Rh/Rt"

! ************ Read in loop ****************
! *******************************************
! do l=1,10

! By activating this loop, you can define exact time in terms of output interval (DELTAT in NBODY6 or dtscreen in mcluster [N-body units (Myr in Nbody6 custom)]).


   DO WHILE ( .TRUE. ) ! repeating the loop until the end of input file

! *******************************************
! define_arrays_size ( code, input_file, loop_index, NTOT, MODEL, NRUN, NK, N, INIT_NTOT ) ! Read number of particles in each loop.
      if (code == 1 .or. code == 3 .or. code == 4) then
         read(1, iostat=IO)iiNTOT,NK,N
         if (debug == 1) write (*,'(3i10)') iiNTOT,NK,N
         if (loop_index == 0) INIT_NTOT = iiNTOT
      elseif (code == 2 .or. code == 5) then
         read(1, iostat=IO)iiNTOT, MODEL, NRUN, NK, N
         if (loop_index == 0) INIT_NTOT = iiNTOT
         if (debug == 1) write (*,'(5i10)') iiNTOT, MODEL, NRUN, NK, N
      endif
      if (code == 5) then
         read(33, iostat=IO)iNTAIL
         if (loop_index == 0) INIT_NTAIL = iNTAIL
         if (debug == 1) write (*,'(2i10)') iNTAIL
      endif
      if (IO /= 0) then
         call termination (IO, err)
      endif
      if (res_INIT_NTOT /= 0) then
         INIT_NTOT = res_INIT_NTOT
      endif

! *******************************************
! this condition is for appropriate termination of the program if cluster dissolves.
      if ( iiNTOT > 0 ) then
         if ( diss_check == 1 ) then
            if ( Ndiss > 0  ) then
               if ( iiNTOT <= Ndiss * INIT_NTOT ) then
                  err = 2
                  call termination ( IO, err )
               endif
            endif
            if ( loop_index == 1) imtot = mtot
            if ( imtot > 0 ) then
               if ( Ndiss == 0 .and. Mdiss > 0 ) then
                  if ( mtot <= Mdiss * imtot ) then
                     err = 2
                     call termination ( IO, err )
                  endif
               elseif ( Ndiss == 0 .and. Mdiss == 0 ) then
                  write (*,*) "Warning! Both Ndiss and Mdiss are zero!!!"
                  write (*,*) "Please input 0 and press Enter to continue."
                  read (*,*) diss_check
                  if ( diss_check /= 0 ) then
                     write (*,*) "Progrom is stoped!"
                     STOP
                  endif
               endif
            endif
         endif
      else
         write (*,*) "iiNTOT <= 0. Progrom is stoped!"
         STOP
      endif

!***********************************
! Read input arrays according to the output.f in Ncode directory of NBODY6. Note that there are some differences between NBODY6 custom and NBODY6 input variables.
!************************************
      allocate (AS(NK)); allocate (iiBODYS(iiNTOT)); allocate(iiXS(3,iiNTOT)); allocate(iiVS(3,iiNTOT))
      allocate(iiRADIUS(iiNTOT)); allocate(iiNAME(iiNTOT)); allocate(iiKSTAR(iiNTOT)); allocate(iiZLMSTY(iiNTOT))

      allocate (ASS(NK)); allocate (BODYSS(iiNTOT)); allocate(XSS(3,iiNTOT)); allocate(VSS(3,iiNTOT))

      if (code == 5) then
         allocate (ASt(13)); allocate (iBODYSt(iNTAIL)); allocate(iXSt(3,iNTAIL)); allocate(iVSt(3,iNTAIL))
         allocate (ASSt(13)); allocate (BODYSSt(iNTAIL)); allocate(XSSt(3,iNTAIL)); allocate(VSSt(3,iNTAIL))
         allocate(iNAMEt(iNTAIL))
      endif

      if (debug == 1) print*,"Allocating readin arrays."

      iiRADIUS = 0 ;iiZLMSTY = 0; iiKSTAR = 0
      T6 = 0; Mstar =0; Rstar = 0; Vstar = 0; Tstar = 0

      if (code == 1 .or. code == 3 .or. code == 4) then
         read(1, iostat=IO)(AS(K),K=1,NK),(iiBODYS(J),J=1,iiNTOT),((iiXS(K,J),K=1,3),J=1,iiNTOT),((iiVS(K,J),K=1,3),J=1,iiNTOT),&
            (iiRADIUS(J),J=1,iiNTOT),(iiNAME(J),J=1,iiNTOT),&
            (iiKSTAR(J),J=1,iiNTOT),(iiZLMSTY(J),J=1,iiNTOT)

         if (debug == 1) write (*,*)(k,AS(K),K=1,NK)
         if (code == 1) then
            rt = AS(25) ! tidal Radius.
            T6 = AS(10) ! Astrophysical time (Myr).
            Mstar = AS(4)
            Rstar = AS(3)
            Vstar = AS(12)
         elseif (code == 3 .or. code == 4) then ! new nbody6 version
            if (AS(1) < 0.0001) rt = AS(5) ! tidal Radius.
            Tstar = AS(11)
            T6 = AS(1) * Tstar ! Astrophysical time (Myr), Nbody time * TSCALE
            Mstar = AS(4)
            Rstar = AS(3)
            Vstar = AS(12)
!             RG(1) = AS(21)
!             RG(2) = AS(22)
!             RG(3) = AS(23)
!             VG(1) = AS(24)
!             VG(2) = AS(25)
!             VG(3) = AS(26)
            OMEGA = AS(6)/2.
         endif
         if (IO /= 0) then
            err = 1
            call termination (IO, err)
         endif
      elseif (code == 2 .or. code == 5) then
         read(1, iostat=IO)(ASS(K),K=1,NK),(BODYSS(J),J=1,iiNTOT),((XSS(K,J),K=1,3),J=1,iiNTOT),((VSS(K,J),K=1,3),J=1,iiNTOT),&
                                           (iiNAME(J),J=1,iiNTOT)
         if (debug == 1 ) write (*,'(i3,f20.9)')(k,ASS(K),K=1,NK)
         do J=1, NK ! Convert to double precision
            AS(J) = ASS(J)
         enddo
         do J=1, iiNTOT ! convert to double precision
            iiBODYS(J) = BODYSS(J)
            do K=1, 3
               iiXS(K,J) = XSS(K,J)
               iiVS(K,J) = VSS(K,J)
            enddo
         enddo
         rt = AS(5) ! tidal Radius.
         Tstar = AS(11)
         T6 = AS(1) * Tstar ! astrophysical time (Myr), Nbody time * TSCALE
         Mstar = AS(4)
         Rstar = AS(3)
         Vstar = AS(12)
         if (IO /= 0) then
            err = 1
            call termination (IO, err)
         endif
      endif
      if (code == 5) then
         ! Read tail data
        read(33, iostat=IO)(ASSt(K),K=1,13),(BODYSSt(J),J=1,iNTAIL),((XSSt(K,J),K=1,3),J=1,iNTAIL),((VSSt(K,J),K=1,3),J=1,iNTAIL),&
                                            (iNAMEt(J),J=1,iNTAIL)
        if (debug == 1 ) write (*,'(i3,f20.9)')(k,ASSt(K),K=1,13)
        do J=1, 13 ! Convert to double precision
           ASt(J) = ASSt(J)
        enddo
        do J=1, iNTAIL ! convert to double precision
           iBODYSt(J) = BODYSSt(J)
           do K=1, 3
              iXSt(K,J) = XSSt(K,J)
              iVSt(K,J) = VSSt(K,J)
           enddo
        enddo
         rt = AS(5) ! tidal Radius.
         Tstar = AS(11)
         T6 = AS(1) * Tstar ! astrophysical time (Myr), Nbody time * TSCALE
         Mstar = AS(4)
         Rstar = AS(3)
         Vstar = AS(12)
         if (IO /= 0) then
            err = 1
            call termination (IO, err)
         endif
      endif
!************************************
      if (debug == 1) print*,"Checking dtout option."
      if (code == 1 .or. code == 3) then
         if (mod (int(nint(AS(10))),tout) /= 0)then
            write(*,*) "No calculation for time = ", AS(10)
            deallocate (AS, iiBODYS, iiXS, iiVS, iiRADIUS, iiNAME, iiKSTAR, iiZLMSTY)
            deallocate (ASS, BODYSS, XSS, VSS)
            cycle
         endif
      elseif (code == 2 .or. code == 5 .or. code == 4) then
         if (mod (int(nint(AS(1))),tout) /= 0) then
            write(*,*) "No calculation for time = ", AS(1)
            if (code == 5) then
               write(*,*) "No calculation for time = ", AS(1)
               deallocate (ASt, iBODYSt, iXSt, iVSt, iNAMEt)
               deallocate (ASSt, BODYSSt, XSSt, VSSt)
            endif
            deallocate (AS, iiBODYS, iiXS, iiVS, iiRADIUS, iiNAME, iiKSTAR, iiZLMSTY)
            deallocate (ASS, BODYSS, XSS, VSS)
            cycle
         endif
      endif
!***********************************
! remove stars with index < 1
!***********************************
      i = 0; k = 0; iNTOT = 0; NTAIL = 0
      if (debug == 1 ) write (*,*) "Allocating ir Array in test_name."
      if (debug == 1 ) write (*,*) "INIT_NTOT = ", INIT_NTOT

      do i = 1, iiNTOT
         if ( iiNAME(i) >= 1 .AND. iiNAME(i) <= INIT_NTOT .AND. iiBODYS(i) > 0 ) then
            k = k + 1
         endif
      enddo
      iNTOT = k
      print*,"iNTOT =", iNTOT
      k = 0
      if (code == 5) then
         do i = 1, iNTAIL
            if ( iNAMEt(i) >= 1 .AND. iNAMEt(i) <= INIT_NTOT .AND. iBODYSt(i) > 0 ) then
               k = k + 1
            endif
         enddo
         NTAIL = k
         if (debug == 1 ) write (*,*) "NTAIL =", NTAIL
      endif
      
      allocate (iBODYS(iNTOT)); allocate(iXS(3,iNTOT)); allocate(iVS(3,iNTOT))
      allocate(iRADIUS(iNTOT)); allocate(iNAME(iNTOT)); allocate(iKSTAR(iNTOT)); allocate(iZLMSTY(iNTOT))

      if (debug == 1 ) write (*,*)"Allocating test_name."
      k = 0
      do i = 1, iiNTOT
         if ( iiNAME(i) >= 1 .AND. iiNAME(i) <= INIT_NTOT .AND. iiBODYS(i) > 0 ) then
            k = k + 1
            iBODYS(k) = iiBODYS(i)
            iXS(1,k) = iiXS(1,i); iXS(2,k) = iiXS(2,i); iXS(3,k) = iiXS(3,i)
            iVS(1,k) = iiVS(1,i); iVS(2,k) = iiVS(2,i); iVS(3,k) = iiVS(3,i)
            iRADIUS(k) = iiRADIUS(i)
            iNAME(k) = iiNAME(i)
            iKSTAR(k) = iiKSTAR(i)
            iZLMSTY(k) = iiZLMSTY(i)
         endif
      enddo
      iiNTOT = 0
      deallocate(iiBODYS, iiXS, iiVS, iiRADIUS , iiNAME, iiKSTAR, iiZLMSTY)
      if (code == 5) then
         allocate (BODYSt(NTAIL)); allocate(XSt(3,NTAIL)); allocate(VSt(3,NTAIL)); allocate(NAMEt(NTAIL))
         k = 0
         do i = 1, iNTAIL
            if ( iNAMEt(i) >= 1 .AND. iNAMEt(i) <= INIT_NTOT .AND. iBODYSt(i) > 0 ) then
               k = k + 1
               BODYSt(k) = iBODYSt(i)
               XSt(1,k) = iXSt(1,i); XSt(2,k) = iXSt(2,i); XSt(3,k) = iXSt(3,i)
               VSt(1,k) = iVSt(1,i); VSt(2,k) = iVSt(2,i); VSt(3,k) = iVSt(3,i)
               NAMEt(k) = iNAMEt(i)
            endif
         enddo
         iNTAIL = 0
         deallocate(iBODYSt, iXSt, iVSt, iNAMEt)
     endif
!************************************
! finding all neighbors and their distances for all stars
! and saving them in list_neighbor_idx and list_neighbor_dis arrays correspondingly
!************************************
      if (debug == 1 ) write (*,*) "Finding neighbor with iNTOT = ", iNTOT
      allocate(mydata(3,iNTOT))
      mydata = iXS ! putting x,y and z of all stars in mydata array to change them to single precision.
      tree => kdtree2_create(mydata,rearrange=.true.,sort=.true.)
      allocate (list_neighbor_idx (iNTOT,n_neighbor + 1))
      allocate (list_neighbor_dis (iNTOT,n_neighbor + 1))

      do i=1,iNTOT ! find star neighbors and their distances to the ith star in XSt array
         allocate (results(n_neighbor + 1)) ! results(1) is associated with the ith star
         !results(2..n_neighbor+1) are associated with n_neighbor nearest neighbors
         ! so for 6 nearest neighbors, size of results array must be 7.

         call kdtree2_n_nearest_around_point(tree,idxin=i,nn=n_neighbor + 1,correltime=0,results=results)

         do ll=1,n_neighbor + 1
            list_neighbor_idx (i,ll) = results (ll)%idx
            list_neighbor_dis (i,ll) = sqrt(results (ll)%dis)
         enddo
         deallocate (results)
      enddo
!************************************
! finding density center of cluster according to the Casterano & Hut (1985) and Aarseth 2001 method
! for more details see center_find subroutine
!************************************
      if (debug == 1 ) write (*,*) "Finding density center."

      allocate(irho(iNTOT)); allocate (idel_r(iNTOT, n_neighbor + 1))
      allocate(idel_v2(iNTOT, n_neighbor + 1)); allocate(inei_bodys(iNTOT, n_neighbor + 1))
      allocate(inei_NAME(iNTOT, n_neighbor + 1)); allocate(iL2(iNTOT, n_neighbor + 1))

      irho = 0; irho_tot = 0; irho_times_dis_tot = 0; Xd = 0
      do i=1,iNTOT
         m_tot_5th=0; sixth_neighbor_distance=0
         sixth_neighbor_distance=sqrt(list_neighbor_dis (i,n_neighbor + 1))
         do ll=2, n_neighbor
            m_tot_5th=m_tot_5th+iBODYS(list_neighbor_idx (i,ll)) ! total mass of fifth nearest neibors
         enddo
         irho(i) = m_tot_5th/(sixth_neighbor_distance ** 3.)
         irho_tot = irho_tot + irho(i)
         irho_times_dis_tot(1) = irho_times_dis_tot(1) + irho(i) * iXS(1,i)
         irho_times_dis_tot(2) = irho_times_dis_tot(2) + irho(i) * iXS(2,i)
         irho_times_dis_tot(3) = irho_times_dis_tot(3) + irho(i) * iXS(3,i)
      enddo
      Xd(1) = ( irho_times_dis_tot(1)/irho_tot )
      Xd(2) = ( irho_times_dis_tot(2)/irho_tot )
      Xd(3) = ( irho_times_dis_tot(3)/irho_tot )
      if (debug == 1 ) then
         write (*,*) "Xd:", (Xd(i),i=1,3)
         write (*,*) "Done."
      endif
! **********************************
      idel_r = 0; idel_v2 =0; inei_bodys = 0; inei_NAME = 0; iL2 = 0
      do i =1, iNTOT
         ll = 0
         do ll=2, n_neighbor+1
            nid = list_neighbor_idx(i,ll) ! Index of nearest neighbor.
            inei_BODYS(i,ll) = iBODYS(nid)
            inei_NAME(i,ll) = iNAME(nid)
            rx = iXS(1,i) - iXS (1,nid)
            ry = iXS(2,i) - iXS (2,nid)
            rz = iXS(3,i) - iXS (3,nid)
            idel_r(i,ll) = sqrt ( rx*rx + ry*ry + rz*rz )
            vx = iVS(1,i) - iVS (1,nid)
            vy = iVS(2,i) - iVS (2,nid)
            vz = iVS(3,i) - iVS (3,nid)
            idel_v2(i,ll) = ( vx*vx + vy*vy + vz*vz )
            Lx = ry*vz - rz*vy
            Ly = rz*vx - rx*vz
            Lz = rx*vy - ry*vx
            iL2(i,ll) = lx*lx + ly*ly + lz*lz
         enddo
      enddo
!***********************************
! correcting the origin of coordinate system by subtracting all star coordinates from measured density center
!***********************************
      do j = 1, iNTOT
         i = 0
         do i = 1, 3
            iXS(i,j) = iXS(i,j) - Xd(i)
         enddo
      enddo
!***********************************
! test rt
!***********************************
      i = 0; k = 0; NTOT = 0
      allocate(ir(iNTOT))
      if (debug == 1 ) write (*,*) "Testing rt."
      if ( AS(1) > 0.0001) then ! This condition is for tidal radius to be correct at t=0 and includes all stars.
         do i = 1, iNTOT
            ir(i) = sqrt( iXS(1,i)*iXS(1,i) + iXS(2,i)*iXS(2,i) + iXS(3,i)*iXS(3,i) )
            if ( ir(i) <= 1e20*rt ) then
               k = k + 1
            endif
         enddo
         NTOT = k
      else
         NTOT = iNTOT
      endif

      if (debug == 1 ) write (*,*)"	Allocating test_rt arrays."
      allocate (BODYS(NTOT)); allocate(XS(3,NTOT)); allocate(VS(3,NTOT))
      allocate(RADIUS(NTOT)); allocate(NAME(NTOT)); allocate(KSTAR(NTOT)); allocate(ZLMSTY(NTOT))
      allocate(r(NTOT)); allocate (rho(NTOT))
      allocate (del_r(NTOT, n_neighbor + 1)); allocate(del_v2(NTOT, n_neighbor +1 ))
      allocate(nei_bodys(NTOT, n_neighbor + 1)); allocate(nei_NAME(NTOT, n_neighbor + 1))
      allocate(L2(NTOT, n_neighbor + 1))
      if (debug == 1 ) write (*,*) "	Done."

      del_r = 0; del_v2 = 0; nei_bodys = 0; nei_NAME = 0
      if ( AS(1) > 0.0001) then
         k = 0
         do i = 1, iNTOT
            ir(i) = sqrt( iXS(1,i)*iXS(1,i) + iXS(2,i)*iXS(2,i) + iXS(3,i)*iXS(3,i) )
            if ( ir(i) <= 1e20*rt ) then
               k = k + 1
               r(k) = ir(i)
               BODYS(k) = iBODYS(i)
               XS(1,k) = iXS(1,i); XS(2,k) = iXS(2,i); XS(3,k) = iXS(3,i)
               VS(1,k) = iVS(1,i); VS(2,k) = iVS(2,i); VS(3,k) = iVS(3,i)
               RADIUS(k) = iRADIUS(i)
               NAME(k) = iNAME(i)
               KSTAR(k) = iKSTAR(i)
               ZLMSTY(k) = iZLMSTY(i)
               rho(k) = irho(i)
               ll = 0
               do ll = 2, n_neighbor + 1
                  del_r(k, ll) = idel_r(i, ll)
                  del_v2(k, ll) = idel_v2(i, ll)
                  L2(k, ll) = iL2(i, ll)
                  nei_bodys(k, ll) = inei_bodys(i , ll)
                  nei_NAME(k, ll) = inei_NAME(i, ll)
               enddo
            endif

         enddo

      else
         BODYS = iBODYS; XS = iXS; VS = iVS; RADIUS = iRADIUS
         NAME = iNAME; KSTAR = iKSTAR; ZLMSTY = iZLMSTY; rho = irho
         del_r = idel_r; del_v2 = idel_v2; nei_bodys = inei_bodys
         nei_NAME = inei_NAME; L2 = iL2
         do i=1, NTOT
            r(i) = sqrt( iXS(1,i)*iXS(1,i) + iXS(2,i)*iXS(2,i) + iXS(3,i)*iXS(3,i) )
         enddo
      endif

      deallocate  ( iBODYS, iXS, iVS, iRADIUS , iNAME, iKSTAR, iZLMSTY, ir )
      deallocate ( idel_r, idel_v2, inei_bodys, inei_NAME, iL2 )
      if (debug == 1 ) write (*,*) "Done."
!***********************************
! Converting from Nbody units to Astrophysical units
!************************************
      if (debug == 1 ) write (*,*) "Converting from Nbody units to Astrophysical units."
      mtot = 0
      do i = 1, NTOT ! Converting from Nbody units to Astrophysical units:
         BODYS(i) = BODYS(i) * Mstar ! Converting masses to Solar Mass.
         mtot = mtot + BODYS(i)
         r(i) = r(i)  *Rstar ! Converting distances to pc.
         rho(i) = rho(i) * Mstar / Rstar / Rstar / Rstar
         do k = 1, 3 ! Converting coordinates to pc & velocities to km/s.
            XS(k,i)=XS(k,i)*Rstar; VS(k,i)=VS(K,i)*Vstar
         enddo
      enddo
      if (AS(1) < 0.0001) rt = rt * Rstar
      OMEGA = OMEGA*Vstar/Rstar ! Convert to Astrophysical unit (kms^-1/pc)
      ! compute new tidal radius
      if (code == 4 .AND. AS(1) > 0.0001) then
      if (debug == 1 ) write (*,*) "Computing new tidal radius."
         mtot_temp = 0.
         rt_temp = 0.
         do while (abs(rt-rt_temp)>0.01)
            if (debug == 1 ) print*,rt,rt_temp,mtot_temp,"###########"
            rt_temp = rt
            mtot_temp = 0.
            ntot_temp = 0
            do i = 1, NTOT
               if (r(i) <= rt_temp) then
                  mtot_temp = mtot_temp + BODYS(i)
                  ntot_temp = ntot_temp + 1
               endif
            enddo
            rt = (G*mtot_temp/(2.*OMEGA*OMEGA))**0.3333 ! In Astrophysical unit (pc)
         enddo
         NTOT = ntot_temp
         mtot = mtot_temp ! This will rewrite mtot for the second and its subsequent snapshots
         if (debug == 1 ) write(*,*)"omega, rt = ",OMEGA,rt
      endif
      if (code == 5) then
         do i = 1, NTAIL ! Converting from Nbody units to Astrophysical units (tidal tail):
            BODYSt(i) = BODYSt(i) * Mstar
            do k = 1, 3 ! Converting coordinates to pc & velocities to km/s.
               XSt(k,i)=XSt(k,i)*Rstar; VSt(k,i)=VSt(K,i)*Vstar
            enddo
         enddo
      endif

!***********************************
! compute Lagrangian radii and changing units
!************************************
      rtemp = 0; mtemp = 0; check = 1
      if (debug == 1 ) write (*,*) "Calculating LRs."
      if ( AS(1) < 0.0001) mtot0 = mtot
      if ( res_mtot0 /= 0. ) mtot0 = res_mtot0

      do while ( check /= 6 )
         mtemp = 0
         do i = 1, NTOT
            if ( r(i) <= rtemp ) then
               mtemp = mtemp + BODYS(i)
               if  ( mtemp >= mtot * 0.01 .AND. check == 1)  then
                  LR(1) = rtemp;! print*,"LR1 = ",LR(1)
                  check = 2; exit
               elseif ( mtemp >= ( mtot * 0.10 ) .AND. check == 2 ) then
                  LR(2) = rtemp;! print*,"LR5 = ",LR(2)
                  check = 3; exit
               elseif ( mtemp >= ( mtot * 0.30 ) .AND. check == 3 ) then
                  LR(3) = rtemp;! print*,"LR50 = ",LR(3)
                  check = 4; exit
               elseif ( mtemp >= ( mtot * 0.50 ) .AND. check == 4 ) then
                  LR(4) = rtemp;! print*,"LR75 = ",LR(4)
                  check = 5; exit
               elseif ( mtemp >= mtot * 0.70 .AND. check == 5 ) then
                  LR(5) = rtemp;! print*,"LR95 = ",LR(5)
                  check = 6; exit
               endif
            endif
         enddo
         rtemp = rtemp + rstep
      enddo
      if (debug == 1 ) write (*,*)"Calculating LRs. Done."
!***********************************
! compute core radius (Trenti et al. 2007).
!************************************
      allocate (rho2(NTOT))
      rc = 0; rho2 = 0; rho2_tot_in_rh = 0
      do i=1, NTOT
         if ( r(i) <= LR(3) ) then
            rho2(i) = rho(i) * rho(i)
            rc = rc + r(i) * r(i) * rho2(i)
            rho2_tot_in_rh = rho2_tot_in_rh + rho2(i)
         endif
      enddo
      rc = sqrt( rc / rho2_tot_in_rh )
      if (debug == 1 ) write (*,*)"Calculating Core Radius. Done."
!***********************************
! finding binaries.
!************************************
      if ( find_binary == 1 ) then
         if (debug == 1 ) write (*,*)"Finding binaries."
         Ebin = 0; a = 0; ecc = 0; ecc2 = 0; Nbin = 0; fbin0 = 0; fbint = 0; j = 0
         do i = 1, NTOT
            ll=0
            do ll=2, n_neighbor+1
               nei_BODYS(i, ll) = nei_BODYS(i, ll) * Mstar
               del_r(i, ll) = del_r(i, ll) * Rstar
               del_v2(i, ll) = del_v2(i, ll) * Vstar * Vstar
               L2(i, ll) = L2(i, ll) * Rstar * Vstar * Rstar * Vstar
               Ebin = - ( G * BODYS(i) * nei_BODYS(i, ll) / del_r(i, ll) ) &
               & + 0.5 *  BODYS(i) * nei_BODYS(i, ll) / ( BODYS(i) + nei_BODYS(i ,ll) )&
               & * del_v2(i, ll)
               a = - 0.5 * G * BODYS(i) * nei_BODYS(i, ll) / Ebin
               ecc2 = 1. - L2(i, ll) / ( G * a * ( BODYS(i) + nei_BODYS(i, ll) ) )
               if ( ecc2 >= 0 ) then
                  ecc = sqrt(ecc2)
               else
                  ecc = 0
               endif
               if ( Ebin < binary_energy_criterion ) then
                  if ( NAME(i) == nei_NAME(i, ll) ) cycle
                  Nbin = Nbin + 1
                  write (10,'(2f10.3, 2i10, 6f12.3)') AS(1), T6, name(i), nei_NAME(i ,ll), BODYS(i),&
                  & nei_BODYS(i ,ll), ecc, a*AU, del_r(i, ll)*AU, Ebin
               endif
            enddo
         enddo

         Nbin = Nbin / 2  ! because each binary counted twice

         if (debug == 1 ) then
            write(*,*)Nbin, "binaries were found."
            write(*,*)"Done."
         endif
      endif !  end of main binary finding condition

      fbin0 = float( Nbin ) / N
      fbint =  float( Nbin ) / NTOT
!************************************
! writing out outputs
!************************************
      if ( major_output == 1 ) then
         Write(output_file , '( f7.1 )' ) T6
         if ( nint(T6) < 10 ) then
            output_file = '0000' // trim(output_file) ! // ".txt"
         else if ( nint(T6) < 100 ) then
            output_file = '000' // trim(output_file)
         else if ( nint(T6) < 1000 ) then
            output_file = '00' // trim(output_file)
         else if ( nint(T6) < 10000 ) then
            output_file = '0' // trim(output_file)
         else
            output_file = trim(output_file)
         endif
         if (code == 5 .or. code == 4) output_file_tail = 'tail/' // trim(output_file)
         output_file = 'snapshot/' // trim(output_file)
         output_file = sweep_blanks(output_file)
         if (code == 5 .or. code == 4) output_file_tail = sweep_blanks(output_file_tail)
         open(4,file=output_file)
         if (code == 5 .or. code == 4) open(8,file=output_file_tail)
         do i=1,NTOT
            if (code == 1 .or. code == 3) then
               write(4,'(i7, f20.9, 6f20.9 , i3, 2f20.9)')NAME(i), BODYS(i),(XS(K,i),K=1,3),&
                    &(VS(K,i),K=1,3), KSTAR(i), ZLMSTY(i), RADIUS(i)
            elseif (code == 4) then
               if ( AS(1) < 0.0001) then
                  rrt = 1e20*rt
               else
                  rrt = rt
               endif
               if (r(i) <= rrt) then
                  write(4,'(i7, f20.9, 6f20.9 , i3, 2f20.9)')NAME(i),&
                       & BODYS(i),(XS(K,i),K=1,3),&
                       &(VS(K,i),K=1,3), KSTAR(i), ZLMSTY(i), RADIUS(i)
               else
                  write(8,'(i7, f20.9, 6f20.9 , i3, 2f20.9)')NAME(i),&
                       & BODYS(i),(XS(K,i),K=1,3),&
                       &(VS(K,i),K=1,3), KSTAR(i), ZLMSTY(i), RADIUS(i)
               endif
            elseif (code == 2 .or. code == 5) then
               write(4,'(i7, f20.9, 6f20.9 )')NAME(i), BODYS(i),(XS(K,i),K=1,3),&
               &(VS(K,i),K=1,3)
            endif
         enddo
         if (code == 5) then
            do i=1,NTAIL
               if ( code == 5 ) then
                  write(8,'(i7, f10.3, 6f35.4 )')NAMEt(i), BODYSt(i),(XSt(K,i),K=1,3),&
                  &(VSt(K,i),K=1,3)
               endif
            enddo
         endif
      endif !End of major_output loop

! write to "overview.txt" file
      write(2,'(2f8.1, i7, f11.2, f8.3, 8f9.3)')AS(1), T6, NTOT, mtot, mtot/mtot0,&
      & LR(4), rt, rc, rc/LR(4), AS(13)*Rstar, LR(4)/rt, fbin0, fbint
      if ( mod (int(nint(T6)),tscreen) == 0 ) then
         write(*,*)
         write(*,*) "    T_NB    T_Myr   N        M      M_ratio    Rh       Rt       Rc&
         &     Rc_O_Rh     RC    Rh_O_Rt   fbin0    fbint"
         write(*,'(2f8.1, i7, f11.2, f8.3, 8f9.3)')AS(1), T6, NTOT, mtot, mtot/mtot0,&
         & LR(4), rt, rc, rc/LR(4), AS(13)*Rstar, LR(4)/rt, fbin0, fbint
      endif
      if (debug == 1) then
         write (*,*)"Writing out arrays. Done."
      endif
! write to radii file
      write(7,'(2f9.2, 10f10.5)') AS(1), T6, ( LR(i),i=1,5 ), rt, rc, rc/LR(3), AS(13)*Rstar, LR(3)/rt
      FLUSH(2)
      FLUSH(3)
      FLUSH(7)
      FLUSH(8)
      FLUSH(10)
!************************************
      call kdtree2_destroy(tree)
      if (debug == 1) write (*,*) "Final deallocation "
      deallocate (AS, BODYS, XS, VS, RADIUS , NAME, KSTAR, ZLMSTY, ASS, BODYSS, XSS, VSS,&
      & mydata, list_neighbor_idx, list_neighbor_dis, r)
      deallocate (irho, rho, rho2)

      deallocate (del_r, del_v2, nei_BODYS, nei_NAME, L2)
      if (code == 5) deallocate (ASt, BODYSt, XSt, VSt, NAMEt, ASSt, BODYSSt, XSSt, VSSt)
      print*,"*************************************************************"

      loop_index = loop_index + 1

   ENDDO

   close(8)
!****************************************************************
!****************************************************************
contains

   subroutine termination ( IO, err)
      implicit none
      integer,intent(in)::IO, err
      real*8::start, finish

      if ( err ==0 ) then
         write(*,*)
         write(*,*) ":-)"
         write(*,*) "Normal end of processing of ", input_file
         call cpu_time(finish)
         write(*,*)
         write(*,*)"Time spent was ",finish-start, "sec"
         STOP
      else if ( err == 1 ) then
         write(*,*)
         write(*,*) ":-("
         write(*,*) "Error in reading main arrays in ", input_file
         STOP
      else if ( err == 2 ) then
         write(*,*)
         write(*,*) "Cluster dissolved befor termination time."
         write(*,*)
         write(*,*) ":-)"
         write(*,*) "Normal end of processing of ", input_file
         call cpu_time(finish)
         write(*,*)
         write(*,*)"Time spent was ",finish-start, "sec"
         STOP
      endif
   end subroutine termination

   character(100) function sweep_blanks(in_str)
      character(100), intent(in) :: in_str
      character(100) :: out_str
      character :: ch
      integer :: j

      out_str = " "
      do j=1, len_trim(in_str)
!      get j-th char
         ch = in_str(j:j)
         if (ch .ne. " ") then
            out_str = trim(out_str) // ch
         endif
         sweep_blanks = out_str
      end do
   end function sweep_blanks


   subroutine get_model_name ( arg1, res_init_ntot, res_mtot0 )
      implicit none
      character(len=100) :: arg1
      character(len=10) :: arg2, arg3
      integer :: res_init_ntot
      real*8 :: res_mtot0


      call getarg(1, arg1)
      call getarg(2, arg2)
      call getarg(3, arg3)
      read(arg2,*) res_init_ntot ! Converts arg2 to integer
      read(arg3,*) res_mtot0 ! Convert arg3 to double precision
      write (*,*) arg2, arg3, res_mtot0

   end subroutine get_model_name

end program read_nbody
