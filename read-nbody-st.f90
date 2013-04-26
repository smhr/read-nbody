!**************************************************************
! read_nbody (st version)
! Author: S.Mohammad Hoseini Rad
! smhr313@gmail.com
! Nov 2012, IASBS, Zanjan
! Last modification: 10 Apr 2013
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
! use read_mod
use kdtree2_module
implicit none
! ******************************************
real*8::start, finish ! For cpu_time subroutine.
integer*4::i,j,l,k,ll ! loop counter variables.
integer*4::IO,code, err, debug
integer::tout
integer,dimension(13)::buff ! This related to measuring size of input file.
integer*4::status ! This is related to measuring size of input file.
integer::loop_index 
integer::INIT_NTOT
! ******************************************
! NBODY6 outputs (here inputs). See output.f in Ncode directory of NBODY6,
! before removing stars with NAME < 1 and outside tidal rdius.
integer*4::iNTOT
integer,dimension(:),allocatable::iNAME,iKSTAR
real*8,dimension(:),allocatable::iBODYS,iRADIUS,iZLMSTY
real*8,dimension(:,:),allocatable::iXS,iVS
! ******************************************
! NBODY6 outputs (here inputs), after removing stars with NAME < 1 and outside tidal rdius.
integer*4::NTOT,NK,N,MODEL,NRUN
integer,dimension(:),allocatable::NAME,KSTAR
real*8,dimension(:),allocatable::AS,BODYS,RADIUS,ZLMSTY
real*8,dimension(:,:),allocatable::XS,VS
real*4,dimension(:,:),allocatable::XSS,VSS
real*4,dimension(:),allocatable::BODYSS,ASS  ! For converting single precesion to double precesion ( if code = 2). 
! ******************************************
real*8::mtot,imtot
integer,dimension(:,:),allocatable::list_neighbor_idx ! Neighbor index array.
real,dimension(:,:),allocatable::list_neighbor_dis ! Neighbor distance array.
character(len=60)::output_file, input_file, fprefix
integer::n_neighbor ! Number of nearest neighbors in Casterano & Hut method.
real*8,dimension(3)::Xd ! Density center coordinates.
real*8,dimension(5)::LR ! Lagrangian radii array.
real*8,dimension(:),allocatable::r,ir ! Position vector of stars after and befor tidal radius detemination.
! ******************************************
type(kdtree2), pointer::tree
real, allocatable:: mydata(:,:)
type(kdtree2_result), allocatable :: results(:)
! ******************************************
real*8::sixth_neighbor_distance, m_tot_5th, rho, rho_tot
real*8,dimension(3)::rho_times_dis_tot
! ******************************************
real*8::mtemp, rtemp, rstep
integer::check, diss_check
! ******************************************
real::Ndiss, Mdiss
! ******************************************
character(len=50)::model_name, dir_command

 call cpu_time(start)

input_file='OUT3' ! Input file.
open(1,file=input_file, form='unformatted')
open(2,file='overview.txt') ! This file includes some cluster characteristics.
open(3,file='density_center.txt') ! This file includes density center coordinate of cluster during evolution.
open(7,file='radii.txt') ! This file includes Lagrangian radii.

! ******************************************
!** Initialization of internal variables ***
! Please don't change this variable.
IO = 0;
err = 0
status = 0
n_neighbor = 6 ! Number of nearest neighbors in Casterano & Hut (1985) and Aarseth (2001) method.
loop_index = 0
INIT_NTOT = 0; iNTOT = 0; imtot = 0
rstep = 0.01 ! Bin length in astronomical unit for calculating lagrangian radii.

! ************ Code options ****************
! ******************************************
 code = 2  	! 1:NBODY6 custom, 2: NBODY6
 tout = 1		! Time interval of output on screen. (NBODY6 custom: Myr, NBODY6: N-body unit)
 debug = 0		! 1: Debug mode.
 diss_check = 1		! 1: Check dissolution of cluster; 0: No check.
 Ndiss = 0.01		! By reaching to this fraction of initial number of stars, the cluster is considered as a dissolved cluster; 0: Suppress this option.
 Mdiss = 0		! By reaching to this fraction of initial mass, the cluster is considered as a dissolved cluster; 0: Suppress this option.
 ! Please note just use one of Ndiss or Mdiss options.
 model_name = 'N100000-RG6.0-Rh6.2'
! ******************************************
 dir_command = 'mkdir '// model_name 
 call system ( dir_command )

 call STAT (input_file, buff,status) ! Calculating size of input file.
 if (status == 0) then
	write(*,*)
	write(*,*)"Size of file is", buff(8), "Byte."
 endif
! ******************************************
write(7,*)" T(NBODY)  Time(Myr)   LR_0.01  LR_0.05  LR_0.50  LR_0.75  LR_0.95"
write(2,*) "      T(NBODY)   T(Myr)   NTOT(in Rt)    MTOT(in Rt)       Rh            Rt"

! ************ Read in loop ****************
! *******************************************
! do l=1,10 

! By activating this loop, you can define exact time in terms of output interval (DELTAT in NBODY6 or dtout in mcluster [N-body units (Myr in Nbody6 custom)]).


DO WHILE ( .TRUE. ) ! Repeating the loop until the end of input file.
	
! *******************************************
! define_arrays_size ( code, input_file, loop_index, NTOT, MODEL, NRUN, NK, N, INIT_NTOT ) ! Read number of particles in each loop.
if ( code == 1 ) then
	read(1, iostat=IO)iNTOT,NK,N
	if (debug == 1 ) write (*,'(3i10)') iNTOT,NK,N
	if ( IO /= 0 ) then
		call termination ( IO, err)
	endif
else if ( code == 2 ) then
	read(1, iostat=IO)iNTOT, MODEL, NRUN, NK
	if ( loop_index == 0 ) INIT_NTOT = iNTOT
	if (debug == 1 ) write (*,'(5i10)') iNTOT, MODEL, NRUN, NK
	if ( IO /= 0 ) then
		call termination ( IO, err)
	endif
endif
! *******************************************
! This condition is for appropriate termination of the program if cluster dissolves.
if ( iNTOT > 0 ) then
	if ( diss_check == 1 ) then 
		if ( Ndiss > 0  ) then
! 			print*,iNTOT , Ndiss , INIT_NTOT, Ndiss * INIT_NTOT, iNTOT - Ndiss * INIT_NTOT
			if ( iNTOT <= Ndiss * INIT_NTOT ) then
				err = 2
				call termination ( IO, err )
			endif
		endif
		if ( loop_index == 1) imtot = mtot
		if ( imtot > 0 ) then
! 			 print*,Ndiss,Mdiss
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
	write (*,*) "iNTOT <= 0. Progrom is stoped!"
	STOP
endif
			
!***********************************
! 	call readin ( code, input_file, iNTOT, NK, AS, iBODYS, iXS, iVS, iRADIUS, iNAME, iKSTAR, iZLMSTY ) ! Read input arrays according to the output.f in Ncode directory of NBODY6. Note that there are some differences between NBODY6 custom and NBODY6 input variables.
!************************************
allocate (AS(NK)); allocate (iBODYS(iNTOT)); allocate(iXS(3,iNTOT)); allocate(iVS(3,iNTOT))
allocate(iRADIUS(iNTOT)); allocate(iNAME(iNTOT)); allocate(iKSTAR(iNTOT)); allocate(iZLMSTY(iNTOT))

allocate (ASS(NK)); allocate (BODYSS(iNTOT)); allocate(XSS(3,iNTOT)); allocate(VSS(3,iNTOT))

if (debug == 1 ) print*,"Allocating readin arrays."

if ( code == 1 ) then

	read(1, iostat=IO)(AS(K),K=1,NK),(iBODYS(J),J=1,iNTOT),((iXS(K,J),K=1,3),J=1,iNTOT),((iVS(K,J),K=1,3),J=1,iNTOT),&
		(iRADIUS(J),J=1,iNTOT),(iNAME(J),J=1,iNTOT),&
		(iKSTAR(J),J=1,iNTOT),(iZLMSTY(J),J=1,iNTOT)
	if (debug == 1 ) write (*,*)(k,AS(K),K=1,NK)
	if ( IO /= 0 ) then
		err = 1
		call termination ( IO, err)
	endif

else if ( code == 2 ) then

	read(1, iostat=IO)(ASS(K),K=1,NK),(BODYSS(J),J=1,iNTOT),((XSS(K,J),K=1,3),J=1,iNTOT),((VSS(K,J),K=1,3),J=1,iNTOT),&
		(iNAME(J),J=1,iNTOT)
	if (debug == 1 ) write (*,*)(k,ASS(K),K=1,NK)

	do J=1, NK	! Convert to double precesion
		AS(J) = ASS(J)
	enddo

	do J=1, iNTOT	! Convert to double precesion
		iBODYS(J) = BODYSS(J)
		do K=1, 3
			iXS(K,J) = XSS(K,J)
			iVS(K,J) = VSS(K,J)
		enddo
	enddo

	if ( IO /= 0 ) then
		err = 1
		call termination ( IO, err)
	endif
	
endif
!************************************
!	write(*,'(i10, f11.1, i 10)')NTOT,AS(10),N
!************************************
! 	call neighbor_find ( n_neighbor, iNTOT, iXS, list_neighbor_idx, list_neighbor_dis ) ! finding all neighbors and their distances for all stars and saving them in list_neighbor_idx and list_neighbor_dis arrays correspondingly.
!************************************
allocate(mydata(3,iNTOT))
mydata = iXS ! putting x,y and z of all stars in mydata array to change them to single precision. 
tree => kdtree2_create(mydata,rearrange=.true.,sort=.true.)
allocate (list_neighbor_idx (iNTOT,n_neighbor + 1))
allocate (list_neighbor_dis (iNTOT,n_neighbor + 1))

do i=1,iNTOT ! Find star neighbors and their distances to the ith star in XSt array.
	allocate (results(n_neighbor + 1)) ! The results(1) is associated with the ith star. The results(2..n_neighbor+1) are associated with n_neighbor nearest neighbors. So for 6 nearest neighbors, size of results array must be 7.

	call kdtree2_n_nearest_around_point(tree,idxin=i,nn=n_neighbor + 1,correltime=0,results=results)
	do ll=1,n_neighbor + 1

		list_neighbor_idx (i,ll) = results (ll)%idx
		list_neighbor_dis (i,ll) = results (ll)%dis
	enddo
	deallocate (results)
enddo
!************************************
	
! 	call center_find ( n_neighbor, iNTOT, iBODYS, iXS, list_neighbor_idx, list_neighbor_dis, Xdt )  ! Finding density center of cluster according to the Casterano & Hut (1985) and Aarseth 2001 method. For more details see center_find subroutine.
!************************************
rho_tot=0; rho_times_dis_tot=0; Xd=0
! tree => kdtree2_create(mydata,rearrange=.true.,sort=.true.)
do i=1,iNTOT
	m_tot_5th=0; rho=0; sixth_neighbor_distance=0
	sixth_neighbor_distance=sqrt(list_neighbor_dis (i,n_neighbor + 1))
	
! 	write(*,*)results(7)%idx, sixth_neighbor_distance
	do ll=2, n_neighbor
		m_tot_5th=m_tot_5th+iBODYS(list_neighbor_idx (i,ll)) ! Total mass of fifth nearest neibors.
	enddo
	rho=m_tot_5th/(sixth_neighbor_distance**3.)
	rho_tot = rho_tot + rho
! 	print*,m_tot_5th, sixth_neighbor_distance ,rho, rho_tot
	rho_times_dis_tot(1) = rho_times_dis_tot(1) + rho * iXS(1,i)
	rho_times_dis_tot(2) = rho_times_dis_tot(2) + rho * iXS(2,i)
	rho_times_dis_tot(3) = rho_times_dis_tot(3) + rho * iXS(3,i)
! 	write(*,*)i, m_tot_5th
	
enddo
Xd(1) = ( rho_times_dis_tot(1)/rho_tot )
Xd(2) = ( rho_times_dis_tot(2)/rho_tot )
Xd(3) = ( rho_times_dis_tot(3)/rho_tot )
write(3,*) (Xd(i),i=1,3)
!***********************************
! 	
! 	call coordinate_correction ( iNTOT, Xdt, iXS ) ! Correcting the origin of coordinate system by subtracting all star coordinates from measured density center.
!***********************************
do j = 1, iNTOT
	do i = 1, 3
		iXS(i,j) = iXS(i,j) - Xd(i)
	enddo
enddo

!***********************************
! call test_name_and_rt ( code, iNTOT, AS, iBODYS, iXS, iVS, iRADIUS, iNAME, iKSTAR, iZLMSTY,&
! 			     & NTOT, BODYS, XS, VS, RADIUS, NAME, KSTAR, ZLMSTY, r )
!print*,intot,ntot,"intot,ntot"
i = 0; k = 0; NTOT = 0
!print*,intot,ntot,"intot,ntot"
allocate(ir(iNTOT))
if (debug == 1 ) write (*,*) "Allocating ir Array in test_name_and_rt subroutine."
do i = 1, iNTOT
		if ( iNAME(i) >= 1) then
			ir(i) = sqrt( iXS(1,i)*iXS(1,i) + iXS(2,i)*iXS(2,i) + iXS(3,i)*iXS(3,i) )
			if ( code == 1 ) then
				if ( ir(i) <= AS(25) ) then
					k = k + 1
				endif
			elseif ( code == 2 ) then
				if ( ir(i) <= AS(5) ) then
					k = k + 1
				endif
			endif
		endif
enddo

NTOT = k
!print*,ntot,"nnnntot"

allocate (BODYS(NTOT)); allocate(XS(3,NTOT)); allocate(VS(3,NTOT))
allocate(RADIUS(NTOT)); allocate(NAME(NTOT)); allocate(KSTAR(NTOT)); allocate(ZLMSTY(NTOT))
allocate(r(NTOT))

if (debug == 1 ) write (*,*)"Allocating test_name_and_rt arrays."
k = 0
do i = 1, iNTOT
		if ( iNAME(i) >= 1) then
			ir(i) = sqrt( iXS(1,i)*iXS(1,i) + iXS(2,i)*iXS(2,i) + iXS(3,i)*iXS(3,i) )
			if ( code == 1 ) then
				if ( ir(i) <= AS(25) ) then
					k = k + 1
					r(k) = ir(i)
					BODYS(k) = iBODYS(i)
					XS(1,k) = iXS(1,i); XS(2,k) = iXS(2,i); XS(3,k) = iXS(3,i)
					VS(1,k) = iVS(1,i); VS(2,k) = iVS(2,i); VS(3,k) = iVS(3,i)
					RADIUS(k) = iRADIUS(i)
					NAME(k) = iNAME(i)
					KSTAR(k) = iKSTAR(i)
					ZLMSTY(k) = iZLMSTY(i)
				endif
			elseif ( code == 2 ) then
				if ( ir(i) <= AS(5) ) then
					k = k + 1
					r(k) = ir(i)
					BODYS(k) = iBODYS(i)
					XS(1,k) = iXS(1,i); XS(2,k) = iXS(2,i); XS(3,k) = iXS(3,i)
					VS(1,k) = iVS(1,i); VS(2,k) = iVS(2,i); VS(3,k) = iVS(3,i)
					RADIUS(k) = iRADIUS(i)
					NAME(k) = iNAME(i)
				endif
			endif
		endif
enddo



!***********************************
! 	
! 	call radii ( NTOT, BODYS, XS, AS, mtot, LR ) ! Computing Lagrangian radii.
!************************************
rtemp = 0; mtemp = 0; check = 1
mtot = 0
do i = 1, NTOT
	BODYS(i) = BODYS(i) * AS(4)
	mtot = mtot + BODYS(i)
	r(i)=r(i)*AS(3)
	do k = 1, 3
		XS(k,i)=XS(k,i)*AS(3); VS(k,i)=VS(K,i)*AS(12)
	enddo
enddo; i = 0
!  print*,mtot
!allocate(r(NTOT))

do while ( check /= 6 )
	mtemp = 0
	rtemp = rtemp + rstep
 	
	do i = 1, NTOT
		!r(i) = sqrt( XS(1,i)*XS(1,i) + XS(2,i)*XS(2,i) + XS(3,i)*XS(3,i) )
		if ( r(i) <= rtemp ) then 
			mtemp = mtemp + BODYS(i)
!  			print*,mtemp
			if  ( mtemp >= mtot * 0.01 .AND. check == 1)  then
				LR(1) = rtemp; !print*,"LR1 = ",LR(1)
				check = 2; exit
			elseif ( mtemp >= ( mtot * 0.05 ) .AND. check == 2 ) then
				LR(2) = rtemp; !print*,"LR5 = ",LR(2)
				check = 3; exit
			elseif ( mtemp >= ( mtot * 0.50 ) .AND. check == 3 ) then
				LR(3) = rtemp; !print*,"LR50 = ",LR(3)
				check = 4; exit
			elseif ( mtemp >= ( mtot * 0.75 ) .AND. check == 4 ) then
				LR(4) = rtemp; !print*,"LR75 = ",LR(4)
				check = 5; exit
			elseif ( mtemp >= mtot * 0.95 .AND. check == 5 ) then
				LR(5) = rtemp; !print*,"LR95 = ",LR(5)
				check = 6; exit
			endif
		endif
	enddo
enddo

if (debug == 1 ) write (*,*)"Calculating LRs. Done."
if  ( code == 1 ) write(7,'(2f9.2, 5f10.5)') AS(1), AS(10), ( LR(i),i=1,5 )
if  ( code == 2 ) write(7,'(2f9.2, 5f10.5)') AS(1), AS(11)*AS(1), ( LR(i),i=1,5 )


!************************************
! 	
! 	call writeout ( code, output_file, mtot, LR, tout, NTOT, NK, AS, BODYS, XS, VS, RADIUS, NAME, KSTAR, ZLMSTY ) ! Writing out outputs.
!************************************
if ( code == 1 ) then

	Write(output_file , '( f7.1 )' )AS(1)
	! write (*,*)nint(AS(10))
	if (nint(AS(10))<10) then
		output_file = trim(model_name) // '-0000' // trim(output_file) // '.txt'
! 		output_file = output_file // ".txt"
! 		print*,output_file,len_trim(output_file)
		else if (nint(AS(10))<100) then
			output_file = trim(model_name) // '-000' // trim(output_file) // '.txt'
		else if (nint(AS(10))<1000) then
			output_file = trim(model_name) // '00' // trim(output_file) // '.txt'
		else if (nint(AS(10))<10000) then
			output_file = trim(model_name) // '0' // trim(output_file) // '.txt'
		else 
		output_file = trim(model_name) // trim(output_file) // '.txt'
	endif

elseif ( code == 2 ) then

	Write(output_file , '( f7.1 )' )AS(1)
	! write (*,*)nint(AS(10))
	if (nint(AS(1))<10) then
		output_file = trim(model_name) //'-0000'// trim(output_file) //".txt"
! 		output_file = output_file // ".txt"
! 		print*,output_file,len_trim(output_file)
		else if (nint(AS(1))<100) then
			output_file = trim(model_name) // '-000' // trim(output_file) // '.txt'
		else if (nint(AS(1))<1000) then
			output_file = trim(model_name) // '00' // trim(output_file) // '.txt'
		else if (nint(AS(1))<10000) then
			output_file = trim(model_name) // '0' // trim(output_file) // '.txt'
		else 
		output_file = trim(model_name) // trim(output_file) // '.txt'
	endif
endif

output_file = trim(model_name) // '/' // trim(output_file)
output_file = sweep_blanks(output_file)
open(4,file=output_file)
! TOT_MASS=0
do i=1,NTOT
! 	BODYS(i)=BODYS(i)*AS(4)
! 	do k=1,3
! 		XS(k,i)=XS(K,i)*AS(3); VS(k,i)=VS(K,i)*AS(12)
! 	enddo
! 	TOT_MASS=TOT_MASS+BODYS(i)
	if ( code == 1 ) then
		write(4,'(i6, f13.9, 6f13.6 , i3, 2f15.5)')NAME(i), BODYS(i),(XS(K,i),K=1,3),&
		&(VS(K,i),K=1,3), KSTAR(i), ZLMSTY(i), RADIUS(i)
	elseif ( code == 2 ) then
		write(4,'(i6, f13.9, 6f13.6 )')NAME(i), BODYS(i),(XS(K,i),K=1,3),&
		&(VS(K,i),K=1,3)
	endif
enddo
	if ( code == 1 ) then
		write(2,'(2f11.1, i12, 3f15.5)') AS(1), AS(10), NTOT, mtot, LR(3), AS(25)*AS(3)
		if ( mod (int(nint(AS(10))),tout) == 0 ) then 
			write(*,*) 
			write(*,*) "Time(NBODY)   Time(Myr)   Total number    Total mass       Rh            Rt"
			write(*,'(2f11.1, i12, 3f15.2)')AS(1), AS(10), NTOT, mtot, LR(3), AS(25)*AS(3)
		endif
	elseif ( code == 2 ) then
		write(2,'(2f11.1, i12, 3f15.5)')AS(1), AS(11)*AS(1), NTOT, mtot, LR(3), AS(5)*AS(3)
! 		AS(11)*AS(1)=TSCALE*TTOT, 
		if ( mod (int(nint(AS(1))),tout) == 0 ) then 
			write(*,*) 
			write(*,*) "Time(NBODY)   Time(Myr)   Total number    Total mass       Rh            Rt"
			write(*,'(2f11.1, i12, 3f15.2)')AS(1), AS(11)*AS(1), NTOT, mtot, LR(3), AS(5)*AS(3)
		endif
	endif
	if (debug == 1 ) then 
	write (*,*)"Writing out arrays. Done."
	print*,"****************************"
endif
!************************************

deallocate ( AS,BODYS, XS, VS, RADIUS , NAME, KSTAR, ZLMSTY, ASS, BODYSS, XSS, VSS,&
           & mydata, list_neighbor_idx, list_neighbor_dis, r )
deallocate  (  iBODYS, iXS, iVS, iRADIUS , iNAME, iKSTAR, iZLMSTY, ir)

loop_index = loop_index + 1

ENDDO
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

character(60) function sweep_blanks(in_str)
character(60), intent(in) :: in_str
character(60) :: out_str
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


end program read_nbody
