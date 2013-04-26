!**************************************************************
! read_nbody (st ver)
! Author: S.Mohammad Hoseini Rad
! smhr313@gmail.com
! Nov 2012, IASBS, Zanjan
! Last modification: 10 Apr 2013
!**************************************************************

To compile, please use "compile.sh" script.

AS arrays members according to "output.f" file in "Ncode" directory (Aarseth's NBODY 6).

      AS(1) = TTOT
      AS(2) = FLOAT(NPAIRS)
      AS(3) = RBAR
      AS(4) = ZMBAR
      AS(5) = RTIDE
      AS(6) = TIDAL(4)
      AS(7) = RDENS(1)
      AS(8) = RDENS(2)
      AS(9) = RDENS(3)
      AS(10) = TTOT/TCR
      AS(11) = TSCALE
      AS(12) = VSTAR
      AS(13) = RC
      AS(14) = NC
      AS(15) = VC
      AS(16) = RHOM
      AS(17) = CMAX
      AS(18) = RSCALE
      AS(19) = RSMIN
      AS(20) = DMIN1
***************************************************************

Changes:
***************
29 Mar 2013: For solving incompatibility problem with old gfortran compilers (e.g 4.1.2),
I re-inputted subroutines of read-mod.f90 module into main file (read-nbody.f90).
***************
8 Apr 2013: Calculation of tidal radius was added. Now just stars within tidal radius are considered 
for calculation of total number of stars and other main specifications of the cluster.
***************
10 Apr 2013: Directory creation with "model_name" name was added. Now, previous "All-" files
are made in the "model_name" directory with the same prefix.

Optimization options "-O3" and "-ffast-math" were added to "compile.sh" file.
