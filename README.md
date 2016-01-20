#### Read-nbody  
 read_nbody (st version)  
 Author: S.Mohammad Hoseini Rad  
 smhr313@gmail.com  
 Nov 2012, IASBS, Zanjan  
 Last modification: 3 September 2015 / 12 Shahrivar 1394
## Compile:   
```sh
$ compile.sh
```
## Usage example:
Common run:  
```sh
$ read-nbody-st.exe OUT3 my_simualtion 0 0
```
Restarted run: 
```sh
$ read-nbody-st.exe OUT3 my_simualtion 150000 100000.0
```
===============================

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
===============================

## Changes:

29 Mar 2013: For solving incompatibility problem with old gfortran compilers (e.g 4.1.2), I re-inputted subroutines of read-mod.f90 module into main file (read-nbody.f90).

8 Apr 2013: Calculation of tidal radius was added. Now just stars within tidal radius are considered for calculation of total number of stars and other main specifications of the cluster.

10 Apr 2013: Directory creation with "model_name" name was added. Now, previous "All-" files are made in the "model_name" directory with the same prefix.

Optimization options "-O3" and "-ffast-math" were added to "compile.sh" file.

30 June 2013: "tout" option was added. It is an option for time interval of output on screen & harddisk.   
It Also prevent neighbor arrays and kdtree2 pointers to be allocated. So if you have many time snapshots, by increasing 'tout', you can pass the memory overflow problem. (NBODY6 custom: Myr, NBODY6: N-body unit).  

24 September 2013: Testing "NAME" of stars and inclusion of them within tidal radius were separated. Procedure of testing "NAME" was corrected. The code now is able to calculate core radius (Trenti et al. 2007). The code now is also able to find binaries (experimental).

15 February 2014 / 26 Bahman 1392: Minor change in "overview.txt" header to coincide with my python script which creates "all_summary.txt" file.

17 February 2015 / 28 Bahman 1393: "Add call kdtree2_destroy(tree)".

3 September 2015 / 12 Shahrivar 1394: Add two input arguments "res_INIT_NTOT, res_mtot0" to deal with restarted simulations.
20 January 2016 / 30 Dey 1394: Correct ecc2 formula
