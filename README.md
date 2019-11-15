#### Read-nbody  
 read-nbody
 Author: S.Mohammad Hoseini Rad
 smhr313@gmail.com
 Nov 2012, IASBS, Zanjan
 Last modification: 15 November 2019 / 24 Aban 1398
## make:   
```sh
$ make READ-NBODY-TAIL
$ make READ-NBODY-RT
$ make READ-NBODY-2RT
```
## Usage example:
For stars inside tidal radius:
```sh
$ READ-NBODY-RT OUT3 my_simulation 0 0
```
Restarted run example:
```sh
$ READ-NBODY-RT OUT3 my_simulation 150000 100000.0
```
For stars inside two*tidal radius:
```sh
$ READ-NBODY-2RT OUT3 my_simulation 0 0
```
Restarted run example:
```sh
$ READ-NBODY-2RT OUT3 my_simulation 150000 100000.0
```
For stars inside tidal radius and tidal tail:
```sh
$ READ-NBODY-TAIL OUT3 0 0
```
Restarted run example:
```sh
$ READ-NBODY-TAIL OUT3 150000 100000.0
```
In this new version, snapshots are saved in the `snapshot` folder and tidal tail stars are saved in the `tail` folder.
This requires an nbody6 input parameter file with `KZ(23)=3`.
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




