#gfortran ./kdtree2.f90 ./read_mod.f90 ./read_nbody.f90 -o read_nbody.exe
gfortran -O3 -ffast-math ./kdtree2.f90 ./read_nbody_st.f90 -o read_nbody_st.exe