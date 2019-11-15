#gfortran ./kdtree2.f90 ./read_mod.f90 ./read_nbody.f90 -o read_nbody.exe
#gfortran -O3 -ffast-math ./kdtree2.f90 ./read-nbody-st.f90 -o read-nbody-st.exe
#gfortran -O3 -ffast-math -g -fbounds-check -ffree-line-length-none -Wall ./kdtree2.f90 ./read-nbody-st.f90 -o read-nbody
gfortran -O3 -ffast-math -g -fbounds-check -ffree-line-length-none -Wall ./kdtree2.f90 ./read-nbody-tail.f90 -o read-nbody-tail
