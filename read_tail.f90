program read_tail

   implicit none
   ! ******************************************
   real*8::start, finish ! For cpu_time subroutine.
   integer*4::i,j,l,k,ll ! loop counter variables.
   integer*4::IO,code, err, debug
   integer::tscreen ! Frequency of output on screen.
   integer::tout ! Frequency of output on harddisk.
   integer,dimension(13)::buff ! This related to measuring size of input file.
   integer*4::status ! This is related to measuring size of input file.
   integer::loop_index
   integer::INIT_NTOT
   ! ******************************************
   ! NBODY6 outputs (here inputs). See output.f in Ncode directory of NBODY6,
   ! before removing stars with NAME < 1 .
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
   real*4,dimension(:),allocatable::BODYSS,ASS  ! For converting single precision to double precision ( if code = 2).
   ! ******************************************
   character(len=40)::model_name, dir_command
   ! *****************************************
   real*8::Mstar, Tstar, Rstar, Vstar
   real,parameter::G = 0.0043009211 ! In pc*km^2/s^2*M_sun
   real,parameter::AU = 206264.806 ! pc to AU
   ! ******************************************
   integer:: res_INIT_NTOT
   real*8:: res_mtot0
   ! ******************************************

   call cpu_time(start)
