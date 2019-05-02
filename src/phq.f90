!**************************************************************************
program phq    
!**************************************************************************
!
!     Program:     phq
!                       
!     Purpose:     This program performs phonon quasiparticle calculations from 
!                  harmonic phonon calculations and molecular dynamics simulations.
!                                              
!**************************************************************************
!  used modules
  
   use main
   use parameter
   use text
   
!**************************************************************************

   implicit none

!**************************************************************************
!  Local variables

   logical :: debug

   character*8   :: date
   character*10  :: time
   character*5   :: zone
   character ( len = 80 ) :: line
   character ( len = 30 ), dimension ( 10 ) :: words
   character ( len = 30 ), dimension ( 7 ), parameter ::                 &
       command = (/                             &
       'dt                            ',     & ! 1
       'step_md_use                   ',     & ! 2
       'correlation_time              ',     & ! 3
       'pole                          ',     & ! 4
       'supercell                     ',     & ! 5
       'temperature                   ',     & ! 6
       'method                        '/)      ! 7

  integer       :: value1(8)
  integer       :: value2(8)
  integer :: i, j, n, n_line, status  
  integer :: n_integers
  integer :: n_words
  integer :: n_reals  
  integer :: integer_numbers(10)
  
  double precision :: real_numbers(10)

!**************************************************************************
!  Start of subroutine

   call DATE_AND_TIME(date, time, zone, value1)

!
!  read from input file
!
   method = 0

   open(unit=5)
   
   n_line = 0

   do 
    
      read(5, "(a80)", iostat = status ) line
  
      if ( status < 0 ) exit
     
      n_line = n_line + 1
    
      if ( line(1:1) == '#' ) cycle
    
      call split_line( line, n_words, words,                              &
            n_integers, integer_numbers,                                  &
            n_reals, real_numbers )

      if (words(1) == command(1)) then

           d_t = real_numbers (1)  

      else if (words(1) == command(2) )  then

           n_step_use = integer_numbers (1)

      else if (words(1) == command(3) )  then

           n_corr =  integer_numbers (1)
           n_window = n_corr

      else if (words(1) == command(4) )  then
    
           pole = integer_numbers (1) 

      else if (words(1) == command(5) )  then

           super(1) = integer_numbers (1)
           super(2) = integer_numbers (2)
           super(3) = integer_numbers (3)

           super_size = super(1) * super(2) * super(3)

           allocate( r_point (super_size,3) )   ! in cartesian coordinate
           allocate( q_point (super_size,3) )   ! in cartesian coordinate

      else if (words(1) == command(6) )  then
           temperaturemd = real_numbers (1)

      else if (words(1) == command(7) )  then
           method = integer_numbers (1)

      end if

   end do

!
!  end of reading input
!
   write (6,*) "Reading input"
   write (6,*) "dt= ", d_t 

   call control ( debug )
 
   call DATE_AND_TIME(date, time, zone, value2)
   write(6,*) "Simulation time =", value2(5)-value1(5), "H",         &
                                   value2(6)-value1(6), "M",         &
                                   value2(5)-value1(7), "S"

   deallocate ( r_point )
   deallocate ( q_point )
   

end program phq
