!*****************************************************************************
module quantity
!*****************************************************************************
!
!  Program:    phq
!
!  Module:     quantity
!
!  Purpose:    This module defines all the physical quantities.

!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = 20 ), dimension ( 10 ) :: element

   integer :: n_atoms
   integer :: n_atom1
   integer :: n_species
   integer :: n_step           !  Total md steps 
   integer :: pole             !  Number of poles for maximum entropy
   integer :: n_step_use       !  Number of MD steps used in frouier transform and maximum entropy when getting frequency 
   integer :: n_step_initial   !  The first n_step_initial MD will not be output by read_md.f90
   integer :: n_corr
   integer :: n_window
   integer :: method
   integer :: super_size
   integer, dimension(3) :: supercell
   integer, dimension(3) :: super

   double precision :: lambda         ! adiabatic and reversible scaling lambda
   double precision :: lattice_parameter
   double precision :: d_t                              ! MD time length
   double precision :: temperaturemd                    ! MD temperature
   double precision :: mass(10)
   double precision, dimension(3,3) :: celldm
   double precision, dimension(3,3) :: recip
   double precision, allocatable :: alat_md_position(:,:,:)
   double precision, allocatable :: alat_position(:,:)
   double precision, allocatable :: ammode (:)
   double precision, allocatable :: atom_mass(:)
   double precision, allocatable :: cartesian_md_position(:,:,:)
   double precision, allocatable :: cartesian_position(:,:)
   double precision, allocatable :: displacement ( :,:,:)
   double precision, allocatable :: kinetic_energy(:)
   double precision, allocatable :: md_force (:,:,:) 
   double precision, allocatable :: omega_corr (:)
   double precision, allocatable :: omega_corr_fit (:)
   double precision, allocatable :: omega_lorentzian (:)
   double precision, allocatable :: omega_mem (:)
   double precision, allocatable :: omega_phonon (:)
   double precision, allocatable :: primitive_cell (:,:)
   double precision, allocatable :: primitive_position (:,:)
   double precision, allocatable :: q_point (:,:)
   double precision, allocatable :: r_point (:,:)
   double precision, allocatable :: real_vector (:,:)
   double precision, allocatable :: tau_fourier (:)
   double precision, allocatable :: tau_mem (:)
   double precision, allocatable :: taumode (:)
   double precision, allocatable :: temperature(:)
   double precision, allocatable :: total_energy(:) 

   complex ( kind = kind( 0.0d0 ) ), allocatable :: eigen_vector(:,:)
   complex ( kind = kind( 0.0d0 ) ), allocatable :: temp_vector(:,:)
   complex ( kind = kind( 0.0d0 ) ), allocatable :: vector_q (:,:,:)

!*****************************************************************************

end module quantity


!*****************************************************************************
module parameter
!*****************************************************************************
!
!  Program:    phq
!
!  Module:     parameter
!
!  Purpose:    This module contains all the parameters used.
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************

! numerical constants

   double precision, parameter :: zero = 0.0d0
   double precision, parameter :: one = 1.0d0
   double precision, parameter :: two = 2.0d0
   double precision, parameter :: three = 3.0d0
   double precision, parameter :: four = 4.0d0
   double precision, parameter :: five = 5.0d0
   double precision, parameter :: six = 6.0d0
   double precision, parameter :: seven = 7.0d0
   double precision, parameter :: eight = 8.0d0
   double precision, parameter :: nine = 9.0d0
   double precision, parameter :: ten = 10.0d0
   double precision, parameter :: eleven = 11.0d0
   double precision, parameter :: twelve = 12.0d0
   double precision, parameter :: half = 0.5d0
   double precision, parameter :: fourth = 0.25d0
   double precision, parameter :: threequaters = 0.75d0
   double precision, parameter :: threehalf = 1.5d0
   double precision, parameter :: sixth = 1.0d0 / 6.0d0
   double precision, parameter :: pi = 3.1415927d0
   double precision, parameter :: thertz = 1.0d12

! physical constants

   double precision, parameter :: adu = 0.52917715d-10
   double precision, parameter :: amu = 1.660538921d-27
   double precision, parameter :: ryd = 13.60569253 * 1.60217646d-19
   double precision, parameter :: atu = 4.8378d-17
   double precision, parameter :: e_mass = 5.4858d-4   !  unit = amu
   double precision, parameter :: n_emass = 911.44424213227d0
   double precision, parameter :: atm = 2.418884326505d-17
   double precision, parameter :: boltzmann_k = 8.617332478d-5   ! ev K-1  
   double precision, parameter :: H_planck_SI = 6.62606896E-34   ! J s
   double precision, parameter :: Hartree_Si = 4.35974394E-18    ! J 

! unit conversion

   double precision, parameter :: rad_to_deg = 57.29578d0
   double precision, parameter :: au_to_joules = 4.3598d-18
   double precision, parameter :: GPa_to_au = 3.398923d-5
   double precision, parameter :: s_to_atu = 2.418884326d-17 
   double precision, parameter :: au_thz = H_planck_SI / Hartree_Si / two / pi * 1.0d12                                            
   double precision, parameter :: ry_to_thz = 1.0d0 / au_thz / four / pi
   double precision, parameter :: thz_to_cm = 33.35641d0

!*****************************************************************************

end module parameter


!*****************************************************************************
module text
!*****************************************************************************
!
!  Program:    phq
!
!  Module:     text
!
!  Purpose:    This module contains all the text processing methods.
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  subroutines internal to this module

   private compress
   private low_case
   private string_to_integer
   private string_to_real

contains

!*****************************************************************************
subroutine compress( line )
!*****************************************************************************
!
!  Purpose:    This subroutine substitutes multiple blanks for a single blank.
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = 80 ), intent (INOUT) :: line

!*****************************************************************************
!  local variables

   integer :: i

!*****************************************************************************
!  begin subroutine 

   do 

      i = index( trim( line ), "  " )
      if ( i == 0 ) exit
      line (i:) = line(i+1:)

   end do

end subroutine compress

!*****************************************************************************
subroutine low_case( line )
!*****************************************************************************
!
!  Purpose:    This subroutine converts all upper to lower case.
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = 80 ), intent (INOUT) :: line

!*****************************************************************************
!  local variables

   integer :: i, length, n

!*****************************************************************************
!  begin subroutine 

   length = len( line )

   do i=1, length
      n = ichar( line(i:i) )
      if ( ( n >= 65 ) .and. ( n <= 91 ) ) n = n + 32
      line(i:i) = char( n )
   end do

end subroutine low_case

!*****************************************************************************
subroutine split_line( line, n_commands, commands,                           &
                             n_integers, integer_numbers,                    &
                             n_reals, real_numbers )
!*****************************************************************************
!
!  Purpose:    This subroutine splits a line of text into its constituent words, 
!              and these are classified as commands, integer numbers and reals.
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = 80 ), intent (INOUT) :: line
   character ( len = 30 ), intent (OUT) :: commands(10)

   integer, intent (OUT) :: n_commands
   integer, intent (OUT) :: n_integers
   integer, intent (OUT) :: n_reals

   integer, intent (OUT) :: integer_numbers(10)

   double precision, intent(OUT) :: real_numbers(10)

!*****************************************************************************
!  local variables

   integer :: end_of_word
   integer :: i
   integer :: n_words

   character ( len = 40 ) :: words( 30 )

   character ( len = 13 ), parameter :: integer_set = " +-0123456789"
   character ( len = 16 ), parameter :: real_set = " +-de.0123456789"

!*****************************************************************************
!  begin subroutine

   n_words = 0

   line = adjustl( line )

! get rid of any multiple-blanks in the text

   call compress( line )

! convert to lower case

   call low_case( line )

! now split line into words

   if ( len_trim( line ) == len( line ) ) then

      n_words = n_words + 1
      words(n_words) = line

   else

      do 

         if ( len_trim( line ) == 0 ) exit

         end_of_word = index( line, " " ) - 1
         n_words = n_words + 1
         words(n_words) = line(1:end_of_word)
         line = line(end_of_word+2: )
     
      end do

   end if

! now check the nature of each word

   n_integers = 0
   n_reals = 0
   n_commands = 0

   commands = (/'empty','empty','empty','empty','empty',                  &
                'empty','empty','empty','empty','empty'/)

   do i=1, n_words
        
      if ( verify( words(i), integer_set ) == 0 ) then
        
         n_integers = n_integers + 1
         integer_numbers( n_integers ) = string_to_integer( words(i) )

      else if ( verify( words(i), real_set ) == 0 ) then

         n_reals = n_reals + 1 
         real_numbers( n_reals ) = string_to_real( words(i) )

      else
   
         n_commands = n_commands + 1
         commands( n_commands ) = words(i) 

      end if

   end do


   if ( n_commands == 0 ) then

      commands(1) = "empty"

   end if

end subroutine split_line

!*****************************************************************************
function string_to_integer( string ) result( integer_number )
!*****************************************************************************
!
!  Purpose:    This subroutine converts the string representation of an integer 
!              number to the actual integer number. 
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = * ) :: string
   integer :: integer_number

!*****************************************************************************
!  local variables

!*****************************************************************************
!  begin subroutine 

   if ( len( string ) > 0 ) then
      
      read( string, * ) integer_number

   else

      print *, "Sorry, this is not an integer: ", string
      stop

   end if

end function string_to_integer

!*****************************************************************************
function string_to_real( string ) result( real_number )
!*****************************************************************************
!
!  Purpose:    This subroutine converts the string representation of a real 
!              number to the actual integer number. 
!
!*****************************************************************************
   
   implicit none

!*****************************************************************************
!  shared variables

   character ( len = * ) :: string
   double precision :: real_number

!*****************************************************************************
!  local variables

   integer :: status

!*****************************************************************************
!  begin subroutine 

   if ( len( string ) > 0 ) then
      
      read( string, *, iostat = status ) real_number

   else

      print *, "Sorry, this is not a real: ", string
      stop

   end if

   if ( status < 0 ) then
 
      print *, "Sorry, this is not a real: ", string
      stop
 
   end if

end function string_to_real

end module text

