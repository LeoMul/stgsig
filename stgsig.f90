! stgsig version 4.0 5/5/2024

! calculates and prints cross sections, collision strengths.

!example input: &stgsig state_to_be_ionized=3 ntran=2 cs_unit=1/
!in file: sigput.

!all input namelist variables explained here:

!cs_unit               - if equal to 0, print cross sections in units of 10-18 cm2
!                        if equal to 1, print cross sections in units of 10-16 cm2
!                        if equal to 2, print cross sections in units of pi*a0**2
!                        if equal to 3, print cross sections in units of a0**2
!                        if equal to 4, plot cross sections in units of 10-21 cm2
!                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.
!                        default value is 0.

!state_to_be_ionized   - bound index to ionize out of. if non-zero, calls the ionization
!                        routines. zero = ignores (default)

!debug                 - boolean if = .true. ,only reads in 10 energy pts 
!                                if = .false. (default) reads in full omega.
!development notes
!i might for fun implement upsilon also

! originally by XXX
! updated by Leo Patrick Mulholland
! with the added advantage of ionization cross sec


!goal is to make it all dynamically allocated
!with the removal of common blocks and gotos
!i.e - bring the code into 2024

!additionally - goal is to have an option to 
!directly calculate ionization cross section
!this has been implemented

!right now only reads formatted OMEGA file
!should ask Connor how to approach unformatted
!additionally - assumes no elastic transitions

program main 

    !I hope these variables are self-explanatory.
    implicit none

    integer :: i
    integer :: ii,jj
    logical :: debug, input_exists, omega_exists
    integer*8 :: NZED , NELEC
    integer*8 :: NAST , NUM_POINTS,NUM_TRAN
    integer*8,allocatable :: stat_weight(:)
    real*8,allocatable :: energies_bound(:)
    real*8,allocatable :: omega(:,:)
    real*8,allocatable :: energies_incident(:)
        
    real*8,allocatable :: omega_tran(:)
    real*8,allocatable :: sigma_tran(:)
    
    !array to be output - num energies x num transitions 
    real*8 ,allocatable :: output_array(:,:)

    real*8 ,allocatable :: ionization_csa(:,:)


    integer * 8,allocatable :: lower_indices_array(:),upper_indices_array(:)
    integer*8 ::current_transition_index
    integer * 8 :: state_to_be_ionized,cs_unit
    real*8,parameter :: pi = 3.141592
    real*8,parameter :: bohr_metres = 5.29177e-11
    real*8,parameter :: atomic_area_m2 = pi * bohr_metres**2

    real*8 ::units
    real*8,parameter :: ip = 0.66    
    integer * 8 ::ionindex,ntran,state_output
    
    character*20::outputfile

    character*20::unit_string,e_unit_string

    namelist/stgsig/ state_to_be_ionized,debug,ntran,cs_unit,state_output,outputfile

    ntran = 0
    state_to_be_ionized = -1
    cs_unit = 1
    units=1.0d0
    debug = .false.
    outputfile = 'sg.dat'
    e_unit_string = 'Ry'
    INQUIRE(file='sigput',EXIST=input_exists)
    INQUIRE(file='OMEGA',EXIST=omega_exists)

    if (.not.(input_exists)) then 
        stop "no input sigput found"
    end if
    if (.not.(omega_exists)) then 
        stop "no omega found"
    end if 

    open(5,file='sigput')
    read(5,stgsig)

    allocate(lower_indices_array(ntran),upper_indices_array(ntran))


    
    !read in omega, get the appropriate data.
    call read_omega(NZED,NELEC,NAST,omega,stat_weight,energies_bound,energies_incident,debug,num_points)
    
    if (ntran .ge. 0) then 
        do i = 1,ntran
            read(5,*) ii,jj
            lower_indices_array(i) = min(ii,jj)
            upper_indices_array(i) = max(ii,jj)
            if (ii.eq.jj) then 
                print*,'one of the requested transitions is between same level '
                stop
            end if 
            if ((ii.ge.nast).or.(jj.ge.nast)) then 
                print*,'one of the requested transitions has a level above max level'
                stop
            end if 

            if ((ii.lt.1).or.(jj.lt.1)) then 
                print*,'one of the requested transitions has a negative level index'
                stop
            end if 

        end do 
    end if 

    if(state_to_be_ionized .ne. -1) then 
        ionization_csa = ionization_cross_section(state_to_be_ionized)
    end if 
    
    allocate(output_array(num_points,ntran))

    do i = 1,ntran
        current_transition_index = transition_index(lower_indices_array(i),upper_indices_array(i),nast)
        output_array(:,i) = omega(:,current_transition_index)
    end do 
    
    !allocate(lower_indices_array(2))
    !lower_indices_array(1) = 1
    !lower_indices_array(2) = 1

    call convert_units(output_array,cs_unit,lower_indices_array,unit_string)

    call write_out(energies_incident,output_array,outputfile,unit_string,e_unit_string)

    contains 

    subroutine read_omega(NZED,NELEC,NAST,OMEGA,STAT_WEIGHT,energies_bound,energies_incident,debug,num_points)
    !reads OMEGA, funny that.
    !returns, charge, num electrons,
    !omega, stat weight
    !Bound and incident energies in Rydbergs. 
    !the charge scaling is corrected in this routine 
    !so it's in real Rydbergs.
        integer*8,intent(inout) :: NZED , NELEC
        integer*8 :: NAST , NUM_POINTS,NUM_TRAN
        integer*8,allocatable :: bigL(:),mult(:)
        integer*8,allocatable :: stat_weight(:)
        logical :: debug
        real*8,allocatable :: energies_bound(:)
        real*8,allocatable :: omega(:,:)
        real*8,allocatable :: energies_incident(:)

        real*8 :: charge_scaling_energy

        integer*8 :: i
        integer*8 :: j


        real*8 :: read_time_1,read_time_2

        print*,'-------------------------------------'
        print*,'Entering routine read_omega()'

        call cpu_time(read_time_1)

        open(1,file='OMEGA',status='old',form='formatted')

        !reading in basic data.
        !this allows allocataion of energies, OMEGA
        read(1,*) NZED,NELEC 
        read(1,*) NAST , NUM_POINTS,NUM_TRAN 

        print*,'found nuclear charge: ',NZED
        print*,'found num electrons : ',NELEC
        print*,'found num bound stat: ',NAST
        print*,'found num mesh point: ',NUM_POINTS
        print*,'found num transition: ',NUM_TRAN

        if (debug) then 
            print*,'DEBUG IS ON, ONLY READING A FEW POINTS'
            num_points = 10
        end if


        !allocating arrays.
        allocate(energies_bound(NAST))  
        allocate(bigL(NAST))  
        allocate(mult(NAST))  

        allocate(energies_incident(NUM_POINTS))  
        allocate(omega(NUM_POINTS,NUM_TRAN))  


        !reading in data.
        read(1,*) (mult(i),bigL(i),i=1,nast)
        read(1,*) (energies_bound(i),i=1,nast)

        !for now I will read in the whole omega file
        !in principle one could just read in the requested transitions...

        !regular omega format
 

        do i = 1, num_points
            read(1,*) energies_incident(i),(omega(i,j),j=1,NUM_TRAN)
        end do     
        close(1)

      

        call cpu_time(read_time_2)

        !I should update this to an actual fortran format.
        print*, 'omega read time ', read_time_2-read_time_1 ,'sec'

        stat_weight = (2*bigL + 1) *abs( mult)

        if ((NZED-NELEC) .gt. 1) then 
            charge_scaling_energy = (NZED-NELEC)**2
            energies_bound = energies_bound * charge_scaling_energy
            energies_incident = energies_incident *charge_scaling_energy
        end if
        print*,'-------------------------------------'

    end subroutine

    function transition_index(lower_index,upper_index,nast)
        !converts lower and upper index unto transition number
        !the inverse of this function is tbd
        integer * 8 :: upper_index,lower_index,nast 
        integer * 8 :: transition_index
        transition_index = upper_index - lower_index &
        + int(0.5*(lower_index-1) * (2*nast - lower_index))
    end function

    function ionizing_index(ionization_potential,bound_energies)
        real * 8 :: ionization_potential,bound_energies(:)
        integer * 8 :: i ,ionizing_index
        
        do i = 1,size(bound_energies)
            if (bound_energies(i) > ionization_potential) then 
                exit
            end if 
        end do 
        ionizing_index = i

    end function

    subroutine extract_transition(omega_tran,omega,lower_index,upper_index,nast)
        !extracts a particular transition
        real * 8 :: omega(:,:)
        real * 8 ,allocatable:: omega_tran(:)
        
        integer,intent(in) :: lower_index,upper_index
        integer*8 :: nast
        integer *8 :: transition_index 
        !this might need checked - might be minus ii. infact it probably is
        transition_index = (lower_index - 1) * nast + upper_index-1
        omega_tran = omega(:,transition_index)


    end subroutine

    subroutine convert_units(output_array,cs_unit,lower_indices_array,unit_string)
    !if equal to 0, print cross sections in units of 10-18 cm2
!                        if equal to 1, print cross sections in units of 10-16 cm2
!                        if equal to 2, print cross sections in units of pi*a0**2
!                        if equal to 3, print cross sections in units of a0**2
!                        if equal to 4, plot cross sections in units of 10-21 cm2
!                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.
        real*8,intent(inout) :: output_array(:,:)
        integer*8 :: cs_unit,lower_indices_array(:)
        real*8,allocatable :: output_again(:,:)
        integer*8 :: ii ,lower_index
        integer * 8 :: lower_stat_weight 
        real * 8 :: lower_bound_energy 
        character*20 ::unit_string
        print*, 'Entering routine convert_units()'

        !if cs_unit is -1, do nothing
        if (cs_unit .ne. -1) then 
    
            !for each transition, convert to cross-section.
            !this naturally involves dividing by 
            !energy of incident electron rel to 
            !lower energy state of transition.

            do ii = 1,size(lower_indices_array)
                lower_index        = lower_indices_array(ii)
                lower_stat_weight  = stat_weight(lower_index)
                lower_bound_energy = energies_bound(lower_index)
                

                output_array(:,ii) = output_array(:,ii) / &
                ( real(lower_stat_weight) * (energies_incident - lower_bound_energy))

            end do 
        end if 

        !if equal to 0, print cross sections in units of 10-18 cm2
!                        if equal to 1, print cross sections in units of 10-16 cm2
!                        if equal to 2, print cross sections in units of pi*a0**2
!                        if equal to 3, print cross sections in units of a0**2
!                        if equal to 4, plot cross sections in units of 10-21 cm2
!                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.

                SELECT CASE (cs_unit)
                     CASE (-1) !collision strength
                      unit_string = "dimensionless"
                      WRITE(*,*)  "collision strength"
                      !output_array(:,ii) = output_array(:,ii)
                   CASE (0) !10-18 cm2
                      WRITE(*,*)  "cross sec in units of 10^{-18} cm^2 "
                      unit_string = "10^{-18} cm^2"
                      output_array(:,:) = output_array(:,:)* atomic_area_m2 * 1e22
                   CASE (1)!10-16 cm2
                      WRITE(*,*)  "cross sec in units of 10^{-16} cm^2 "
                      unit_string = "10^{-16} cm^2"

                      output_array(:,:) = output_array(:,:)* atomic_area_m2 * 1e20
                   CASE (2) !pi*a0^2 - which is the atomic unit of area, we are in already.
                      unit_string = "pi a_0^2"
                      output_array(:,:) = output_array(:,:)
                   CASE (3)!a^0 ^2 - so just multiply prev by pi
                      unit_string = "a_0^2"
                      output_array(:,:) = output_array(:,:)* pi * 1e22
                   CASE (4)
                      unit_string = "10^{-21} cm^2"
                      output_array(:,:) = output_array(:,:)* atomic_area_m2 * 1e25
                   CASE DEFAULT
                      unit_string = "10^{-18} cm^2"
                      WRITE(*,*)  "CS unit ambigious - cross sec in units of 10^{-18} cm^2 "
                      output_array(:,:) = output_array(:,:)* atomic_area_m2 * 1e22
                END SELECT

    end subroutine

    function omega_2_sigma(omega_tran,lower_index)
        !converts a single transition omega to sigma in units.
        real * 8 :: omega_tran(:)
        real * 8, allocatable :: omega_2_sigma(:)
        integer :: lower_index

        integer * 8 :: lower_stat_weight 
        real * 8 :: lower_bound_energy 

        real * 8 :: unit  
        unit = 1.0d0
        allocate(omega_2_sigma(size(omega_tran)))
        lower_stat_weight = stat_weight( lower_index)
        lower_bound_energy = energies_bound(lower_index)
        omega_2_sigma = omega_tran /( real(lower_stat_weight) * (energies_incident - lower_bound_energy))
        omega_2_sigma = omega_2_sigma* atomic_area_m2 * 1e22

    end function

!if equal to 0, print cross sections in units of 10-18 cm2
!                        if equal to 1, print cross sections in units of 10-16 cm2
!                        if equal to 2, print cross sections in units of pi*a0**2
!                        if equal to 3, print cross sections in units of a0**2
!                        if equal to 4, plot cross sections in units of 10-21 cm2
!                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.

    subroutine write_out(energies_incident,output_measure,outputfile,cs_unit_string,e_unit_string)
        !writes out, i'll eventually make this a pretty ouptut.
        !one must also consider if many transitions are output... might have to edit.
        !assumes output_measure is in collision strength, and we edit the units according 
        !to cs_unit
        real *8 energies_incident(:) , output_measure(:,:)
        character*200:: transition_header
        character*20::outputfile,cs_unit_string,e_unit_string

        integer * 8 :: num_points  , i , cs_unit
        num_points = size(energies_incident)


        open(2,file = outputfile)
        write(2,220) e_unit_string,unit_string

        if (ntran .ge. 0) then
            write(2,230) (lower_indices_array(i),upper_indices_array(i),i=1,ntran)
        end if

        do i = 1,num_points 
            write(2,210) energies_incident(i) , output_measure (i,:)
        end do 

        close(2)
210     format(200es16.6)

220     format('#e- energy(',a4,')  CS (',a14,')')
230     format('   #Transition: '20(i6,' - ',i6,' '))
    end subroutine


    function ionization_cross_section(index)
    !calculates the (naive) ionization csa 
    !by summing over all transitions from level 
    !index to the levels above the ionization threshold.
    
    !since i assume you only ionize out of a single level (or term in LS)
    !i just add the omegas first then convert to sigma
    !this saves on a large number of divisions
        real *8,allocatable :: ionization_cross_section(:,:)
        integer *8  :: index , trans_index, ion_index
        integer * 8 :: ii

        print* , "Calling routine ionization_cross_section" 
        print* , "Ionizing out of level",index

        allocate(ionization_cross_section(size(energies_incident),1))

        ion_index = ionizing_index(ip,energies_bound)

        !print*, 'found index',ion_index,nast
        trans_index = transition_index(index,ion_index,nast)

        !print*, 'found index',trans_index
        print*, 'adding trans from',trans_index,'to',trans_index+nast-ion_index 
        ionization_cross_section = 0.0d0
        do ii = trans_index,trans_index+ nast-ion_index 
            ionization_cross_section(:,1) = ionization_cross_section(:,1) + omega(:,ii)
        end do
        print*, 'omega added, calculating csa'

        
    
        ionization_cross_section(:,1) = omega_2_sigma(ionization_cross_section(:,1),int(index))

    end function

end program main

