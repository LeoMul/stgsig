!stgsig version 4.0 5/5/2024

!calculates and prints cross sections, collision strengths.

!example input: &stgsig state_to_be_ionized=3 ntran=2 cs_unit=1/
!in file: sigput.

!all input namelist variables explained here:

!gfortran stgsig.f90 -o stgsig-new.x


!cs_unit               - if equal to 0, print cross sections in units of 10-18 cm2
!                        if equal to 1, print cross sections in units of 10-16 cm2
!                        if equal to 2, print cross sections in units of pi*a0**2
!                        if equal to 3, print cross sections in units of a0**2
!                        if equal to 4, plot cross sections in units of 10-21 cm2
!                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.
!                        default value is 0, which it will revert to if it's not set to one of these options i.e > 4 or < -1
!                        I kept all of the units from the previous version of stgsig, should one want them.
!                        note however the numbers are slightly different from stgsig v3

!en_unit               - if equal to 0, print energies in units of Ry
!                        if equal to 1, print energies in units of eV
!                        if equal to 2, print energies in units of cm^{-1}
!                        default value is 0, which it will revert to if it's not set to one of these options i.e > 3 or < 0           





!outputting transitions:

!ntran                 - int*8: number of transitions to be output.
                            !if non-zero - the namelist must be followed by the lower and upper indicies
                            ! of the requested transitions. e.g:

                            !&stgsig ntran=4 cs_unit=0 outputfile='output_mb'
                            !&end 
                            !1 4
                            !1 6
                            !1 8
                            !2 10
!                       note that in using this output mode, it takes the energies directly from omega
!                       and scales them according to the effective charge squared. 
!                       it does NOT also subtract the lower bound energy - the user must do this later.
!                       this is the caveat with outputting multiple transitions.
!                       in principle i could do an energy column for each transition, this would be a bit ugly.

!output_tran           - character*20: name of output file for transitions, default is sg.dat


!ionization stuff:

!state_to_be_ionized   - bound index to ionize out of. if > 0, calls the ionization
!                        routines. <= 0 = ignores (default)
!                      - It continues to vex me that ionization is the correct spelling in academia.
!ion_pot_ryd           - ionization energy (from ground, in Rydbergs) this is required if you use the above.            

!ion_ignore            - number of level indices to be ignored in the ionization sum.
!                        the indices themselves follow the namelist, after the ntran transitions are ented.
!                        NOTE: not! transition indices. rather the final-state index of the transition
!                        state_to_be_ionzied to (index)
!                        i.e, if state_to_be_ionized = 1, and I wish to ignore transitions to levels 3 4 5
!                        my namelist would be followed by 
!                        3 4 5

!output_ion           - character*20: name of output file for ionization, default is sg_ion.dat

!misc:

!debug                 - boolean if = .true. ,only reads in 10 energy pts 
!                                if = .false. (default) reads in full omega.
!                      - this is for testing only, you should just leave it out of your namelist.


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
    integer * 8 :: state_to_be_ionized,cs_unit,en_unit
    real*8,parameter :: pi = 3.141592
    real*8,parameter :: bohr_metres = 5.29177e-11
    real*8,parameter :: atomic_area_m2 = pi * bohr_metres**2

    real*8 ::ion_pot_ryd    

    integer * 8 ::ionindex,ntran
    
    character*20::output_tran,output_ion

    character*20::unit_string,e_unit_string

    namelist/stgsig/ state_to_be_ionized,debug,ntran,cs_unit,en_unit,output_tran,ion_pot_ryd,output_ion

    ntran = 0
    
    !default ionization params.
    state_to_be_ionized = -1
    ion_pot_ryd = -1.0d0

    cs_unit = 1
    debug = .false.

    output_tran = 'sg.dat'
    output_ion = 'sg_ion.dat'

    e_unit_string = 'Ry'
    INQUIRE(file='sigput',EXIST=input_exists)
    INQUIRE(file='OMEGA',EXIST=omega_exists)

    if (.not.(input_exists)) then 
        call print_input()
        stop "no input sigput found, writing example to sigput.ex"
    end if
    if (.not.(omega_exists)) then 
        stop "no omega found"
    end if 

    open(5,file='sigput')
    read(5,stgsig)

    allocate(lower_indices_array(ntran),upper_indices_array(ntran))


    
    !read in omega, get the appropriate data.
    call read_omega(NZED,NELEC,NAST,omega,stat_weight,energies_bound,energies_incident,debug,num_points)
    
    !if outputting transitions, collect and do that
    !can this be moved to a function? who knows.

    if (ntran .gt. 0) then 
        print*,'Outputting ',ntran, 'transitions:'
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
        !once the transitions have been determined
        !collect them
        allocate(output_array(num_points,ntran))
        do i = 1,ntran
            current_transition_index = transition_index(lower_indices_array(i),upper_indices_array(i),nast)
            output_array(:,i) = omega(:,current_transition_index)
        end do 
        
        call convert_units(output_array,cs_unit,lower_indices_array,unit_string)

        call write_out(energies_incident,output_array,output_tran,unit_string,en_unit)

    end if 


    if(state_to_be_ionized .ne. -1) then 
        ionization_csa = get_ionization_collision_strength(state_to_be_ionized,ion_pot_ryd)
        call convert_units(ionization_csa,cs_unit,[state_to_be_ionized],unit_string)
        call write_out(energies_incident - energies_bound(state_to_be_ionized),ionization_csa,output_ion,unit_string,en_unit)
    end if 
    
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
                
                !print*,lower_index,lower_stat_weight,lower_bound_energy

                output_array(:,ii) = output_array(:,ii) / &
                ( real(lower_stat_weight) * (energies_incident - lower_bound_energy))
                !print*,output_array(1,1)
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

    subroutine write_out(energies_incident,output_measure,outputfile,cs_unit_string,en_unit)
        !writes out, i'll eventually make this a pretty ouptut.
        !one must also consider if many transitions are output... might have to edit.
        !assumes output_measure is in collision strength, and we edit the units according 
        !to cs_unit
        real *8 energies_incident(:) , output_measure(:,:)
        real *8 ,allocatable:: energies_incident_output(:)
        integer*8 :: en_unit
        character*200:: transition_header
        character*20::outputfile,cs_unit_string

        integer * 8 :: num_points  , i , cs_unit
        num_points = size(energies_incident)
        energies_incident_output = energies_incident

        SELECT CASE (en_unit)
                   CASE (0) !Ry
                      WRITE(*,*)  "Energies in Rydbergs (Ry)"
                      e_unit_string = "Ry"
                      energies_incident_output = energies_incident
                   CASE (1) !eV
                      WRITE(*,*)  "Energies in eV (eV)"
                      e_unit_string = "eV"
                      energies_incident_output = energies_incident * 13.605693009d0
                   CASE (2)!cm^{-1}
                      WRITE(*,*)  "Energies in wavenumber (cm^{-1}) "
                      e_unit_string = "cm^-1"
                      energies_incident_output = energies_incident * 109677.57d0
                   CASE DEFAULT
                      WRITE(*,*)  "En unit ambigious - Energies in Rydbergs (Ry) "
                      e_unit_string = "Ry"
                      energies_incident_output = energies_incident
                END SELECT


        open(2,file = outputfile)
        write(2,220) e_unit_string,unit_string

        if (ntran .ge. 0) then
            write(2,230) (lower_indices_array(i),upper_indices_array(i),i=1,ntran)
        end if

        do i = 1,num_points 
            write(2,210) energies_incident_output(i) , output_measure (i,:)
        end do 

        close(2)
210     format(200es16.6)
220     format('#e- energy(',a6,')  CS (',a14,')')
230     format('   #Transition: '20(i6,' - ',i6,' '))
    end subroutine


    function get_ionization_collision_strength(index,ion_pot_ryd)
    !calculates the (naive) ionization csa 
    !by summing over all transitions from level 
    !index to the levels above the ionization threshold.
    
    !since i assume you only ionize out of a single level (or term in LS)
    !i just add the omegas first then convert to sigma
    !this saves on a large number of divisions
        real*8,allocatable :: get_ionization_collision_strength(:,:)
        real*8 :: ion_pot_ryd

        integer *8  :: index , trans_index, ion_index
        integer * 8 :: ii
        print* , "---------------------------------------"
        print* , "Calling routine get_ionization_collision_strength()" 
        write(*,110) index
        
        !obligatory checks the user is using my code correctly
        if (index.ge.nast) then 
            stop 'requested level to ionize out of above nast'
        end if 

        allocate(get_ionization_collision_strength(size(energies_incident),1))
        
        !given the ioniztion potential, find the level just above this
        ion_index = ionizing_index(ion_pot_ryd,energies_bound)

        !obligatory checks the user is using my code correctly
        if (index.ge.ion_index) then 
            stop 'requested level to ionize out of above ionization potential'
        end if 


        trans_index = transition_index(index,ion_index,nast)

        write(*,100) trans_index,trans_index+nast-ion_index 

        get_ionization_collision_strength = 0.0d0

        do ii = trans_index,trans_index+ nast-ion_index 
            get_ionization_collision_strength(:,1) = get_ionization_collision_strength(:,1) + omega(:,ii)
        end do


100     format(' adding transisitons from ',i4,' to ',i4)
110     format(" Ionizing out of level ",i4)

        !get_ionization_cross_section(:,1) = omega_2_sigma(get_ionization_cross_section(:,1),int(index))

    end function

    subroutine trapezium_rule(x_array,y_array,integral)
        real*8:: x_array(:), y_array(:), integral
        integer*8 :: ii , N 

        N = size(x_array)
        if ((size(x_array) .ne. size(y_array))) then 
            stop 'failure in routine trapezium_rule() - x and y array not same size'
        end if 

        integral = 0.0d0 
        do ii = 1,N-1
            integral = integral &
            + 0.5d0 * (x_array(ii+1) - x_array(ii)) &
            * (y_array(ii+1)+y_array(ii))
        end do 

    end subroutine 

    subroutine maxwellian_average_collision_strength(electron_temperature_kelvin,&
                                                     maxwellian_av_array,&
                                                     electron_energy_ryd,&
                                                     collision_strength,&
                                                     initial_level_ryd)
        !calculates a maxwellian average.
        !unfinished and untested.
        !the idea was to develop this further for the calculation of 
        !ionization rates as well as bound-bound transitions.
        !recall that bound-free is a bit more different 
        !..

        real*8,intent(inout) :: maxwellian_av_array(:)
        real*8 :: electron_temperature_kelvin(:)
        real*8,allocatable :: electron_temperature_ryd(:)
        
        real*8 ,allocatable:: exponent_array(:)

        real*8 :: electron_energy_ryd(:)
        real*8 :: initial_level_ryd
        real*8 :: collision_strength(:) 
        real*8 :: integral
        !parameters
        real*8,parameter :: ryd_to_ev = 13.606d0
        real*8,parameter :: kel_to_ev = 8.6173303e-5
        real*8,parameter :: fine_struc= 0.0072973525693d0
        
        integer*8 :: ii

        electron_temperature_ryd = electron_temperature_kelvin * kel_to_ev / ryd_to_ev

        do ii = 1,size(electron_temperature_ryd)
            exponent_array = (electron_energy_ryd - initial_level_ryd) / electron_temperature_ryd(ii)

            call trapezium_rule(exponent_array,exp(-exponent_array)*collision_strength,integral)

            maxwellian_av_array(ii) = integral

        end do 

    end subroutine
    
    subroutine print_input()
        open(1,file='sigput.ex')
        write(1,*) '&stgsig' 
        write(1,*) 'state_to_be_ionized=-1'
        write(1,*) 'ion_pot_ryd = -1.0'
        write(1,*) 'ntran=1 '
        write(1,*) 'cs_unit=0 '
        write(1,*) 'en_unit=1 '
        write(1,*) "outputfile='output_mb'"
        write(1,*) '&end '
        write(1,*) '1  2 '
        close(1)
    end subroutine



end program main

