!stgsig version 4.0 5/5/2024

!calculates and prints cross sections, collision strengths.

!example input: &stgsig state_to_be_ionized=3 ntran=2 cs_unit=1/
!in file: sigput.

!all input namelist variables explained here:

!gfortran stgsig.f90 -o stgsig-new.x

!cs_unit               - if equal to  0, print cross sections in units of 10-18 cm2
!                        if equal to  1, print cross sections in units of 10-16 cm2
!                        if equal to  2, print cross sections in units of pi*a0**2
!                        if equal to  3, print cross sections in units of a0**2
!                        if equal to  4, prubt cross sections in units of 10-21 cm2
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

!ion_include           - number of level indices to be specified in the ionization sum.
!                        the indices themselves follow the namelist, after the ntran transitions are entered.
!                        NOTE: not! transition indices. rather the final-state index of the transition
!                        state_to_be_ionzied to (index)
!                        e.g::
!                         1 0
!                         2 1
!                         3 0
!                         4 1

!                        i.e, if state_to_be_ionized = 1, and I wish to include transitions to levels 4 5
!                        sin the sum.
!                      - otherwise! all transitions above the ion threshold are included.


!output_ion           - character*20: name of output file for ionization, default is sg_ion.dat

!misc:

!debug                 - boolean if = .true. ,only reads in 10 energy pts 
!                                if = .false. (default) reads in full omega.
!                      - this is for testing only, you should just leave it out of your namelist.


!development notes
!i might for fun implement upsilon also

! originally by XXXXXX, 
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

program stgsig24
    implicit none 

    !booleans for user error
    logical               :: input_exists, omega_exists
    
    !OMEGA data
    integer*8             :: NZED , NELEC
    integer*8             :: NAST , NUM_POINTS
    integer*8,allocatable :: stat_weight(:)
    real*8   ,allocatable :: energies_bound(:)
    real*8   ,allocatable :: energies_incident(:)    
    real*8   ,allocatable :: omega(:,:)
    

    integer*8,parameter   :: num_steps = 10
    !Input and input adjacent:

    !Debugging:
    logical               :: debug
    logical               :: istop

    !Ionization
    integer*8             :: state_to_be_ionized
    real*8                :: ion_pot_ryd    
    integer*8             :: ion_include
    character*20          :: output_ion

    !Outputting transitions
    character*20          :: omega_path
    character*20          :: output_tran
    integer*8             :: ntran
    integer*8,allocatable :: lower_indices_array(:),upper_indices_array(:)
    !Units 
    integer*8             :: cs_unit,en_unit
    character*20          :: unit_string,e_unit_string

    !namelist
    namelist/stgsig/ state_to_be_ionized,debug,ntran,cs_unit,en_unit,output_tran,ion_pot_ryd,output_ion,ion_include,omega_path
    !default namelist arguments
    ntran = 0
    state_to_be_ionized = -1
    ion_pot_ryd = -1.0d0
    cs_unit = 1
    ion_include = 0
    debug = .false.
    omega_path = 'OMEGA'
    output_tran = 'sg.dat'
    output_ion = 'sg_ion.dat'
    e_unit_string = 'Ry'

    INQUIRE(file='sigput',EXIST=input_exists)

    INQUIRE(file=omega_path,EXIST=omega_exists)

    istop = .False.
    if (.not.(input_exists)) then 
        call print_input()
        istop = .True.
        print*, "no input sigput found, writing example to sigput.ex"
    end if
    if (.not.(omega_exists)) then 
        istop = .True.
        print*, "no omega found"
    end if 

    if (istop) then 
        stop 'Stopping abnormally'
    end if

    open(5,file='sigput')
    read(5,stgsig)

    call read_omega(NZED,NELEC,NAST,omega,stat_weight,energies_bound,energies_incident,debug,num_points,omega_path)

    call outputting_transitions(NTRAN,STAT_WEIGHT,energies_bound,energies_incident)
    print*,'-------------------------------------'


    
    call ionization_cross_section(state_to_be_ionized,&
                                  ion_include,&
                                  nast,&
                                  ion_pot_ryd,&
                                  cs_unit,&
                                  en_unit,&
                                  stat_weight,&
                                  energies_bound,&
                                  energies_incident)
    print*,'-------------------------------------'
    print*,'STOPPING NORMALLY'
    
    contains 

    subroutine read_omega(NZED, &
                          NELEC,&
                          NAST ,&
                          OMEGA,&
                          STAT_WEIGHT,&
                          energies_bound,&
                          energies_incident,&
                          debug,&
                          num_points,&
                          omega_path)

        !reads OMEGA, funny that.
        !returns, charge, num electrons,
        !omega, stat weight
        !Bound and incident energies in Rydbergs. 
        !the charge scaling is corrected in this routine 
        !so it's in real Rydbergs.

            !Parameters determined from OMEGA
            integer*8,intent(inout) :: NZED , NELEC
            integer*8               :: NAST , NUM_POINTS,NUM_TRAN
            real*8                  :: charge_scaling_energy

            !Allocatables 
            integer*8,allocatable :: bigL(:),mult(:)
            integer*8,allocatable :: stat_weight(:)
            real*8   ,allocatable :: energies_bound(:)
            real*8   ,allocatable :: omega(:,:)
            real*8   ,allocatable :: energies_incident(:)

            !Input options
            logical      :: debug
            character*20 :: omega_path

            !Iterators
            integer*8 :: i
            integer*8 :: j
    
            !Read in times
            real*8 :: read_time_1,read_time_2
    
            print*,'-------------------------------------'
            print*,'Entering routine read_omega()'
    
            call cpu_time(read_time_1)
    
            open(1,file=omega_path,status='old',form='formatted')

            print*,'reading file,' ,omega_path

            !reading in basic data.
            !this allows allocataion of energies, OMEGA

            read(1,*) NZED,NELEC 
            read(1,*) NAST,NUM_POINTS,NUM_TRAN 
    
            print*,'found nuclear charge: ',NZED
            print*,'found num electrons : ',NELEC
            print*,'found num bound stat: ',NAST
            print*,'found num mesh point: ',NUM_POINTS
            print*,'found num transition: ',NUM_TRAN
    
            if (debug) then 
                print*,'DEBUG IS ON, ONLY READING A FEW POINTS'
                num_points = 1
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
            !should probably add an option for unformatted
            !but i think it's as simple as changing the open statement
     
            do i = 1, num_points
                read(1,*) energies_incident(i),(omega(i,j),j=1,NUM_TRAN)
            end do     
            close(1)
    
            call cpu_time(read_time_2)
    
            write(*,100) read_time_2-read_time_1

            !Stat weight = (2L+1)*(2S+1)
            !I think in fine structure this is also taken care of, as multiplicity
            !is also set to zero in OMEGA - and L-> J in OMEGA
            stat_weight = (2*bigL + 1) *abs( mult)
            
            !Determining the charge scaling, so scaled Rydbergs never harms another
            !my hairline again.
            if ((NZED-NELEC) .gt. 1) then 
                charge_scaling_energy = (NZED-NELEC)**2
                energies_bound = energies_bound * charge_scaling_energy
                energies_incident = energies_incident *charge_scaling_energy
            end if

            print*,'-------------------------------------'
            100 format (' OMEGA read time: ',F5.1,' seconds.')
    
        end subroutine

        subroutine outputting_transitions(NTRAN,STAT_WEIGHT,energies_bound,energies_incident)
            
            !inputs
            integer*8    :: NTRAN
            integer*8    :: STAT_WEIGHT(:)
            real*8       :: energies_bound(:)
            real*8       :: energies_incident(:)

            integer*8    :: current_transition_index
            character*20 :: unit_string

            !Alloctables
            integer*8,allocatable    :: lower_indices_array(:)
            integer*8,allocatable    :: upper_indices_array(:)
            real*8, allocatable      :: output_array(:,:)

            !iterators - this is awful coding
            integer*8:: i,ii,jj
            

            allocate(lower_indices_array(ntran),upper_indices_array(ntran))

            if (ntran .lt. 1) then 
                print*,'ntran less than 1 - skipping transition output'
                return  
            end if 
            print*,'Outputting ',abs(ntran), 'transitions:'

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
            
            call convert_units(output_array,cs_unit,lower_indices_array,unit_string,stat_weight,energies_bound,energies_incident)
    
            call write_out(energies_incident,output_array,output_tran,&
                    unit_string,en_unit,lower_indices_array,upper_indices_array,ntran)


        end subroutine

        subroutine ionization_cross_section(state_to_be_ionized,&
                                            ion_include,&
                                            nast,&
                                            ion_pot_ryd,&
                                            cs_unit,&
                                            en_unit,&
                                            stat_weight,&
                                            energies_bound,&
                                            energies_incident)

            !Outputs the ionization cross section.

                !done here:
            !Step 0: determine states with which transitions are to be included in the summation.
            !Step 1: determine ionizing index

                !done in subrotuines
            !Step 2: perform summation
            !Step 3: convert units
            !Step 4: ouput

            !inputs:
            integer*8 :: state_to_be_ionized
            integer*8 :: ion_include
            integer*8:: nast
            real*8  :: ion_pot_ryd

            integer*8 :: cs_unit,en_unit
            integer*8 :: stat_weight(:)
            real*8    :: energies_bound(:)
            real*8    :: energies_incident(:)
            
            !allocatables
            integer*8,allocatable :: ion_included_levels(:)
            real*8,allocatable :: ionization_collision_strength(:,:)
            
            !indicators:
            integer :: levelnum, pseudo_indicator
            
            !parameters:
            integer*8 :: index_of_first_unbound 
            integer*8,parameter :: zero = 0 
            !iterators
            integer*8 :: i

            !state_to_be_ionized = 0 or less is the indicator to skip this routine.
                if (state_to_be_ionized .lt. 1) then 
                    print*, 'state_to_be_ionized less than 1 - exiting ionization_cross_section()'
                    return 
                end if 
            
            print* , "Entering routine ionization_cross_section()" 
    
            !obligatory checks the user is using my code correctly
                if (state_to_be_ionized.ge.nast) then 
                    stop 'requested level to ionize out of above nast'
                end if 
            
            !determine where we are ionizing from.
                index_of_first_unbound = ionizing_index(ion_pot_ryd,energies_bound)
            
            !is the requested bound state above the ionization limit?
            !obligatory checks the user is using my code correctly
                if (state_to_be_ionized.ge.index_of_first_unbound) then 
                    stop 'requested level to ionize out of above ionization potential'
                end if 


            !determine how we are performing the summation.

                allocate(ion_included_levels(nast))

                !if ion include > 0, then we are restricting the ionization summatuon.
                if (ion_include.gt.0) then 
                    ion_included_levels = 0
                
                    do i = 1,ion_include
                        read(5,*) levelnum,pseudo_indicator
                        ion_included_levels(i) = pseudo_indicator
                    end do  

                else
                    !else, all the levels are included.
                    ion_included_levels = 1
                end if

            !calculate ionization:
                call calculate_ionization(ion_included_levels,&
                                          state_to_be_ionized,&
                                          index_of_first_unbound,&
                                          nast,&
                                          ionization_collision_strength,&
                                          energies_incident)

                call convert_units(ionization_collision_strength,&
                                   cs_unit,&
                                   [state_to_be_ionized],&
                                   unit_string,&
                                   stat_weight,&
                                   energies_bound,&
                                   energies_incident)

                call write_out(energies_incident - energies_bound(state_to_be_ionized),& 
                               ionization_collision_strength,&
                               output_ion,&
                               unit_string,&
                               en_unit,&
                               lower_indices_array,&
                               upper_indices_array,&
                               zero)
                               

            
        end subroutine

        subroutine calculate_ionization(ion_included_levels,&
                                        state_to_be_ionized,&
                                        index_of_first_unbound,&
                                        nast,&
                                        ionization_collision_strength,&
                                        energies_incident)
            
            !inputs
            integer*8 :: ion_included_levels(:)
            integer*8 :: state_to_be_ionized
            integer*8 :: index_of_first_unbound
            integer*8 :: nast
            real*8    :: energies_incident(:)

            !allocatables
            real*8,allocatable :: ionization_collision_strength(:,:)

            !parameters:
            integer*8 :: level_index
            integer*8 :: trans_index

            allocate(ionization_collision_strength(size(energies_incident),1))
            ionization_collision_strength  = 0.0d0

            !todo: make this a proper fortran format.
            write(*,100) index_of_first_unbound

            write(*,110) state_to_be_ionized

            if(all(ion_included_levels .eq. 1)) then 
                write(*,120)
            else 
                write(*,130)
            end if 

            trans_index =  transition_index(state_to_be_ionized,index_of_first_unbound,nast)

            !loop through all the levels, and keep track of the transition index.
            !if the current level isnt included in the sum, move on to the next one.
            do level_index = index_of_first_unbound,nast

                if (ion_included_levels(level_index) .eq. 1) then
                    ionization_collision_strength(:,1) = &
                    ionization_collision_strength(:,1) + omega(:,trans_index)
                end if 

                trans_index = trans_index + 1
            end do

            100     format(" First unbound state is:  ",i4,'.')
            110     format(" Ionizing out of level ",i4)
            120     format(" Including all terms in summation.")
            130     format(" Excluding some terms in summation.")

        end subroutine
        
        function transition_index(lower_index,upper_index,nast)
            !converts lower and upper index unto transition number
            !the inverse of this function is tbd, not needed yet idk

            integer * 8 :: upper_index,lower_index,nast 
            integer * 8 :: transition_index

            transition_index = upper_index - lower_index &
                             + (lower_index-1) * (2*nast - lower_index) /2 

        end function
    
        function ionizing_index(ionization_potential,bound_energies)
            !locates the term index of the first free term
            real * 8 :: ionization_potential,bound_energies(:)
            integer * 8 :: i ,ionizing_index
            do i = 1,size(bound_energies)
                if (bound_energies(i) > ionization_potential) then 
                    exit
                end if 
            end do 
            ionizing_index = i
        end function

        subroutine convert_units(output_array,&
                                 cs_unit,&
                                 lower_indices_array,&
                                 unit_string,&
                                 stat_weight,&
                                 energies_bound,&
                                 energies_incident)
            !if equal to 0, print cross sections in units of 10-18 cm2
        !                        if equal to 1, print cross sections in units of 10-16 cm2
        !                        if equal to 2, print cross sections in units of pi*a0**2
        !                        if equal to 3, print cross sections in units of a0**2
        !                        if equal to 4, plot cross sections in units of 10-21 cm2
        !                        if equal to -1, print COLLISION STRENGTHS, i.e do not convert.
                
                !Constants, physical parameters
                real*8,parameter      :: pi = 3.141592d0
                real*8,parameter      :: bohr_metres = 5.29177e-11
                real*8,parameter      :: atomic_area_m2 = pi * bohr_metres**2

                !Inputs:
                real*8    :: energies_bound(:)
                real*8    :: energies_incident(:)
                integer*8 :: stat_weight(:)
                integer*8 :: cs_unit,lower_indices_array(:)
                !Outputs:
                character*20 ::unit_string
                real*8,intent(inout) :: output_array(:,:)

                !Iterators and temp variables
                integer*8 :: ii ,lower_index
                integer * 8 :: lower_stat_weight 
                real * 8 :: lower_bound_energy 

                print*,'-------------------------------------'

                print*, 'Entering routine convert_units()'
        
                !if cs_unit is -1, do nothing
                if (cs_unit .ne. -1) then 
            
                    !for each transition, convert to cross-section.
                    !this naturally involves dividing by 
                    !energy of incident electron rel to 
                    !lower energy state of transition.
        
                    do ii = 1,size(lower_indices_array)
                        !print*,'converting col',ii
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

        subroutine write_out(energies_incident,&
                             output_measure,&
                             outputfile,&
                             cs_unit_string,&
                             en_unit,&
                             lower_indices_array,&
                             upper_indices_array,&
                             ntran)

            !writes out, i'll eventually make this a pretty ouptut.
            !one must also consider if many transitions are output... might have to edit.
            !assumes output_measure is in collision strength, and we edit the units according 
            !to cs_unit

            integer*8 :: lower_indices_array(:)
            integer*8 :: upper_indices_array(:)
            integer*8 :: ntran
            real *8 energies_incident(:) , output_measure(:,:)
            real *8 ,allocatable:: energies_incident_output(:)
            integer*8 :: en_unit
            character*20::outputfile,cs_unit_string
    
            integer * 8 :: num_points  , i 
            num_points = size(energies_incident)
            energies_incident_output = energies_incident
            
            !this should be moved inside the convert units routine.
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
            write(2,220) e_unit_string,cs_unit_string
    
            if (ntran .gt. 0) then
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



        subroutine print_input()
            open(1,file='sigput.ex')
            write(1,*) '&stgsig' 
            write(1,*) 'state_to_be_ionized=-1'
            write(1,*) 'ion_pot_ryd = -1.0'
            write(1,*) 'ntran=0 '
            write(1,*) 'cs_unit=0 '
            write(1,*) 'en_unit=0 '
            write(1,*) "output_tran='output_mb'"
            write(1,*) '&end '
            close(1)
        end subroutine

end program stgsig24