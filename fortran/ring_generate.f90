program ring_generate

  use params
  use octree
  use energy_terms
  use ring_polymer
  use ring_objects
  use ring_grow
  use read_write
  use omp_lib

  implicit none

  type(ring) :: ring_init, ring_rep
  type(obstacle_distribution) :: obstacle_distribution_init, obstacle_distribution_rep
  type(grow_params) :: grow_params_sim
  type(deform_params) :: deform_params_sim
  type(energy_params) :: energy_params_grow_sim, energy_params_move_sim
  
  logical :: sim_completed
  logical :: valid_tree
  logical :: input_specified, log_specified
  logical :: run_move
  logical :: file_exists

  integer :: run_flag
  integer :: i_arg, num_cl_args, num_threads
  integer :: i_node, i_rep
  integer :: unit_rep
  integer :: config_test(1:2), node_obst_test(1:2)
  integer, allocatable :: temp_ids(:)

  character(len=120) :: output_file, restart_file, output_label, output_dir, input_file, log_file
  character(len=120) :: format_string, temp_output_file, arg_specifier, arg



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  input_specified=.false.
  log_specified=.false.
  num_threads=max_num_threads

  num_cl_args=command_argument_count()

  if (num_cl_args .gt. 0) then
     do i_arg=1,num_cl_args
        call getarg(i_arg,arg)
        arg_specifier=trim(adjustl(arg(:index(arg,"=")-1)))
        arg=trim(adjustl(arg(index(arg,"=")+1:)))
        if (arg_specifier .eq. "input") then
           input_file=trim(arg)
           write(6,*)"input_file = ",trim(input_file)
           input_specified=.true.
        elseif (arg_specifier .eq. "log") then
           log_file=trim(arg)
           write(6,*)"log_file = ",trim(log_file)
           log_specified=.true.
        elseif (arg_specifier .eq. "num_threads") then
           read(arg,"(I3)")num_threads
           num_threads=min(num_threads,max_num_threads)
        end if
     end do
  end if

  !if the log unit is not stdout, open a log file at that unit number
  if (unit_log .ne. 6) then
     if (log_specified .eqv. .false.) then
        log_file="./run.log"
     end if
     open(UNIT=unit_log,FILE=trim(log_file),ACTION='write',STATUS='replace')
  end if

  if (input_specified .eqv. .false.) then
     !find input location
     input_file="../input/base_input.txt"
     open(UNIT=offset_units,FILE=trim(input_file),ACTION='read',STATUS='old')
     read(offset_units,"(A)")input_file
     input_file=trim(adjustl(input_file(index(input_file,"=")+1:)))
     close(offset_units)
  end if

  write(unit_log,*)'[ring_generate] BEGIN'

  !read detailed input
  call load_input_file(input_file,ring_init,obstacle_distribution_init,&
       deform_params_sim,grow_params_sim,&
       energy_params_move_sim,energy_params_grow_sim,&
       output_dir,output_label,run_flag)

  !write the system information to the log file
  write(unit_log,*)' [ring_generate] ------ PROGRAM DIRECTIVES -----'
  if (unit_log .ne. 6) then
     write(unit_log,*)' [ring_generate] log_file =',trim(log_file)
  end if
  write(unit_log,*)' [ring_generate] input_file =',trim(input_file)
  write(unit_log,*)' [ring_generate] run_flag =',run_flag
  write(unit_log,*)' [ring_generate] num_threads =',num_threads
  write(unit_log,*)' [ring_generate] output_dir =',trim(output_dir)
  write(unit_log,*)' [ring_generate] output_label =',trim(output_label)
  write(unit_log,*)' [ring_generate] ------ GROWTH PARAMETERS -----'
  write(unit_log,*)' [ring_generate] grow_flag =',grow_params_sim%grow_flag
  write(unit_log,*)' [ring_generate] growth_target =',grow_params_sim%growth_target
  write(unit_log,*)' [ring_generate] final_species =',trim(grow_params_sim%species_file)
  write(unit_log,*)' [ring_generate] N_final_species =',grow_params_sim%N_species_final
  write(unit_log,*)' [ring_generate] grow_energy_file =',trim(grow_params_sim%energy_file)
  write(unit_log,*)'  [ring_generate] V_nn_type =',energy_params_grow_sim%V_nn_type
  write(unit_log,*)'  [ring_generate] V_bend_type =',energy_params_grow_sim%V_bend_type
  write(unit_log,*)'  [ring_generate] V_twist_type =',energy_params_grow_sim%V_twist_type
  write(unit_log,*)'  [ring_generate] V_coord_type =',energy_params_grow_sim%V_coord_type
  write(unit_log,*)'  [ring_generate] V_obst_type =',energy_params_grow_sim%V_obst_type
  write(unit_log,*)'  [ring_generate] V_spectwist_type =',energy_params_grow_sim%V_spectwist_type
  write(unit_log,*)'  [ring_generate] V_pairwise_type =',energy_params_grow_sim%V_pairwise_type
  write(unit_log,*)' [ring_generate] max_growth =',grow_params_sim%max_growth
  write(unit_log,*)' [ring_generate] initial_growth =',grow_params_sim%initial_growth
  write(unit_log,*)' [ring_generate] rand_init =',grow_params_sim%rand_init
  write(unit_log,"(A,ES20.10)")" [ring_generate] alpha_1 = ",grow_params_sim%alpha(1)
  write(unit_log,"(A,ES20.10)")" [ring_generate] eta_1 = ",grow_params_sim%eta(1)
  write(unit_log,*)" [ring_generate] min_iter_move = ",grow_params_sim%min_iter_move
  write(unit_log,"(A,ES20.10)")" [ring_generate] alpha_2 = ",grow_params_sim%alpha(2)
  write(unit_log,"(A,ES20.10)")" [ring_generate] eta_2 = ",grow_params_sim%eta(2)
  write(unit_log,*)" [ring_generate] max_wait = ",grow_params_sim%max_wait
  write(unit_log,*)' [ring_generate] max_iter =',grow_params_sim%max_iter
  write(unit_log,*)' [ring_generate] err_check_freq =',grow_params_sim%err_check_freq
  write(unit_log,*)' [ring_generate] ------ DEFORMATION PARAMETERS -----'
  write(unit_log,*)' [ring_generate] min_rep =',deform_params_sim%min_rep
  write(unit_log,*)' [ring_generate] max_rep =',deform_params_sim%max_rep
  write(unit_log,*)' [ring_generate] N_reps =',deform_params_sim%N_reps
  write(unit_log,*)' [ring_generate] ergodic_flag =',deform_params_sim%ergodic_flag
  write(unit_log,*)' [ring_generate] move_energy_file =',trim(deform_params_sim%energy_file)
  write(unit_log,*)'  [ring_generate] V_nn_type =',energy_params_move_sim%V_nn_type
  write(unit_log,*)'  [ring_generate] V_bend_type =',energy_params_move_sim%V_bend_type
  write(unit_log,*)'  [ring_generate] V_twist_type =',energy_params_move_sim%V_twist_type
  write(unit_log,*)'  [ring_generate] V_coord_type =',energy_params_move_sim%V_coord_type
  write(unit_log,*)'  [ring_generate] V_obst_type =',energy_params_move_sim%V_obst_type
  write(unit_log,*)'  [ring_generate] V_spectwist_type =',energy_params_grow_sim%V_spectwist_type
  write(unit_log,*)'  [ring_generate] V_pairwise_type =',energy_params_move_sim%V_pairwise_type
  write(unit_log,*)' [ring_generate] num_iter =',deform_params_sim%num_iter
  write(unit_log,*)' [ring_generate] ref_init =',deform_params_sim%ref_init
  write(unit_log,*)' [ring_generate] ref_wait =',deform_params_sim%ref_wait
  write(unit_log,*)' [ring_generate] ref_freq =',deform_params_sim%ref_freq
  write(unit_log,*)' [ring_generate] ref_mult =',deform_params_sim%ref_mult
  write(unit_log,*)' [ring_generate] kink_mult =',deform_params_sim%kink_mult
  write(unit_log,*)' [ring_generate] kink_repeat =',deform_params_sim%kink_repeat
  write(unit_log,*)' [ring_generate] obst_freq =',deform_params_sim%obst_freq
  write(unit_log,*)' [ring_generate] obst_mult =',deform_params_sim%obst_mult
  write(unit_log,*)' [ring_generate] max_iter =',deform_params_sim%max_iter
  write(unit_log,*)' [ring_generate] err_check_freq =',deform_params_sim%err_check_freq
  write(unit_log,*)' [ring_generate] ------ SYSTEM DESCRIPTION -----'
  write(unit_log,*)' [ring_generate] dim =',ring_init%dim
  write(unit_log,*)' [ring_generate] N_nodes =',ring_init%N_nodes
  write(unit_log,*)' [ring_generate] fixed =',ring_init%fixed
  write(unit_log,*)' [ring_generate] N_fixed =',ring_init%N_fixed
  write(unit_log,*)' [ring_generate] N_species =',ring_init%N_species
  write(unit_log,*)' [ring_generate] bound_type =',ring_init%bound_type
  write(unit_log,*)' [ring_generate] N_bounds =',ring_init%N_bounds
  write(unit_log,*)' [ring_generate] obst_type =',obstacle_distribution_init%obst_type
  write(unit_log,*)' [ring_generate] N_obst =',obstacle_distribution_init%N
  write(unit_log,*)' [ring_generate] N_placed =',obstacle_distribution_init%N_placed
  write(unit_log,*)' [ring_generate] ring length =',ring_length(ring_init)

  write(unit_log,*)' [ring_generate] ------------------------------'
  write(unit_log,*)' [ring_generate] -------- MAIN ROUTINE --------'
  write(unit_log,*)' [ring_generate] ------------------------------'
  write(unit_log,*)' [ring_generate] testing initial configuration'


  call convert_obstacle_distribution(obstacle_distribution_init)  
  
  !test the initial configuration before proceeding
  !prematurely kill the program if any of the tests are failed

  !test the ring length to guarantee connectedness
  if (ring_length(ring_init) .eq. ring_init%N_nodes) then
     write(unit_log,*)'   [ring_generate] ring is closed correctly'
  else
     write(unit_log,*)'   [ring_generate] ERROR - invalid configuration - ring is not closed'
     stop
  end if

  !test for self-intersections
  config_test=valid_config(ring_init)
  
  if (sum(config_test) .eq. -2) then
     write(unit_log,*)'   [ring_generate] no self-intersections detected'
  else
     write(unit_log,*)'   [ring_generate] ERROR - invalid configuration - self-intersection was found'
     write(unit_log,"(a,i6,a,i6)")'    [ring_generate] intersection at nodes',&
          config_test(1),' and ',config_test(2)
     stop
  end if

  if (grow_params_sim%rand_init .eq. 0) then

     !test for obstacle intersections
     node_obst_test=obstacle_test(ring_init,obstacle_distribution_init)

     if (sum(node_obst_test) .eq. -2) then
        write(unit_log,*)'   [ring_generate] no obstacle intersections detected - slow'
     else
        write(unit_log,*)'   [ring_generate] ERROR - invalid configuration - obstacle intersection was found'
        write(unit_log,"(a,i6,a,i6)")'    [ring_generate] intersection of obstacle',&
             node_obst_test(2),' by node ',node_obst_test(1)
        stop
     end if

  end if
  

  write(unit_log,*)'  [ring_generate] TESTS PASSED - proceed with ring generation'



  !if (obstacle_distribution_init%obst_type .ne. 0)  then
  if (obstacle_distribution_init%N .ne. 0)  then
     write(unit_log,*)' [ring_generate] begin building obstacle tree'
     call octree_init(obstacle_distribution_init%obst_tree,&
          init_bbox(obstacle_distribution_init%x))

     allocate(temp_ids(1:obstacle_distribution_init%N))

     do i_node=1,obstacle_distribution_init%N
        temp_ids(i_node)=i_node
     end do

     call octree_build(obstacle_distribution_init%N,&
          temp_ids,&
          obstacle_distribution_init%x,&
          obstacle_distribution_init%placed,&
          obstacle_distribution_init%obst_tree%root_node)


     write(unit_log,*)' [ring_generate] Test the tree'
     
     call validate_tree(valid_tree,obstacle_distribution_init)
     
     if (valid_tree .eqv. .true.) then
        write(unit_log,*)'  [ring_generate] VALID tree'
     else
        write(unit_log,*)'  [ring_generate] INVALID tree'
        stop
     end if
     
     deallocate(temp_ids)
     write(unit_log,*)' [ring_generate] end building obstacle tree'

     ! call print_tree(obstacle_distribution_init%obst_tree%root_node,&
     !      obstacle_distribution_init%N,&
     !      obstacle_distribution_init%x)

     if (grow_params_sim%rand_init .eq. 0) then
     
        node_obst_test=obstacle_test_fast(ring_init,obstacle_distribution_init)

        if (sum(node_obst_test) .eq. -2) then
           write(unit_log,*)'   [ring_generate] no obstacle intersections detected in octree'
        else
           write(unit_log,*)'   [ring_generate] ERROR - invalid configuration - obstacle intersection found in octree'
           write(unit_log,"(a,i6,a,i6)")'    [ring_generate] intersection of obstacle',&
                node_obst_test(2),' by node ',node_obst_test(1)
           stop
        end if

     end if

  end if

  run_move=.true.
  
  !if the run flag was set to 1, then proceed with serial execution
  if (run_flag .eq. 1) then
     
     write(unit_log,*)' [ring_generate] ----- SERIAL EXECUTION: BEGIN -----'

     !set the number of OpenBLAS threads for the optimized BLAS library to match the threads
     !specified in the params module
     ! call openblas_set_num_threads(num_threads) !uncomment if using OpenBLAS

     !if the number of replicates is greater than 1, then iterate over the replicate numbers
     if (deform_params_sim%N_reps .gt. 1) then

        format_string="(a,i5.5)"
        
        do i_rep=deform_params_sim%min_rep,deform_params_sim%max_rep

           ring_rep=new_ring(ring_init%dim,ring_init%N_nodes,ring_init%N_fixed,ring_init%fixed,&
                ring_init%N_species,ring_init%bound_type,ring_init%N_bounds)
           ring_rep%x=ring_init%x
           ring_rep%fixed_nodes=ring_init%fixed_nodes
           ring_rep%species=ring_init%species
           ring_rep%bounds=ring_init%bounds

           call copy_obst_dist(obstacle_distribution_rep,obstacle_distribution_init)

           if (grow_params_sim%rand_init .eq. 1) then

              call rand_init_ring(ring_rep,obstacle_distribution_rep)

           end if

           
           !write(unit_log,*)' [ring_generate] --- Replicate',i_rep,' ---'
           
           !write(unit_log,*)' [ring_generate] ring before moving, length=',ring_length(ring_rep)

           unit_rep=offset_units+(i_rep-deform_params_sim%min_rep)

           output_file=trim(output_dir)//trim(output_label)//"_rep"
           write(output_file,format_string)trim(output_file),i_rep
           restart_file=trim(output_file)//"_restart.inp"


           file_exists=.false.
           inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
           if (file_exists .eqv. .true.) then
              open(UNIT=unit_rep,FILE=trim(output_file)//"_energy.dat",STATUS="old")
              close(UNIT=unit_rep,STATUS="delete")
           end if

           file_exists=.false.
           inquire(FILE=trim(output_file)//"_energy_fin.dat",EXIST=file_exists)
           if (file_exists .eqv. .true.) then
              open(UNIT=unit_rep,FILE=trim(output_file)//"_energy_fin.dat",STATUS="old")
              close(UNIT=unit_rep,STATUS="delete")
           end if
           
           if ((grow_params_sim%grow_flag .ne. 0) .and.&
              (grow_params_sim%growth_target .gt. ring_rep%N_nodes)) then

              run_move=.false.
              
              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   (obstacle_distribution_rep%N_placed .gt. 0)) then

                 call grow_ring_obstacles(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_grow_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,i_rep)

              else

                 call grow_ring(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_grow_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,i_rep)

              end if

              !update the species
              if (sim_completed .eqv. .true.) then
                 ring_rep%N_species=grow_params_sim%N_species_final
                 ring_rep%species=grow_params_sim%species_final
              end if

              if (grow_params_sim%grow_flag .ge. 3) then
                 run_move=.true.

                 temp_output_file=trim(output_file)//"_grow.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_x(ring_rep,unit_rep)

                 close(unit_rep)

                 temp_output_file=trim(output_file)//"_species_grow.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_species_x(ring_rep,unit_rep)

                 close(unit_rep)
                 
              end if

           else
              run_move=.true.
           end if

           if (run_move .eqv. .true.) then

              !conditional for MH sampling if there are potentials
              if ((energy_params_move_sim%V_nn_type .ne. 0) .or.&
                   (energy_params_move_sim%V_bend_type .ne. 0) .or.&
                   (energy_params_move_sim%V_twist_type .ne. 0) .or.&
                   (energy_params_move_sim%V_coord_type .ne. 0) .or.&
                   (energy_params_move_sim%V_obst_type .ne. 0) .or.&
                   (energy_params_move_sim%V_pairwise_type .ne. 0)) then

                 if ((energy_params_move_sim%V_pairwise_type .ne. 0) .and.&
                      (ring_rep%N_species .gt. 1)) then
                    call minimize_ring_V_pairwise(ring_rep,energy_params_move_sim)
                 end if
                 
                 !conditional for testing obstacles
                 if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                      ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                      (deform_params_sim%obst_freq .gt. 0))) then

                    call move_ring_obstacles_MH(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         energy_params_move_sim,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 else

                    call move_ring_MH(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         energy_params_move_sim,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 end if

              else

                 !conditional for testing obstacles
                 if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                      ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                      (deform_params_sim%obst_freq .gt. 0))) then

                    call move_ring_obstacles(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 else

                    call move_ring(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 end if

              end if

           end if

           !write(unit_log,*)' [ring_generate] ring after moving, length=',ring_length(ring_rep)

           if (sim_completed .eqv. .true.) then
           
              open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

              call write_restart_file(ring_rep,obstacle_distribution_rep,&
                   deform_params_sim,grow_params_sim,unit_rep,.true.)

              close(unit_rep)

              temp_output_file=trim(output_file)//".dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_x(ring_rep,unit_rep)

              close(unit_rep)

              temp_output_file=trim(output_file)//"_species.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_species_x(ring_rep,unit_rep)

              close(unit_rep)

              if (obstacle_distribution_rep%obst_base%N_placed .gt. 0) then

                 temp_output_file=trim(output_file)//"_obst.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_obst_placed_x(obstacle_distribution_rep,unit_rep)

                 close(unit_rep)

              end if

              if (track_energy .eqv. .true.) then
                 file_exists=.false.
                 inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
                 if (file_exists .eqv. .true.) then
                    call rename(trim(output_file)//"_energy.dat",trim(output_file)//"_energy_fin.dat")
                 end if
              end if

           end if

           if (obstacle_distribution_rep%obst_type .ne. 0) then
              call octree_final(obstacle_distribution_rep%obst_tree)
           end if

        end do

     else

        ring_rep=new_ring(ring_init%dim,ring_init%N_nodes,ring_init%N_fixed,ring_init%fixed,&
             ring_init%N_species,ring_init%bound_type,ring_init%N_bounds)
        ring_rep%x=ring_init%x
        ring_rep%fixed_nodes=ring_init%fixed_nodes
        ring_rep%species=ring_init%species
        ring_rep%bounds=ring_init%bounds

        call copy_obst_dist(obstacle_distribution_rep,obstacle_distribution_init)

        if (grow_params_sim%rand_init .eq. 1) then

           call rand_init_ring(ring_rep,obstacle_distribution_rep)

        end if

        !write(unit_log,*)' [ring_generate] ring before moving, length=',ring_length(ring_rep)

        unit_rep=offset_units+deform_params_sim%min_rep

        output_file=trim(output_dir)//trim(output_label)
        restart_file=trim(output_file)//"_restart.inp"

        file_exists=.false.
        inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
        if (file_exists .eqv. .true.) then
           open(UNIT=unit_rep,FILE=trim(output_file)//"_energy.dat",STATUS="old")
           close(UNIT=unit_rep,STATUS="delete")
        end if

        file_exists=.false.
        inquire(FILE=trim(output_file)//"_energy_fin.dat",EXIST=file_exists)
        if (file_exists .eqv. .true.) then
           open(UNIT=unit_rep,FILE=trim(output_file)//"_energy_fin.dat",STATUS="old")
           close(UNIT=unit_rep,STATUS="delete")
        end if

        if ((grow_params_sim%grow_flag .ne. 0) .and.&
              (grow_params_sim%growth_target .gt. ring_rep%N_nodes)) then

           run_move=.false.
           
           !conditional for testing obstacles
           if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                (obstacle_distribution_rep%N_placed .gt. 0)) then

              call grow_ring_obstacles(sim_completed,&
                   ring_rep,obstacle_distribution_rep,&
                   energy_params_grow_sim,&
                   deform_params_sim,grow_params_sim,&
                   restart_file,deform_params_sim%min_rep)

           else

              call grow_ring(sim_completed,&
                   ring_rep,obstacle_distribution_rep,&
                   energy_params_grow_sim,&
                   deform_params_sim,grow_params_sim,&
                   restart_file,deform_params_sim%min_rep)

           end if

           !update the species
           if (sim_completed .eqv. .true.) then
              ring_rep%N_species=grow_params_sim%N_species_final
              ring_rep%species=grow_params_sim%species_final
           end if

           if (grow_params_sim%grow_flag .ge. 3) then
              run_move=.true.

              temp_output_file=trim(output_file)//"_grow.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_x(ring_rep,unit_rep)

              close(unit_rep)

              temp_output_file=trim(output_file)//"_species_grow.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_species_x(ring_rep,unit_rep)

              close(unit_rep)
              
           end if

        else
           run_move=.true.
        end if

        if (run_move .eqv. .true.) then

           !conditional for MH sampling if there are potentials
           if ((energy_params_move_sim%V_nn_type .ne. 0) .or.&
                (energy_params_move_sim%V_bend_type .ne. 0) .or.&
                (energy_params_move_sim%V_twist_type .ne. 0) .or.&
                (energy_params_move_sim%V_coord_type .ne. 0) .or.&
                (energy_params_move_sim%V_obst_type .ne. 0) .or.&
                (energy_params_move_sim%V_pairwise_type .ne. 0)) then

              if ((energy_params_move_sim%V_pairwise_type .ne. 0) .and.&
                   (ring_rep%N_species .gt. 1)) then
                 call minimize_ring_V_pairwise(ring_rep,energy_params_move_sim)
              end if

              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                   (deform_params_sim%obst_freq .gt. 0))) then

                 call move_ring_obstacles_MH(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_move_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              else

                 call move_ring_MH(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_move_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              end if

           else

              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                   (deform_params_sim%obst_freq .gt. 0))) then

                 call move_ring_obstacles(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              else

                 call move_ring(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              end if

           end if

        end if


        !write(unit_log,*)' [ring_generate] ring after moving, length=',ring_length(ring_rep)

        if (sim_completed .eqv. .true.) then
        
           open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

           call write_restart_file(ring_rep,obstacle_distribution_rep,&
                deform_params_sim,grow_params_sim,unit_rep,.true.)

           close(unit_rep)

           temp_output_file=trim(output_file)//".dat"

           open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

           call write_x(ring_rep,unit_rep)

           close(unit_rep)

           temp_output_file=trim(output_file)//"_species.dat"

           open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

           call write_species_x(ring_rep,unit_rep)

           close(unit_rep)

           if (obstacle_distribution_rep%obst_base%N_placed .gt. 0) then

              temp_output_file=trim(output_file)//"_obst.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_obst_placed_x(obstacle_distribution_rep,unit_rep)

              close(unit_rep)

           end if

           if (track_energy .eqv. .true.) then
              file_exists=.false.
              inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
              if (file_exists .eqv. .true.) then
                 call rename(trim(output_file)//"_energy.dat",trim(output_file)//"_energy_fin.dat")
              end if
           end if

        end if

        if (obstacle_distribution_rep%obst_type .ne. 0) then
           call octree_final(obstacle_distribution_rep%obst_tree)
        end if

     end if
     
     write(unit_log,*)' [ring_generate] ----- SERIAL EXECUTION: END -----'


  !if the run flag was set to 2, then proceed with parallel execution
  elseif (run_flag .eq. 2) then
     
     write(unit_log,*)' [ring_generate] ----- PARALLEL EXECUTION: BEGIN -----'

     
     if (deform_params_sim%N_reps .gt. 1) then

        format_string="(a,i5.5)"

        !set the number of OpenMP threads to match the threads
        !specified in the params module
        call omp_set_num_threads(num_threads)

        !set the number of OpenBLAS threads for the optimized BLAS library to 1
        !call openblas_set_num_threads(1)
        
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(deform_params_sim,grow_params_sim,&
        !$OMP energy_params_move_sim,energy_params_grow_sim,&
        !$OMP ring_init,obstacle_distribution_init,&
        !$OMP output_dir,output_label,format_string)
        
        do i_rep=deform_params_sim%min_rep,deform_params_sim%max_rep

           ring_rep=new_ring(ring_init%dim,ring_init%N_nodes,ring_init%N_fixed,ring_init%fixed,&
                ring_init%N_species,ring_init%bound_type,ring_init%N_bounds)
           ring_rep%x=ring_init%x
           ring_rep%fixed_nodes=ring_init%fixed_nodes
           ring_rep%species=ring_init%species
           ring_rep%bounds=ring_init%bounds
           
           call copy_obst_dist(obstacle_distribution_rep,obstacle_distribution_init)

           if (grow_params_sim%rand_init .eq. 1) then

              call rand_init_ring(ring_rep,obstacle_distribution_rep)

           end if
           
           !write(unit_log,*)' [ring_generate] --- Replicate',i_rep,' ---'
           !write(unit_log,*)' [ring_generate] ring before moving, length=',ring_length(ring_rep)

           unit_rep=offset_units+(i_rep-deform_params_sim%min_rep)
           
           output_file=trim(output_dir)//trim(output_label)//"_rep"
           write(output_file,format_string)trim(output_file),i_rep
           restart_file=trim(output_file)//"_restart.inp"

           file_exists=.false.
           inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
           if (file_exists .eqv. .true.) then
              open(UNIT=unit_rep,FILE=trim(output_file)//"_energy.dat",STATUS="old")
              close(UNIT=unit_rep,STATUS="delete")
           end if

           file_exists=.false.
           inquire(FILE=trim(output_file)//"_energy_fin.dat",EXIST=file_exists)
           if (file_exists .eqv. .true.) then
              open(UNIT=unit_rep,FILE=trim(output_file)//"_energy_fin.dat",STATUS="old")
              close(UNIT=unit_rep,STATUS="delete")
           end if
           
           if ((grow_params_sim%grow_flag .ne. 0) .and.&
              (grow_params_sim%growth_target .gt. ring_rep%N_nodes)) then

              run_move=.false.
              
              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   (obstacle_distribution_rep%N_placed .gt. 0)) then

                 call grow_ring_obstacles(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_grow_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,i_rep)

              else

                 call grow_ring(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_grow_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,i_rep)

              end if

              !update the species
              if (sim_completed .eqv. .true.) then
                 ring_rep%N_species=grow_params_sim%N_species_final
                 ring_rep%species=grow_params_sim%species_final
              end if
              
              if (grow_params_sim%grow_flag .ge. 3) then
                 run_move=.true.

                 temp_output_file=trim(output_file)//"_grow.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_x(ring_rep,unit_rep)

                 close(unit_rep)

                 temp_output_file=trim(output_file)//"_species_grow.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_species_x(ring_rep,unit_rep)

                 close(unit_rep)
                 
              end if
              
           else
              run_move=.true.
           end if

           if (run_move .eqv. .true.) then

              !conditional for MH sampling if there are potentials
              if ((energy_params_move_sim%V_nn_type .ne. 0) .or.&
                   (energy_params_move_sim%V_bend_type .ne. 0) .or.&
                   (energy_params_move_sim%V_twist_type .ne. 0) .or.&
                   (energy_params_move_sim%V_coord_type .ne. 0) .or.&
                   (energy_params_move_sim%V_obst_type .ne. 0) .or.&
                   (energy_params_move_sim%V_pairwise_type .ne. 0)) then


                 if ((energy_params_move_sim%V_pairwise_type .ne. 0) .and.&
                   (ring_rep%N_species .gt. 1)) then
                    call minimize_ring_V_pairwise(ring_rep,energy_params_move_sim)
                 end if

                 !conditional for testing obstacles
                 if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                      ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                      (deform_params_sim%obst_freq .gt. 0))) then

                    call move_ring_obstacles_MH(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         energy_params_move_sim,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 else

                    call move_ring_MH(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         energy_params_move_sim,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 end if

              else

                 !conditional for testing obstacles
                 if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                      ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                      (deform_params_sim%obst_freq .gt. 0))) then

                    call move_ring_obstacles(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 else

                    call move_ring(sim_completed,&
                         ring_rep,obstacle_distribution_rep,&
                         deform_params_sim,grow_params_sim,&
                         restart_file,i_rep)

                 end if

              end if

           end if
           
           !write(unit_log,*)' [ring_generate] ring after moving, length=',ring_length(ring_rep)

           if (sim_completed .eqv. .true.) then
           
              open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

              call write_restart_file(ring_rep,obstacle_distribution_rep,&
                   deform_params_sim,grow_params_sim,unit_rep,.true.)

              close(unit_rep)

              temp_output_file=trim(output_file)//".dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_x(ring_rep,unit_rep)

              close(unit_rep)

              temp_output_file=trim(output_file)//"_species.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_species_x(ring_rep,unit_rep)

              close(unit_rep)

              if (obstacle_distribution_rep%obst_base%N_placed .gt. 0) then

                 temp_output_file=trim(output_file)//"_obst.dat"

                 open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

                 call write_obst_placed_x(obstacle_distribution_rep,unit_rep)

                 close(unit_rep)

              end if

              if (track_energy .eqv. .true.) then
                 file_exists=.false.
                 inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
                 if (file_exists .eqv. .true.) then
                    call rename(trim(output_file)//"_energy.dat",trim(output_file)//"_energy_fin.dat")
                 end if
              end if

           end if

           if (obstacle_distribution_rep%obst_type .ne. 0) then
              call octree_final(obstacle_distribution_rep%obst_tree)
           end if

        end do

     else

        !set the number of OpenBLAS threads for the optimized BLAS library to match the threads
        !specified in the params module
        ! call openblas_set_num_threads(num_threads) !uncomment if using OpenBLAS

        ring_rep=new_ring(ring_init%dim,ring_init%N_nodes,ring_init%N_fixed,ring_init%fixed,&
             ring_init%N_species,ring_init%bound_type,ring_init%N_bounds)
        ring_rep%x=ring_init%x
        ring_rep%fixed_nodes=ring_init%fixed_nodes
        ring_rep%species=ring_init%species
        ring_rep%bounds=ring_init%bounds

        call copy_obst_dist(obstacle_distribution_rep,obstacle_distribution_init)

        if (grow_params_sim%rand_init .eq. 1) then

           call rand_init_ring(ring_rep,obstacle_distribution_rep)

        end if

        
        !write(unit_log,*)' [ring_generate] ring before moving, length=',ring_length(ring_rep)

        unit_rep=offset_units+deform_params_sim%min_rep

        output_file=trim(output_dir)//trim(output_label)
        restart_file=trim(output_file)//"_restart.inp"

        file_exists=.false.
        inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
        if (file_exists .eqv. .true.) then
           open(UNIT=unit_rep,FILE=trim(output_file)//"_energy.dat",STATUS="old")
           close(UNIT=unit_rep,STATUS="delete")
        end if

        file_exists=.false.
        inquire(FILE=trim(output_file)//"_energy_fin.dat",EXIST=file_exists)
        if (file_exists .eqv. .true.) then
           open(UNIT=unit_rep,FILE=trim(output_file)//"_energy_fin.dat",STATUS="old")
           close(UNIT=unit_rep,STATUS="delete")
        end if
        
        if ((grow_params_sim%grow_flag .ne. 0) .and.&
              (grow_params_sim%growth_target .gt. ring_rep%N_nodes)) then

           run_move=.false.
           
           !conditional for testing obstacles
           if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                (obstacle_distribution_rep%N_placed .gt. 0)) then

              call grow_ring_obstacles(sim_completed,&
                   ring_rep,obstacle_distribution_rep,&
                   energy_params_grow_sim,&
                   deform_params_sim,grow_params_sim,&
                   restart_file,deform_params_sim%min_rep)

           else

              call grow_ring(sim_completed,&
                   ring_rep,obstacle_distribution_rep,&
                   energy_params_grow_sim,&
                   deform_params_sim,grow_params_sim,&
                   restart_file,deform_params_sim%min_rep)

           end if

           !update the species
           if (sim_completed .eqv. .true.) then
              ring_rep%N_species=grow_params_sim%N_species_final
              ring_rep%species=grow_params_sim%species_final
           end if

           if (grow_params_sim%grow_flag .ge. 3) then
              run_move=.true.

              temp_output_file=trim(output_file)//"_grow.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_x(ring_rep,unit_rep)

              close(unit_rep)

              temp_output_file=trim(output_file)//"_species_grow.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_species_x(ring_rep,unit_rep)

              close(unit_rep)
                 
           end if

        else
           run_move=.true.
        end if

        if (run_move .eqv. .true.) then

           !conditional for MH sampling if there are potentials
           if ((energy_params_move_sim%V_nn_type .ne. 0) .or.&
                (energy_params_move_sim%V_bend_type .ne. 0) .or.&
                (energy_params_move_sim%V_twist_type .ne. 0) .or.&
                (energy_params_move_sim%V_coord_type .ne. 0) .or.&
                (energy_params_move_sim%V_obst_type .ne. 0) .or.&
                (energy_params_move_sim%V_pairwise_type .ne. 0)) then

              if ((energy_params_move_sim%V_pairwise_type .ne. 0) .and.&
                   (ring_rep%N_species .gt. 1)) then
                 call minimize_ring_V_pairwise(ring_rep,energy_params_move_sim)
              end if

              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                   (deform_params_sim%obst_freq .gt. 0))) then

                 call move_ring_obstacles_MH(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_move_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              else

                 call move_ring_MH(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      energy_params_move_sim,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              end if

           else

              !conditional for testing obstacles
              if ((obstacle_distribution_rep%obst_type .ne. 0) .and.&
                   ((obstacle_distribution_rep%N_placed .gt. 0) .or.&
                   (deform_params_sim%obst_freq .gt. 0))) then

                 call move_ring_obstacles(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              else

                 call move_ring(sim_completed,&
                      ring_rep,obstacle_distribution_rep,&
                      deform_params_sim,grow_params_sim,&
                      restart_file,deform_params_sim%min_rep)

              end if

           end if

        end if

        !write(unit_log,*)' [ring_generate] ring after moving, length=',ring_length(ring_rep)

        if (sim_completed .eqv. .true.) then
        
           open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

           call write_restart_file(ring_rep,obstacle_distribution_rep,&
                deform_params_sim,grow_params_sim,unit_rep,.true.)

           close(unit_rep)

           temp_output_file=trim(output_file)//".dat"

           open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

           call write_x(ring_rep,unit_rep)

           close(unit_rep)

           temp_output_file=trim(output_file)//"_species.dat"

           open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

           call write_species_x(ring_rep,unit_rep)

           close(unit_rep)

           if (obstacle_distribution_rep%obst_base%N_placed .gt. 0) then

              temp_output_file=trim(output_file)//"_obst.dat"

              open(UNIT=unit_rep,FILE=trim(temp_output_file),ACTION='write',STATUS='replace')

              call write_obst_placed_x(obstacle_distribution_rep,unit_rep)

              close(unit_rep)

           end if

           if (track_energy .eqv. .true.) then
              file_exists=.false.
              inquire(FILE=trim(output_file)//"_energy.dat",EXIST=file_exists)
              if (file_exists .eqv. .true.) then
                 call rename(trim(output_file)//"_energy.dat",trim(output_file)//"_energy_fin.dat")
              end if
           end if

        end if

        if (obstacle_distribution_rep%obst_type .ne. 0) then
           call octree_final(obstacle_distribution_rep%obst_tree)
        end if

     end if
     
     write(unit_log,*)' [ring_generate] ----- PARALLEL EXECUTION: END -----'
     
  end if

  if (obstacle_distribution_init%obst_type .eq. 1) then
     write(unit_log,*)' [ring_generate] begin cleaning obstacle tree'
     call octree_final(obstacle_distribution_init%obst_tree)
     write(unit_log,*)' [ring_generate] end cleaning obstacle tree'
  end if
  
  
  write(unit_log,*)'[ring_generate] END'

  !if not writing to stdout, then close the log file
  if (unit_log .ne. 6) then
     close(unit_log)
  end if

end program ring_generate
