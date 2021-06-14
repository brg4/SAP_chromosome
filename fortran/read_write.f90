module read_write

  !module used to handle reading the input file and
  !writing the restart and output files
  
  use params
  use ring_objects
  use energy_terms

  implicit none

contains

  !routine to load the input file
  subroutine load_input_file(input_file,&
       out_ring,out_obstacle_distribution,&
       out_deform_params,out_grow_params,&
       out_energy_params_move,out_energy_params_grow,&
       output_dir,output_label,run_flag)

    implicit none

    character(len=120), intent(in) :: input_file
    type(ring), intent(out) :: out_ring
    type(grow_params), intent(out) :: out_grow_params
    type(deform_params), intent(out) :: out_deform_params
    type(obstacle_distribution), intent(out) :: out_obstacle_distribution
    type(energy_params), intent(out) :: out_energy_params_grow, out_energy_params_move
    integer, intent(out) :: run_flag
    character(len=120), intent(out) :: output_dir, output_label

    logical :: fixed
    
    integer :: unit_input
    
    integer :: grow_flag, growth_target, initial_growth, max_growth, rand_init, ref_node
    integer :: max_iter_grow, err_check_freq_grow
    integer :: min_iter_move, max_wait
    real(rp) :: alpha(1:2), eta(1:2)
    
    integer :: ergodic_flag, num_iter, max_iter, err_check_freq
    integer :: ref_init, ref_wait, ref_freq, ref_mult, kink_mult, kink_repeat
    integer :: obst_freq, obst_mult
    integer :: min_rep, max_rep
    
    integer :: dim, bound_type, obst_type
    integer :: N_nodes, N_fixed, N_species, N_bounds, N_obst, N_placed
    
    integer :: i_node, i_bound, i_obst, temp_id
    character(len=120) :: temp_line, move_energy_file, grow_energy_file, final_species_file

    integer :: V_nn_type_temp, V_bend_type_temp, V_twist_type_temp
    real(rp) :: nn_param_temp, bend_param_temp, twist_param_temp

    unit_input=offset_units
    
    open(UNIT=unit_input,FILE=trim(input_file),ACTION='read',STATUS='old')
    !skip ##### PROGRAM DIRECTIVES #####
    read(unit_input,*)
    !read the run_flag
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I1)")run_flag
    !read the output_dir
    read(unit_input,"(A)")temp_line
    output_dir=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    !print *,output_dir
    !read the output_label
    read(unit_input,"(A)")temp_line
    output_label=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    !print *,output_label

    !skip ##### GROWTH PARAMETERS #####
    read(unit_input,*)
    !read the grow_flag
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I3)")grow_flag
    !read the grow_flag
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I9)")growth_target
    !read the final_species_file
    read(unit_input,"(A)")temp_line
    final_species_file=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    !read the grow_energy_file
    read(unit_input,"(A)")temp_line
    grow_energy_file=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    !read the max_growth
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I9)")max_growth
    !read the initial_growth
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I9)")initial_growth
    !read the rand_init
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I3)")rand_init
    !read the ref_node
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")ref_node
    !read the first alpha
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,*)alpha(1)
    !read the first eta
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,*)eta(1)
    !read the minimum number of move iterations
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")min_iter_move
    !read the second alpha
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,*)alpha(2)
    !read the second eta
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,*)eta(2)
    !read the maximum wait
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")max_wait
    !read the growth max_iter
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")max_iter_grow
    !read the growth err_check_freq
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")err_check_freq_grow


    
    !skip ##### DEFORMATION PARAMETERS #####
    read(unit_input,*)
    !read the minimum replicate number
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I3)")min_rep
    !read the maximum replicate number
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I3)")max_rep
    !read the ergodicity flag
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I3)")ergodic_flag
    !read the move_energy_file
    read(unit_input,"(A)")temp_line
    move_energy_file=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    !read the number of iterations
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")num_iter
    !read the initial number of reflections
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")ref_init
    !read the number of iterations to wait before allowing reflections
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")ref_wait
    !read the reflection frequency
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")ref_freq
    !read the ref multiplicity
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")ref_mult
    !read the kink multiplicity
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")kink_mult
    !read the kink repeat number
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")kink_repeat
    !read the obstacle frequency
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")obst_freq
    !read the obstacle multiplicity
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")obst_mult
    !read the maximum number of sampling attempts
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")max_iter
    !read the error checking frequency
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")err_check_freq
    
    !skip ##### SYSTEM DESCRIPTION #####
    read(unit_input,*)
    !read the dimensions
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I1)")dim
    !read the number of nodes
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")N_nodes
    !read the number of fixed nodes
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")N_fixed
    !read the number of species
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")N_species
    !read the bound type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I2)")bound_type
    !read the number of bounds
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I2)")N_bounds
    !read the obstacle types
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I2)")obst_type
    !read the number of obstacles
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")N_obst
    !read the number of placed obstacles
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I7)")N_placed

    !initialize ring object and deform_params object
    if (N_fixed .gt. 0) then
       fixed=.true.
    else
       fixed=.false.
    end if

    
    out_ring=new_ring(dim,N_Nodes,N_fixed,fixed,N_species,bound_type,N_bounds)


    out_deform_params=new_deform_params(move_energy_file,&
         ergodic_flag,num_iter,&
         ref_init,ref_wait,ref_freq,ref_mult,&
         kink_mult,kink_repeat,&
         obst_freq,obst_mult,&
         max_iter,err_check_freq,&
         min_rep,max_rep)

    out_grow_params=new_grow_params(grow_energy_file,final_species_file,&
         grow_flag,growth_target,max_growth,initial_growth,rand_init,ref_node,&
         alpha,eta,&
         min_iter_move,max_wait,&
         err_check_freq_grow,max_iter_grow)
    
    out_obstacle_distribution=new_obstacle_distribution(dim,N_obst,obst_type)
    out_obstacle_distribution%N_placed=N_placed


    !skip ##### bounds #####
    read(unit_input,*)

    !read the bounds list
    if (out_ring%bound_type .ne. 0) then
       do i_bound=1,out_ring%N_bounds
          read(unit_input,*)out_ring%bounds(:,i_bound)
       end do
    end if

    !skip ##### placed obstacles #####
    read(unit_input,*)temp_line

    !read the placed obstacles list
    if (out_obstacle_distribution%N_placed .ne. 0) then
       do i_obst=1,out_obstacle_distribution%N_placed
          read(unit_input,"(I8)")out_obstacle_distribution%placed(i_obst)
       end do
    end if

    !skip ##### obstacle coordinates #####
    read(unit_input,*)

    !read the obstacle coordinate list
    if (out_obstacle_distribution%N .ne. 0) then
       do i_obst=1,out_obstacle_distribution%N
          read(unit_input,*)temp_id,out_obstacle_distribution%mapping(i_obst),&
               out_obstacle_distribution%x(:,i_obst)
       end do
    end if

    
    !skip ##### fixed nodes #####
    read(unit_input,*)

    !read the fixed node list
    if (out_ring%fixed .eqv. .true.) then
       do i_node=1,out_ring%N_fixed
          !read(unit_input,"(A)")temp_line
          !temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
          !read(temp_line,"(I4)")out_ring%fixed_nodes(i_node)
          read(unit_input,"(I8)")out_ring%fixed_nodes(i_node)
       end do
    end if
    
    !skip ##### node coordinates #####
    read(unit_input,*)

    !read the node coordinate list
    do i_node=1,out_ring%N_nodes
       !read(unit_input,*)out_ring%x(1:out_ring%dim,i_node)
       read(unit_input,*)temp_id,out_ring%species(i_node),out_ring%x(1:out_ring%dim,i_node)
    end do

    close(unit_input)

    !call a routine to open the species file and update the grow parameters
    call read_species_file(out_grow_params)
    
    call read_energy_file(out_grow_params%energy_file,&
         out_ring%dim,out_grow_params%N_species_final,&
         out_obstacle_distribution,out_energy_params_grow)

    call read_energy_file(out_deform_params%energy_file,&
         out_ring%dim,out_grow_params%N_species_final,&
         out_obstacle_distribution,out_energy_params_move)

    if ((out_grow_params%grow_flag .eq. 4) .or.&
         (out_grow_params%grow_flag .eq. 5)) then

       !nn
       V_nn_type_temp=out_energy_params_grow%V_nn_type
       if (out_energy_params_grow%V_nn_type .ne. 0) then
          nn_param_temp=out_energy_params_grow%nn_params(1)
       end if
       
       !bend
       V_bend_type_temp=out_energy_params_grow%V_bend_type
       if (out_energy_params_grow%V_bend_type .ne. 0) then
          bend_param_temp=out_energy_params_grow%bend_params(1)
       end if

       !twist
       V_twist_type_temp=out_energy_params_grow%V_twist_type
       if (out_energy_params_grow%V_twist_type .ne. 0) then
          twist_param_temp=out_energy_params_grow%twist_params(1)
       end if
       
       call read_energy_file(out_deform_params%energy_file,&
            out_ring%dim,out_grow_params%N_species_final,&
            out_obstacle_distribution,out_energy_params_grow)
       
       !nn
       out_energy_params_grow%V_nn_type=V_nn_type_temp
       out_energy_params_grow%nn_params(1)=nn_param_temp
       !bend
       out_energy_params_grow%V_bend_type=V_bend_type_temp
       out_energy_params_grow%bend_params(1)=bend_param_temp
       !twist
       out_energy_params_grow%V_twist_type=V_twist_type_temp
       out_energy_params_grow%twist_params(1)=twist_param_temp
       
    end if
    
  end subroutine load_input_file

  subroutine read_species_file(in_grow_params)

    implicit none

    type(grow_params), intent(inout) :: in_grow_params

    integer :: unit_input
    integer :: temp_id
    integer :: i_node

    character(len=120) :: temp_line

    unit_input=offset_units

    open(UNIT=unit_input,FILE=trim(in_grow_params%species_file),&
         ACTION='read',STATUS='old')

    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")in_grow_params%N_species_final

    if (allocated(in_grow_params%species_final)) deallocate(in_grow_params%species_final)
    allocate(in_grow_params%species_final(1:in_grow_params%growth_target))

    do i_node=1,in_grow_params%growth_target
       read(unit_input,*)temp_id,in_grow_params%species_final(i_node)
    end do
    
    close(unit_input)

  end subroutine read_species_file

  subroutine read_energy_file(energy_file,dim,N_species,&
       in_obstacle_distribution,out_energy_params)

    implicit none

    character(len=120), intent(in) :: energy_file
    
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    integer, intent(in) :: dim, N_species
    
    type(energy_params), intent(out) :: out_energy_params

    integer :: unit_input
    integer :: V_nn_type, V_bend_type, V_twist_type
    integer :: V_coord_type, V_obst_type, V_spectwist_type, V_pairwise_type
    integer :: N_V_coord, N_V_obst, N_spectwist, N_pairs
    integer :: N_nn_params, N_bend_params, N_twist_params
    integer :: N_coord_params, N_obst_params, N_spectwist_params, N_pair_params

    integer :: i_V, i_specie, N_species_temp, specie_temp_1, specie_temp_2
    integer :: temp_pair_type, temp_N_pair_params
    integer :: temp_spectwist_type, temp_N_spectwist_params
    
    real(rp) :: beta
    real(rp), allocatable :: temp_params(:)
    
    character(len=120) :: temp_line

    unit_input=offset_units
    
    open(UNIT=unit_input,FILE=trim(energy_file),ACTION='read',STATUS='old')
    
    !skip the first line
    read(unit_input,*)temp_line
    
    !read beta
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,*)beta
    !read V_nn_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_nn_type
    !read V_bend_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_bend_type
    !read V_twist_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_twist_type
    !read V_coord_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_coord_type
    !read V_obst_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_obst_type
    !read V_spectwist_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_spectwist_type
    !read V_pairwise_type
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I4)")V_pairwise_type

    !skip "Numbers of Interactions" Label
    read(unit_input,*)

    !read the number of V_coord
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")N_V_coord
    !read the number of V_obst
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I5)")N_V_obst
    !read the number of pairwise interactions
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")N_spectwist
    !read the number of pairwise interactions
    read(unit_input,"(A)")temp_line
    temp_line=trim(adjustl(temp_line(index(temp_line,"=")+1:)))
    read(temp_line,"(I6)")N_pairs

    if (V_nn_type .eq. 0) then
       N_nn_params=1
    elseif (V_nn_type .eq. 1) then !nn cost
       N_nn_params=1
    end if

    if (V_bend_type .eq. 0) then
       N_bend_params=1
    elseif (V_bend_type .eq. 1) then !bending cost
       N_bend_params=1
    end if

    if (V_twist_type .eq. 0) then
       N_twist_params=1
    elseif (V_twist_type .eq. 1) then !twisting cost
       N_twist_params=1
    end if
    
    if (V_coord_type .eq. 0) then
       N_coord_params=1
    elseif (V_coord_type .eq. 1) then !Harmonic potential
       N_coord_params=4
    end if

    if (V_obst_type .eq. 0) then
       N_obst_params=1
    elseif (V_obst_type .eq. 1) then !NN interactions
       N_obst_params=1
    end if

    if (V_spectwist_type .eq. 0) then
       N_spectwist_params=1
    elseif (V_spectwist_type .eq. 1) then !NN interactions
       N_spectwist_params=1
    end if

    if (V_pairwise_type .eq. 0) then
       N_pair_params=1
    elseif (V_pairwise_type .eq. 1) then !NN interactions
       N_pair_params=1
    elseif (V_pairwise_type .eq. 2) then !Harmonic interactions
       N_pair_params=1
    elseif (V_pairwise_type .eq. 3) then !Harmonic interactions with offset
       N_pair_params=2
    end if

    allocate(temp_params(1:N_pair_params))

    out_energy_params=new_energy_params(beta,&
         V_nn_type,N_nn_params,&
         V_bend_type,N_bend_params,&
         V_twist_type,N_twist_params,&
         V_coord_type,N_coord_params,N_V_coord,&
         V_obst_type,N_obst_params,N_V_obst,&
         V_spectwist_type,N_spectwist_params,N_spectwist,&
         V_pairwise_type,N_pair_params,N_species)

    !skip "V_nn" Label
    read(unit_input,*)temp_line

    if (out_energy_params%V_nn_type .ne. 0) then
       read(unit_input,*)out_energy_params%nn_params(1:out_energy_params%N_nn_params)
    end if

    !skip "V_bend" Label
    read(unit_input,*)temp_line
    
    if (out_energy_params%V_bend_type .ne. 0) then
       read(unit_input,*)out_energy_params%bend_params(1:out_energy_params%N_bend_params)
    end if

    !skip "V_twist" Label
    read(unit_input,*)temp_line

    if (out_energy_params%V_twist_type .ne. 0) then
       read(unit_input,*)out_energy_params%twist_params(1:out_energy_params%N_twist_params)
    end if
    
    !skip "V_coord" Label
    read(unit_input,*)temp_line

    if (out_energy_params%V_coord_type .ne. 0) then
       !for each V_coord, read the parameters
       !then read the number of affected species
       !then read the list of affected species
       do i_V=1,out_energy_params%N_V_coord
          read(unit_input,*)N_species_temp,&
               out_energy_params%coord_params(1:out_energy_params%N_coord_params,i_V)
          do i_specie=1,N_species_temp
             read(unit_input,*)specie_temp_1
             out_energy_params%species_coord(specie_temp_1,i_V)=1
          end do
       end do
    end if

    !skip "V_obst" Label
    read(unit_input,*)

    if (out_energy_params%V_obst_type .ne. 0) then
       !for each V_obst, read the parameters
       !then read the number of affected species
       !then read the list of affected species
       do i_V=1,out_energy_params%N_V_coord
          read(unit_input,*)N_species_temp,&
               out_energy_params%obst_params(1:out_energy_params%N_obst_params,i_V)
          do i_specie=1,N_species_temp
             read(unit_input,*)specie_temp_1
             out_energy_params%species_obst(specie_temp_1,i_V)=1
          end do
       end do
    end if

    !skip "V_spectwist" Label
    read(unit_input,*)

    if (out_energy_params%V_spectwist_type .ne. 0) then
       !for each pairwise interaction, read the interacting pair and fill the symmetric entries
       do i_V=1,N_spectwist
          read(unit_input,*)specie_temp_1,temp_spectwist_type
          if (temp_spectwist_type .eq. 0) then
             read(unit_input,*)
             temp_N_spectwist_params=1
             cycle
          elseif (temp_spectwist_type .eq. 1) then
             temp_N_spectwist_params=1
          end if
          !print *,shape(out_energy_params%spectwist_active)
          out_energy_params%spectwist_active(1,specie_temp_1)=temp_spectwist_type
          out_energy_params%spectwist_active(2,specie_temp_1)=temp_N_spectwist_params
          read(unit_input,*)out_energy_params%spectwist_params(1:temp_N_spectwist_params,specie_temp_1)
       end do
    end if

    !print *,out_energy_params%spectwist_params(1,:)

    !skip "V_pairwise" Label
    read(unit_input,*)

    if (out_energy_params%V_pairwise_type .ne. 0) then
       !for each pairwise interaction, read the interacting pair and fill the symmetric entries
       do i_V=1,N_pairs
          read(unit_input,*)specie_temp_1,specie_temp_2,temp_pair_type
          if (temp_pair_type .eq. 0) then
             read(unit_input,*)
          elseif (temp_pair_type .eq. 1) then
             temp_N_pair_params=1
          elseif (temp_pair_type .eq. 2) then
             temp_N_pair_params=1
          elseif (temp_pair_type .eq. 3) then
             temp_N_pair_params=2
          end if
          out_energy_params%pairs_active(1,specie_temp_1,specie_temp_2)=temp_pair_type
          out_energy_params%pairs_active(2,specie_temp_1,specie_temp_2)=temp_N_pair_params
          read(unit_input,*)out_energy_params%pair_params(1:temp_N_pair_params,specie_temp_1,specie_temp_2)
          if (specie_temp_1 .ne. specie_temp_2) then
             out_energy_params%pair_params(:,specie_temp_2,specie_temp_1)=&
                  out_energy_params%pair_params(:,specie_temp_1,specie_temp_2)
             out_energy_params%pairs_active(:,specie_temp_2,specie_temp_1)=&
                  out_energy_params%pairs_active(:,specie_temp_1,specie_temp_2)
          end if
       end do
    end if


    if (allocated(temp_params)) deallocate(temp_params)

    close(unit_input) 
    
  end subroutine read_energy_file

  
  !subroutine to write the coordinates on the integer lattice
  subroutine write_x(in_ring,unit_x)

    implicit none

    type(ring), intent(in) :: in_ring
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_node, i_dim
    
    format_string="(I1,'(A1,(I4))')"
    write(x_string,format_string)(in_ring%dim-1)
    x_string='(I4,'//trim(x_string)//')'

    do i_node=1,in_ring%N_nodes
       write(unit_x,x_string)in_ring%x(1,i_node),(",",in_ring%x(i_dim,i_node),i_dim=2,in_ring%dim)
    end do
    
  end subroutine write_x

  !subroutine to write the species and coordinates on the integer lattice
  subroutine write_species_x(in_ring,unit_x)

    implicit none

    type(ring), intent(in) :: in_ring
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_node, i_dim
    
    format_string="(I1,'(A1,(I6))')"
    write(x_string,format_string)(in_ring%dim+1)
    x_string='(I6,'//trim(x_string)//')'

    do i_node=1,in_ring%N_nodes
       write(unit_x,x_string)i_node,&
            ",",in_ring%species(i_node),&
            (",",in_ring%x(i_dim,i_node),i_dim=1,in_ring%dim)
    end do
    
  end subroutine write_species_x

  !subroutine to write the coordinates on the integer lattice
  subroutine write_obst_x(in_obstacle_distribution,unit_x)

    implicit none

    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_obst, i_dim
    
    format_string="(I1,'(A1,(I6))')"
    write(x_string,format_string)(in_obstacle_distribution%dim)
    x_string='(I6,A1,I3,'//trim(x_string)//')'

    do i_obst=1,in_obstacle_distribution%obst_base%N
       write(unit_x,x_string)i_obst,",",in_obstacle_distribution%obst_base%obst_types(i_obst),&
            (",",in_obstacle_distribution%obst_base%x(i_dim,i_obst),i_dim=1,in_obstacle_distribution%dim)
    end do
    
  end subroutine write_obst_x

    !subroutine to write the coordinates on the integer lattice
  subroutine write_obst_x_old(in_obstacle_distribution,unit_x)

    implicit none

    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_obst, i_dim
    
    format_string="(I1,'(A1,(I6))')"
    write(x_string,format_string)(in_obstacle_distribution%dim)
    x_string='(I6,'//trim(x_string)//')'

    do i_obst=1,in_obstacle_distribution%N
       write(unit_x,x_string)i_obst,&
            (",",in_obstacle_distribution%x(i_dim,i_obst),i_dim=1,in_obstacle_distribution%dim)
    end do
    
  end subroutine write_obst_x_old

  !subroutine to write the obstacle index and coordinates on the integer lattice
  subroutine write_obst_placed_x(in_obstacle_distribution,unit_x)

    implicit none

    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_obst, i_dim
    
    format_string="(I1,'(A1,(I6))')"
    write(x_string,format_string)(in_obstacle_distribution%dim)
    x_string='(I6,A1,I3,'//trim(x_string)//')'

    do i_obst=1,in_obstacle_distribution%obst_base%N
       if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
          write(unit_x,x_string)i_obst,",",in_obstacle_distribution%obst_base%obst_types(i_obst),&
               (",",in_obstacle_distribution%obst_base%x(i_dim,i_obst),i_dim=1,in_obstacle_distribution%dim)
       end if
    end do
    
  end subroutine write_obst_placed_x

    !subroutine to write the obstacle index and coordinates on the integer lattice
  subroutine write_obst_placed_x_old(in_obstacle_distribution,unit_x)

    implicit none

    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    integer, intent(in) :: unit_x

    character(len=120) :: format_string, x_string
    integer :: i_obst, i_dim
    
    format_string="(I1,'(A1,(I6))')"
    write(x_string,format_string)(in_obstacle_distribution%dim)
    x_string='(I6,'//trim(x_string)//')'

    do i_obst=1,in_obstacle_distribution%N_placed
       write(unit_x,x_string)in_obstacle_distribution%placed(i_obst),&
            (",",in_obstacle_distribution%x(i_dim,in_obstacle_distribution%placed(i_obst)),&
            i_dim=1,in_obstacle_distribution%dim)
    end do
    
  end subroutine write_obst_placed_x_old

  !subroutine to write the bounds
  subroutine write_bounds(in_ring,unit_bounds)

    implicit none

    type(ring), intent(in) :: in_ring
    integer, intent(in) :: unit_bounds

    character(len=120) :: format_string, bounds_string
    integer :: N_entries
    integer :: i_bound, i_entry


    if (in_ring%bound_type .ne. 0) then
       format_string="(I1,'(A1,(I6))')"

       N_entries=1
       if (in_ring%bound_type .eq. 1) then
          N_entries=2
       elseif (in_ring%bound_type .eq. 2) then
          N_entries=in_ring%dim
       end if
       
       write(bounds_string,format_string)(N_entries)
       bounds_string='(I6,'//trim(bounds_string)//')'

       do i_bound=1,in_ring%N_bounds
          write(unit_bounds,bounds_string)in_ring%bounds(0,i_bound),&
               (",",in_ring%bounds(i_entry,i_bound),i_entry=1,N_entries)
       end do
    end if
    
  end subroutine write_bounds

  !subroutine to write a restart file
  subroutine write_restart_file(in_ring,in_obstacle_distribution,&
       in_deform_params,in_grow_params,&
       unit_restart,final_flag)

    implicit none

    type(ring), intent(in) :: in_ring
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(grow_params), intent(in) :: in_grow_params
    type(deform_params), intent(in) :: in_deform_params
    integer, intent(in) :: unit_restart
    logical, intent(in) :: final_flag

    integer :: i_node, i_obst

    write(unit_restart,*)"##### PROGRAM DIRECTIVES #####"
    write(unit_restart,*)"run = 0 (0=no,1=serial,2=parallel)"
    write(unit_restart,*)"output_dir = "
    write(unit_restart,*)"output_label = "
    write(unit_restart,*)"##### GROWTH PARAMETERS #####"
    if (final_flag .eqv. .true.) then
       write(unit_restart,*)"grow_flag = ",0
    else
       write(unit_restart,*)"grow_flag = ",in_grow_params%grow_flag
    end if
    write(unit_restart,*)"growth_target = ",in_grow_params%growth_target
    write(unit_restart,*)"final_species = ",trim(in_grow_params%species_file)
    write(unit_restart,*)"grow_energy_file = ",trim(in_grow_params%energy_file)
    write(unit_restart,*)"max_growth = ",in_grow_params%max_growth
    write(unit_restart,*)"initial_growth = ",in_grow_params%initial_growth
    write(unit_restart,*)"rand_init = ",in_grow_params%rand_init
    write(unit_restart,*)"ref_node = ",in_grow_params%ref_node
    write(unit_restart,"(A,ES20.10)")" alpha_1 = ",in_grow_params%alpha(1)
    write(unit_restart,"(A,ES20.10)")" eta_1 = ",in_grow_params%eta(1)
    write(unit_restart,*)" min_iter_move = ",in_grow_params%min_iter_move
    write(unit_restart,"(A,ES20.10)")" alpha_2 = ",in_grow_params%alpha(2)
    write(unit_restart,"(A,ES20.10)")" eta_2 = ",in_grow_params%eta(2)
    write(unit_restart,*)" max_wait = ",in_grow_params%max_wait
    write(unit_restart,*)"max_iter = ",in_grow_params%max_iter
    write(unit_restart,*)"err_check_freq = ",in_grow_params%err_check_freq
    write(unit_restart,*)"##### DEFORMATION PARAMETERS #####"
    if (final_flag .eqv. .true.) then
       write(unit_restart,*)"min_rep = ",in_deform_params%min_rep," (previous run)"
       write(unit_restart,*)"max_rep = ",in_deform_params%max_rep," (previous run)"
       write(unit_restart,*)"ergodic_flag = ",in_deform_params%ergodic_flag," (previous run)"
       write(unit_restart,*)"move_energy_file = ",&
            trim(in_deform_params%energy_file)," (previous run)"
       write(unit_restart,*)"num_iter = ",in_deform_params%num_iter," (previous run)"
       write(unit_restart,*)"ref_init = ",in_deform_params%ref_init," (previous run)"
       write(unit_restart,*)"ref_wait = ",in_deform_params%ref_wait," (previous run)"
       write(unit_restart,*)"ref_freq = ",in_deform_params%ref_freq," (previous run)"
       write(unit_restart,*)"ref_mult = ",in_deform_params%ref_mult," (previous run)"
       write(unit_restart,*)"kink_mult = ",in_deform_params%kink_mult," (previous run)"
       write(unit_restart,*)"kink_repeat = ",in_deform_params%kink_repeat," (previous run)"
       write(unit_restart,*)"obst_freq = ",in_deform_params%obst_freq," (previous run)"
       write(unit_restart,*)"obst_mult = ",in_deform_params%obst_mult," (previous run)"
       write(unit_restart,*)"max_iter = ",in_deform_params%max_iter," (previous run)"
       write(unit_restart,*)"err_check_freq = ",in_deform_params%err_check_freq," (previous run)"
    else
       write(unit_restart,*)"min_rep = ",in_deform_params%min_rep
       write(unit_restart,*)"max_rep = ",in_deform_params%max_rep
       write(unit_restart,*)"ergodicity_flag = ",in_deform_params%ergodic_flag
       write(unit_restart,*)"move_energy_file = ",trim(in_deform_params%energy_file)
       write(unit_restart,*)"num_iter = ",in_deform_params%num_iter," (iterations remaining)"
       write(unit_restart,*)"ref_init = ",in_deform_params%ref_init
       write(unit_restart,*)"ref_wait = ",in_deform_params%ref_wait
       write(unit_restart,*)"ref_freq = ",in_deform_params%ref_freq
       write(unit_restart,*)"ref_mult = ",in_deform_params%ref_mult
       write(unit_restart,*)"kink_mult = ",in_deform_params%kink_mult
       write(unit_restart,*)"kink_repeat = ",in_deform_params%kink_repeat
       write(unit_restart,*)"obst_freq = ",in_deform_params%obst_freq
       write(unit_restart,*)"obst_mult = ",in_deform_params%obst_mult
       write(unit_restart,*)"max_iter = ",in_deform_params%max_iter
       write(unit_restart,*)"err_check_freq = ",in_deform_params%err_check_freq
    end if
    write(unit_restart,*)"##### SYSTEM DESCRIPTION #####"
    write(unit_restart,*)"dim = ",in_ring%dim
    write(unit_restart,*)"N_nodes = ",in_ring%N_nodes
    write(unit_restart,*)"N_fixed = ",in_ring%N_fixed
    write(unit_restart,*)"N_species = ",in_ring%N_species
    write(unit_restart,*)"bound_type = ",in_ring%bound_type
    write(unit_restart,*)"N_bounds = ",in_ring%N_bounds
    write(unit_restart,*)"obst_type = ",in_obstacle_distribution%obst_type
    write(unit_restart,*)"N_obst = ",in_obstacle_distribution%obst_base%N
    write(unit_restart,*)"N_placed = ",in_obstacle_distribution%obst_base%N_placed
    
    write(unit_restart,*)"##### bounds #####"
    call write_bounds(in_ring,unit_restart)

    write(unit_restart,*)"##### placed obstacles #####"
    if (in_obstacle_distribution%N_placed .ne. 0) then
       do i_obst=1,in_obstacle_distribution%obst_base%N
          if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
             write(unit_restart,"(i5)")i_obst !something is wrong with this statement
          end if
       end do
    end if

    write(unit_restart,*)"##### obstacle coordinates #####"
    call write_obst_x(in_obstacle_distribution,unit_restart)
    
    write(unit_restart,*)"##### fixed nodes #####"
    if (in_ring%N_fixed .ne. 0) then
       do i_node=1,in_ring%N_fixed
          write(unit_restart,"(i5)")in_ring%fixed_nodes(i_node) !something is wrong with this statement
       end do
    end if
    
    write(unit_restart,*)"##### id, species, and node coordinates #####"
    call write_species_x(in_ring,unit_restart)
    
    
  end subroutine write_restart_file


end module read_write
