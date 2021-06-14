module ring_grow

  !this module contains the subroutines and functions used to generate
  !configurations of the polymer ring
  
  use params
  use ring_objects
  use read_write
  use ring_procedures
  use ring_config
  use ring_polymer

  implicit none

contains

  pure function equil_wait(alpha,eta,deltaN)

    implicit none

    integer, intent(in) :: deltaN
    real(rp), intent(in) :: alpha, eta

    integer :: equil_wait

    equil_wait=int(alpha*(1.0-exp(-deltaN/eta)))
    
  end function equil_wait

  pure function equil_num_iter(alpha,eta,N)

    implicit none

    integer, intent(in) :: N
    real(rp), intent(in) :: alpha, eta

    integer :: equil_num_iter

    equil_num_iter=int(alpha*exp(-N/eta))

  end function equil_num_iter
  
  subroutine propose_grow(node_pair,x_pairs,max_growth,grow_flag,in_ring)

    implicit none

    type(ring), intent(in) :: in_ring
    integer, intent(in) :: max_growth, grow_flag
    
    integer(ip), intent(out) :: x_pairs(1:in_ring%dim,1:2*max_growth,1:2*(in_ring%dim-1))
    integer, intent(out) :: node_pair(1:2)

    logical :: accepted_first, fixed_flag
    
    integer :: i_fixed, i_pair, i_dim, i_grow
    integer :: dim_count
    integer :: matched_dims(1:in_ring%dim-1), mapping(1:2*(in_ring%dim-1))
    integer(ip) :: temp_x_pairs(1:in_ring%dim,1:2*max_growth,2*(in_ring%dim-1))

    accepted_first=.false.

    do while (accepted_first .eqv. .false.)
    
       !randomly choose first node
       if (grow_flag .eq. 5) then
          node_pair(1)=int_rand(0,1)
          if (node_pair(1) .eq. 0) then
             node_pair(1)=int_rand(1,in_ring%N_nodes/4)
             !node_pair(1)=int_rand(1,in_ring%N_nodes/4-1)
          else
             node_pair(1)=int_rand(in_ring%N_nodes-in_ring%N_nodes/4+1,in_ring%N_nodes)
          end if
       else
          node_pair(1)=int_rand(1,in_ring%N_nodes)
       end if

       !test if the first node is a fixed node or not
       if (in_ring%fixed .eqv. .true.) then
          
          fixed_flag=.false.
          do i_fixed=1,in_ring%N_fixed
             if (node_pair(1) .eq. in_ring%fixed_nodes(i_fixed)) then
                fixed_flag=.true.
                exit
             end if
          end do

          if (fixed_flag .eqv. .false.) then
             accepted_first=.true.
          end if
          
       else
          accepted_first=.true.
       end if

    end do
    
    !randomly choose forward or backward direction for the second node

    if (int_rand(0,1) .eq. 0) then
       if (node_pair(1) .eq. in_ring%N_nodes) then
          node_pair(2)=node_pair(1)
          node_pair(1)=1
       else
          node_pair(2)=node_pair(1)+1
       end if
    else
       if (node_pair(1) .eq. 1) then
          node_pair(2)=in_ring%N_nodes
       else
          node_pair(2)=node_pair(1)
          node_pair(1)=node_pair(1)-1
       end if
    end if

    !find matching dimensions
    dim_count=1
    do i_dim=1,in_ring%dim
       if (in_ring%x(i_dim,node_pair(1)) .eq. in_ring%x(i_dim,node_pair(2))) then
          matched_dims(dim_count)=i_dim
          dim_count=dim_count+1
       end if
    end do

    !print *,matched_dims

    !generate the new positions of the two nodes
    do i_pair=1,2*(in_ring%dim-1)
       x_pairs(:,1,i_pair)=in_ring%x(:,node_pair(1))
       x_pairs(:,2*max_growth,i_pair)=in_ring%x(:,node_pair(2))
       mapping(i_pair)=i_pair
       i_dim=(i_pair+1)/2
       if (mod(i_pair,2) .ne. 1) then
          x_pairs(matched_dims(i_dim),1,i_pair)=x_pairs(matched_dims(i_dim),1,i_pair)+1
          x_pairs(matched_dims(i_dim),2*max_growth,i_pair)=x_pairs(matched_dims(i_dim),2*max_growth,i_pair)+1
          if (max_growth .gt. 1) then
             do i_grow=2,max_growth
                x_pairs(:,i_grow,i_pair)=x_pairs(:,i_grow-1,i_pair)
                x_pairs(:,2*max_growth-i_grow+1,i_pair)=x_pairs(:,2*max_growth-i_grow+2,i_pair)
                x_pairs(matched_dims(i_dim),i_grow,i_pair)=&
                     x_pairs(matched_dims(i_dim),i_grow,i_pair)+1
                x_pairs(matched_dims(i_dim),2*max_growth-i_grow+1,i_pair)=&
                     x_pairs(matched_dims(i_dim),2*max_growth-i_grow+1,i_pair)+1
             end do
          end if
       else
          x_pairs(matched_dims(i_dim),1,i_pair)=x_pairs(matched_dims(i_dim),1,i_pair)-1
          x_pairs(matched_dims(i_dim),2*max_growth,i_pair)=x_pairs(matched_dims(i_dim),2*max_growth,i_pair)-1
          if (max_growth .gt. 1) then
             do i_grow=2,max_growth
                x_pairs(:,i_grow,i_pair)=x_pairs(:,i_grow-1,i_pair)
                x_pairs(:,2*max_growth-i_grow+1,i_pair)=x_pairs(:,2*max_growth-i_grow+2,i_pair)
                x_pairs(matched_dims(i_dim),i_grow,i_pair)=&
                     x_pairs(matched_dims(i_dim),i_grow,i_pair)-1
                x_pairs(matched_dims(i_dim),2*max_growth-i_grow+1,i_pair)=&
                     x_pairs(matched_dims(i_dim),2*max_growth-i_grow+1,i_pair)-1
             end do
          end if
       end if
    end do

    !randomly permute the mapping such that there is no bias
    call generate_permutation(mapping,2*(in_ring%dim-1))

    temp_x_pairs=x_pairs

    do i_pair=1,2*(in_ring%dim-1)
       x_pairs(:,:,i_pair)=temp_x_pairs(:,:,mapping(i_pair))
    end do
       
  end subroutine propose_grow

  subroutine grow_ring(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_energy_params,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params
    type(deform_params), intent(in) :: in_deform_params
    type(grow_params), intent(in) :: in_grow_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(ring_subset) :: ring_growth
    type(deform_params) :: equil_deform_params
    
    logical :: grow_success, move_flag
    logical :: length_error, intersect_error

    integer :: unit_rep
    integer :: prev_equil, wait_count
    integer :: write_update_freq
    integer :: initial_length, test_length, N_nodes_init
    integer :: count, i_growth, grow_size

    integer :: node_pair(1:2)
    integer :: config_test(1:2)
    integer(ip) :: x_pairs(1:in_ring%dim,1:2*in_grow_params%max_growth,1:2*(in_ring%dim-1))

    integer :: ref_node, ref_node_loc, ref_node_sep
    integer :: i_node
    integer :: ref_coord(1:in_ring%dim)

    character(len=120) :: format_string

    sim_completed=.false.

    write_update_freq=1
    N_nodes_init=in_ring%N_nodes
    
    if ((in_grow_params%growth_target-N_nodes_init) .le. 100) then
       write_update_freq=1
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 500) then
       write_update_freq=5
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 1000) then
       write_update_freq=10
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 5000) then
       write_update_freq=50
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 10000) then
       write_update_freq=100
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 50000) then
       write_update_freq=500
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 100000) then
       write_update_freq=1000
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 1000000) then
       write_update_freq=10000
    end if

    !write_update_freq=1
    !print *,"write_update_freq=",write_update_freq

    format_string="(a,i5.5,a,i8,a,i8,a)"
    
    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    equil_deform_params=new_deform_params(in_deform_params%energy_file,&
         in_deform_params%ergodic_flag,1,&
         0,0,in_deform_params%ref_freq,1,&
         1,in_deform_params%kink_repeat,&
         0,0,&
         in_deform_params%max_iter,in_deform_params%err_check_freq,&
         in_deform_params%min_rep,in_deform_params%max_rep)

    prev_equil=N_nodes_init
    wait_count=min(in_grow_params%max_wait,&
         equil_wait(in_grow_params%alpha(2),in_grow_params%eta(2),0))

    initial_length=ring_length(in_ring)

    if ((in_grow_params%grow_flag .eq. 4) .or.&
         (in_grow_params%grow_flag .eq. 5)) then
       if (in_grow_params%ref_node .eq. -1) then
          ref_node=int_rand(1,in_grow_params%growth_target)
       else
          ref_node=in_grow_params%ref_node
       end if
       ref_coord=in_ring%x(:,in_ring%N_nodes/2)
       !write(unit_log,"(a,i5.5,a,i3.3)")'[grow_ring] --- i_rep=',i_rep,' --- ref_specie=',&
       !in_grow_params%species_final(ref_node)
       !print *,ref_node
    end if

    write(unit_log,"(a,i5.5,a)")'[grow_ring] --- i_rep=',i_rep,' --- BEGIN'
    
    ring_growth=new_ring_subset(in_ring%dim,2*in_grow_params%max_growth)
    
    !loop while the ring has not grown to its full size
    do while (in_ring%N_nodes .lt. in_grow_params%growth_target)
       
       count=0
       grow_success=.false.
       move_flag=.false.

       !loop and attempt to find an acceptable growth
       do while ((grow_success .eqv. .false.) .and.&
            (count .lt. in_grow_params%max_iter))

          !propose a grow
          call propose_grow(node_pair,x_pairs,&
               in_grow_params%max_growth,in_grow_params%grow_flag,&
               in_ring)

          !set the ring subset equal to the grow
          ring_growth%node_indices(1:2)=node_pair

          grow_size=min(in_grow_params%max_growth,(in_grow_params%growth_target-in_ring%N_nodes)/2)

          do while ((grow_success .eqv. .false.) .and.&
               (grow_size .gt. 0))

             i_growth=1

             ring_growth%N=2*grow_size

             !loop over the growth directions and test them
             do while ((grow_success .eqv. .false.) .and.&
                  (i_growth .le. 2*(in_ring%dim-1)))

                ring_growth%x(:,1:grow_size)=&
                     x_pairs(:,1:grow_size,i_growth)
                ring_growth%x(:,grow_size+1:2*grow_size)=&
                     x_pairs(:,2*in_grow_params%max_growth-grow_size+1:2*in_grow_params%max_growth,i_growth)

                call perform_grow(grow_success,in_ring,ring_growth)

                i_growth=i_growth+1

                count=count+1

             end do

             grow_size=grow_size-1

          end do

       end do

       !if no growth was possible, then move the ring for a bit before attempting
       !new growth options
       if (grow_success .eqv. .false.) then

          move_flag=.true.
          count=0

       elseif (in_grow_params%grow_flag .gt. 1) then

          count=0

          if (((in_ring%N_nodes-N_nodes_init) .gt. in_grow_params%initial_growth) .and.&
               (((in_ring%N_nodes-prev_equil) .gt. wait_count) .or.&
               (in_ring%N_nodes .eq. in_grow_params%growth_target))) then

             move_flag=.true.
             prev_equil=in_ring%N_nodes
             wait_count=min(in_grow_params%max_wait,&
                  equil_wait(in_grow_params%alpha(2),in_grow_params%eta(2),&
                  in_ring%N_nodes-N_nodes_init))

          end if

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then

             !shift array to match reference node
             ref_node_loc=-1
             i_node=1
             do while ((ref_node_loc .eq. -1) .and.&
                  (i_node .lt. in_ring%N_nodes+1))
                if (sum(abs(ref_coord-in_ring%x(:,i_node))) .eq. 0) then
                   ref_node_loc=i_node
                else
                   i_node=i_node+1
                end if
             end do
             ref_node_sep=in_ring%N_nodes/2-ref_node_loc
             !print *,'ref_coord=',ref_coord
             in_ring%x=cshift(in_ring%x,-ref_node_sep,2)
             !print *,'shifted_coord=',in_ring%x(:,in_ring%N_nodes/2)

          end if
          
       end if

       if (move_flag .eqv. .true.) then

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then

             !assign species
             call update_species(ref_node,in_ring%N_nodes,&
                  in_grow_params%growth_target,in_grow_params%species_final,in_ring%species)

             !write(unit_log,"(a,i5.5,a,i3.3)")'[grow_ring_obstacles] --- i_rep=',i_rep,' --- specie_new=',&
             !in_ring%species(in_ring%N_nodes/2)

             ! write(unit_log,*)in_ring%N_nodes
             ! do i_node=1,in_ring%N_nodes
             !    write(unit_log,*)i_node,in_ring%species(i_node)
             ! end do

          end if

          write(unit_log,format_string)' [grow_ring] --- i_rep=',i_rep,&
               ' --- progress: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes - move'

          equil_deform_params%num_iter=max(in_grow_params%min_iter_move,&
               equil_num_iter(in_grow_params%alpha(1),&
               in_grow_params%eta(1),in_ring%N_nodes))
          equil_deform_params%kink_mult=in_ring%N_nodes/3
          equil_deform_params%ref_mult=max(1,equil_deform_params%kink_mult/10)
          equil_deform_params%err_check_freq=equil_deform_params%ref_mult

          !conditional for MH sampling if there are potentials
          if ((in_energy_params%V_nn_type .ne. 0) .or.&
               (in_energy_params%V_bend_type .ne. 0) .or.&
               (in_energy_params%V_twist_type .ne. 0) .or.&
               (in_energy_params%V_coord_type .ne. 0) .or.&
               (in_energy_params%V_obst_type .ne. 0) .or.&
               (in_energy_params%V_pairwise_type .ne. 0)) then
          
             call move_ring_MH(sim_completed,&
                  in_ring,in_obstacle_distribution,&
                  in_energy_params,&
                  equil_deform_params,in_grow_params,&
                  restart_file,i_rep)

          else

             call move_ring(sim_completed,&
                  in_ring,in_obstacle_distribution,&
                  equil_deform_params,in_grow_params,&
                  restart_file,i_rep)

          end if

          move_flag=.false.

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then
             ref_coord=in_ring%x(:,in_ring%N_nodes/2)
          end if

       end if

       if ((mod(in_ring%N_nodes-N_nodes_init,in_grow_params%err_check_freq) .eq. 0) .or.&
            (in_ring%N_nodes .eq. in_grow_params%growth_target)) then

          length_error=.true.
          intersect_error=.true.

          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if ((test_length-initial_length) .ne. (in_ring%N_nodes-N_nodes_init)) then
             write(unit_log,format_string)' [grow_ring] --- i_rep=',i_rep,&
               ' --- ERROR at: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'
             write(unit_log,"(a,i5.5,a,i6)")'  [grow_ring] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [grow_ring] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [grow_ring] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,format_string)' [grow_ring] --- i_rep=',i_rep,&
               ' --- ERROR at: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'
             write(unit_log,"(a,i5.5,a)")' [grow_ring] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [grow_ring] --- i_rep=',i_rep,&
                  ' --- interesection at nodes',&
                  config_test(1),' and ',config_test(2)
             exit
          else
             intersect_error=.false.
          end if

          !if the error checks are passed, then write a restart file
          if ((length_error .eqv. .false.) .and.&
               (intersect_error .eqv. .false.)) then

             open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  in_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (in_ring%N_nodes .eq. in_grow_params%growth_target) then
                sim_completed=.true.

                if ((in_grow_params%grow_flag .eq. 4) .or.&
                     (in_grow_params%grow_flag .eq. 5)) then

                   ! write(unit_log,*)ref_node
                   ! !write(unit_log,*)in_ring%N_nodes
                   ! do i_node=1,in_ring%N_nodes
                   !    write(unit_log,*)i_node,in_ring%species(i_node)
                   ! end do
                   ref_node_sep=ref_node-in_ring%N_nodes/2+1
                   !write(unit_log,*)'ref_node_sep=',ref_node_sep
                   in_ring%x=cshift(in_ring%x,-ref_node_sep,2)
                   in_ring%species=cshift(in_ring%species,-ref_node_sep,1)
                   ! do i_node=1,in_ring%N_nodes
                   !    write(unit_log,*)i_node,in_ring%species(i_node)
                   ! end do

                end if
                
             end if

          end if

       end if !end conditional for error checking

       if ((mod(in_ring%N_nodes-N_nodes_init,write_update_freq) .eq. 0) .or.&
            ((in_grow_params%growth_target-in_ring%N_nodes) .eq. 0)) then
          write(unit_log,format_string)' [grow_ring] --- i_rep=',i_rep,&
               ' --- progress: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'

       end if

    end do

    write(unit_log,"(a,i5.5,a)")'[grow_ring] --- i_rep=',i_rep,' --- END'

  end subroutine grow_ring

  subroutine grow_ring_obstacles(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_energy_params,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    
    type(obstacle_distribution), intent(inout) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params
    type(deform_params), intent(in) :: in_deform_params
    type(grow_params), intent(in) :: in_grow_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(ring_subset) :: ring_growth
    type(deform_params) :: equil_deform_params
    
    logical :: grow_success, move_flag
    logical :: length_error, intersect_error, obstacle_error

    integer :: unit_rep
    integer :: prev_equil, wait_count
    integer :: write_update_freq
    integer :: initial_length, test_length, N_nodes_init
    integer :: count, i_growth, grow_size

    integer :: node_pair(1:2)
    integer :: config_test(1:2), node_obst_test(1:2)
    integer(ip) :: x_pairs(1:in_ring%dim,1:2*in_grow_params%max_growth,1:2*(in_ring%dim-1))

    integer :: ref_node, ref_node_loc, ref_node_sep
    integer :: i_node
    integer :: ref_coord(1:in_ring%dim)

    character(len=120) :: format_string

    sim_completed=.false.

    write_update_freq=1
    N_nodes_init=in_ring%N_nodes
    
    if ((in_grow_params%growth_target-N_nodes_init) .le. 100) then
       write_update_freq=1
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 500) then
       write_update_freq=5
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 1000) then
       write_update_freq=10
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 5000) then
       write_update_freq=50
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 10000) then
       write_update_freq=100
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 50000) then
       write_update_freq=500
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 100000) then
       write_update_freq=1000
    elseif ((in_grow_params%growth_target-N_nodes_init) .le. 1000000) then
       write_update_freq=10000
    end if

    !write_update_freq=1
    !print *,"write_update_freq=",write_update_freq

    format_string="(a,i5.5,a,i8,a,i8,a)"
    
    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    equil_deform_params=new_deform_params(in_deform_params%energy_file,&
         in_deform_params%ergodic_flag,1,&
         0,0,in_deform_params%ref_freq,1,&
         1,in_deform_params%kink_repeat,&
         0,0,&
         in_deform_params%max_iter,in_deform_params%err_check_freq,&
         in_deform_params%min_rep,in_deform_params%max_rep)

    prev_equil=N_nodes_init
    wait_count=min(in_grow_params%max_wait,&
         equil_wait(in_grow_params%alpha(2),in_grow_params%eta(2),0))

    initial_length=ring_length(in_ring)

    if ((in_grow_params%grow_flag .eq. 4) .or.&
         (in_grow_params%grow_flag .eq. 5)) then
       if (in_grow_params%ref_node .eq. -1) then
          ref_node=int_rand(1,in_grow_params%growth_target)
       else
          ref_node=in_grow_params%ref_node
       end if
       ref_coord=in_ring%x(:,in_ring%N_nodes/2)
       !write(unit_log,"(a,i5.5,a,i3.3)")'[grow_ring_obstacles] --- i_rep=',i_rep,' --- ref_specie=',&
            !in_grow_params%species_final(ref_node)
       !print *,ref_node
    end if

    write(unit_log,"(a,i5.5,a)")'[grow_ring_obstacles] --- i_rep=',i_rep,' --- BEGIN'
    
    ring_growth=new_ring_subset(in_ring%dim,2*in_grow_params%max_growth)
    
    !loop while the ring has not grown to its full size
    do while (in_ring%N_nodes .lt. in_grow_params%growth_target)
       
       count=0
       grow_success=.false.
       move_flag=.false.
       
       !loop and attempt to find an acceptable growth
       do while ((grow_success .eqv. .false.) .and.&
            (count .lt. in_grow_params%max_iter))

          !propose a grow
          call propose_grow(node_pair,x_pairs,&
               in_grow_params%max_growth,in_grow_params%grow_flag,&
               in_ring)
          
          !set the ring subset equal to the grow
          ring_growth%node_indices(1:2)=node_pair

          grow_size=min(in_grow_params%max_growth,(in_grow_params%growth_target-in_ring%N_nodes)/2)

          do while ((grow_success .eqv. .false.) .and.&
               (grow_size .gt. 0))
          
             i_growth=1

             ring_growth%N=2*grow_size

             !loop over the growth directions and test them
             do while ((grow_success .eqv. .false.) .and.&
                  (i_growth .le. 2*(in_ring%dim-1)))

                ring_growth%x(:,1:grow_size)=&
                     x_pairs(:,1:grow_size,i_growth)
                ring_growth%x(:,grow_size+1:2*grow_size)=&
                     x_pairs(:,2*in_grow_params%max_growth-grow_size+1:2*in_grow_params%max_growth,i_growth)

                call perform_grow_obstacles(grow_success,in_ring,ring_growth,&
                     in_obstacle_distribution)

                i_growth=i_growth+1

                count=count+1

             end do

             grow_size=grow_size-1

          end do
          
       end do

       !if no growth was possible, then move the ring for a bit before attempting
       !new growth options
       if (grow_success .eqv. .false.) then

          move_flag=.true.
          count=0

       elseif (in_grow_params%grow_flag .gt. 1) then

          count=0

          if (((in_ring%N_nodes-N_nodes_init) .gt. in_grow_params%initial_growth) .and.&
               (((in_ring%N_nodes-prev_equil) .gt. wait_count) .or.&
               (in_ring%N_nodes .eq. in_grow_params%growth_target))) then

             move_flag=.true.
             prev_equil=in_ring%N_nodes
             wait_count=min(in_grow_params%max_wait,&
                  equil_wait(in_grow_params%alpha(2),in_grow_params%eta(2),&
                  in_ring%N_nodes-N_nodes_init))

          end if

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then

             !shift array to match reference node
             ref_node_loc=-1
             i_node=1
             do while ((ref_node_loc .eq. -1) .and.&
                  (i_node .lt. in_ring%N_nodes+1))
                if (sum(abs(ref_coord-in_ring%x(:,i_node))) .eq. 0) then
                   ref_node_loc=i_node
                else
                   i_node=i_node+1
                end if
             end do
             ref_node_sep=in_ring%N_nodes/2-ref_node_loc
             !print *,'ref_coord=',ref_coord
             in_ring%x=cshift(in_ring%x,-ref_node_sep,2)
             !print *,'shifted_coord=',in_ring%x(:,in_ring%N_nodes/2)
             
          end if
          
       end if

       if (move_flag .eqv. .true.) then

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then
             
             !assign species
             call update_species(ref_node,in_ring%N_nodes,&
                  in_grow_params%growth_target,in_grow_params%species_final,in_ring%species)

             !write(unit_log,"(a,i5.5,a,i3.3)")'[grow_ring_obstacles] --- i_rep=',i_rep,' --- specie_new=',&
                  !in_ring%species(in_ring%N_nodes/2)

             ! write(unit_log,*)in_ring%N_nodes
             ! do i_node=1,in_ring%N_nodes
             !    write(unit_log,*)i_node,in_ring%species(i_node)
             ! end do
             
          end if

          write(unit_log,format_string)' [grow_ring_obstacles] --- i_rep=',i_rep,&
               ' --- progress: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes - move'

          equil_deform_params%num_iter=max(in_grow_params%min_iter_move,&
               equil_num_iter(in_grow_params%alpha(1),&
               in_grow_params%eta(1),in_ring%N_nodes))
          equil_deform_params%kink_mult=in_ring%N_nodes/3
          equil_deform_params%ref_mult=max(1,equil_deform_params%kink_mult/10)
          equil_deform_params%err_check_freq=equil_deform_params%ref_mult

          !conditional for MH sampling if there are potentials
          if ((in_energy_params%V_nn_type .ne. 0) .or.&
               (in_energy_params%V_bend_type .ne. 0) .or.&
               (in_energy_params%V_twist_type .ne. 0) .or.&
               (in_energy_params%V_coord_type .ne. 0) .or.&
               (in_energy_params%V_obst_type .ne. 0) .or.&
               (in_energy_params%V_pairwise_type .ne. 0)) then

             call move_ring_obstacles_MH(sim_completed,&
                  in_ring,in_obstacle_distribution,&
                  in_energy_params,&
                  equil_deform_params,in_grow_params,&
                  restart_file,i_rep)

          else

             call move_ring_obstacles(sim_completed,&
                  in_ring,in_obstacle_distribution,&
                  equil_deform_params,in_grow_params,&
                  restart_file,i_rep)

          end if

          move_flag=.false.

          if ((in_grow_params%grow_flag .eq. 4) .or.&
               (in_grow_params%grow_flag .eq. 5)) then
             ref_coord=in_ring%x(:,in_ring%N_nodes/2)
          end if

       end if

       if ((mod(in_ring%N_nodes-N_nodes_init,in_grow_params%err_check_freq) .eq. 0) .or.&
            (in_ring%N_nodes .eq. in_grow_params%growth_target)) then

          length_error=.true.
          intersect_error=.true.

          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if ((test_length-initial_length) .ne. (in_ring%N_nodes-N_nodes_init)) then
             write(unit_log,format_string)' [grow_ring_obstacles] --- i_rep=',i_rep,&
               ' --- ERROR at: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'
             write(unit_log,"(a,i5.5,a,i6)")'  [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,format_string)' [grow_ring_obstacles] --- i_rep=',i_rep,&
               ' --- ERROR at: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'
             write(unit_log,"(a,i5.5,a)")' [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- interesection at nodes',&
                  config_test(1),' and ',config_test(2)
             exit
          else
             intersect_error=.false.
          end if

          !test for obstacle intersections
          node_obst_test=obstacle_test_fast(in_ring,in_obstacle_distribution)
          
          if (sum(node_obst_test) .ne. -2) then
             write(unit_log,format_string)' [grow_ring_obstacles] --- i_rep=',i_rep,&
               ' --- ERROR at: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'
             write(unit_log,"(a,i5.5,a)")' [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring-obstacle intersection'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [grow_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- interesection of obstacle ',&
                  node_obst_test(2),' by node ',node_obst_test(1)
             exit
          else
             obstacle_error=.false.
          end if

          !if the error checks are passed, then write a restart file
          if ((length_error .eqv. .false.) .and.&
               (intersect_error .eqv. .false.) .and.&
               (obstacle_error .eqv. .false.)) then

             open(UNIT=unit_rep,FILE=trim(restart_file),ACTION='write',STATUS='replace')

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  in_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (in_ring%N_nodes .eq. in_grow_params%growth_target) then
                sim_completed=.true.

                if ((in_grow_params%grow_flag .eq. 4) .or.&
                     (in_grow_params%grow_flag .eq. 5)) then

                   ! write(unit_log,*)ref_node
                   ! !write(unit_log,*)in_ring%N_nodes
                   ! do i_node=1,in_ring%N_nodes
                   !    write(unit_log,*)i_node,in_ring%species(i_node)
                   ! end do
                   ref_node_sep=ref_node-in_ring%N_nodes/2+1
                   !write(unit_log,*)'ref_node_sep=',ref_node_sep
                   in_ring%x=cshift(in_ring%x,-ref_node_sep,2)
                   in_ring%species=cshift(in_ring%species,-ref_node_sep,1)
                   ! do i_node=1,in_ring%N_nodes
                   !    write(unit_log,*)i_node,in_ring%species(i_node)
                   ! end do

                end if

             end if

          end if

       end if !end conditional for error checking

       if ((mod(in_ring%N_nodes-N_nodes_init,write_update_freq) .eq. 0) .or.&
            ((in_grow_params%growth_target-in_ring%N_nodes) .eq. 0)) then
          write(unit_log,format_string)' [grow_ring_obstacles] --- i_rep=',i_rep,&
               ' --- progress: ',&
               in_ring%N_nodes,' of ',in_grow_params%growth_target,' nodes'

       end if

    end do

    write(unit_log,"(a,i5.5,a)")'[grow_ring_obstacles] --- i_rep=',i_rep,' --- END'

  end subroutine grow_ring_obstacles
  

end module ring_grow
