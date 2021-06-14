module ring_config

  !this module contains the subroutines and functions used to analyze the ring configuration

  use params
  use ring_objects
  use ring_procedures
  use energy_terms

  implicit none

contains
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! CONFIGURATION TEST ROUTINES AND FUNCTIONS !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !function to calculate the length of the ring
  pure function ring_length(in_ring)

    implicit none

    type(ring), intent(in) :: in_ring
    integer :: i_node
    integer :: ring_length

    ring_length=0

    !calculate the ring length using a Manhattan distance (L1-norm)
    do i_node=1,in_ring%N_nodes-1
       ring_length=ring_length+sum(abs(in_ring%x(:,i_node+1)-in_ring%x(:,i_node)))
    end do
    ring_length=ring_length+sum(abs(in_ring%x(:,1)-in_ring%x(:,in_ring%N_nodes)))

  end function ring_length

  !function to test if a configuration is valid
  pure function valid_config(in_ring)

    implicit none

    type(ring), intent(in) :: in_ring
    integer :: i_node, j_node, i_dim
    integer :: count
    integer :: valid_config(1:2)

    valid_config=(/-1,-1/)

    !iterate over all possible pairs of nodes and test for intersections

    !iterate up the next to last node
    do i_node=1,(in_ring%N_nodes-1)
       !iterate from i_node+1 up to the last node
       do j_node=(i_node+1),in_ring%N_nodes

          !test for an intersection
          count=0
          do i_dim=1,in_ring%dim
             if (in_ring%x(i_dim,i_node) .eq. in_ring%x(i_dim,j_node)) then
                count=count+1
             end if
          end do

          !if an intersection is detected, then indicate the intersecting nodes
          !and break out of the inner loop
          if (count .eq. in_ring%dim) then
             valid_config=(/i_node,j_node/)
             return
             !exit
          end if

       end do !finish iterating over the N_nodes-i_node nodes

       !if an intersection is detected, then break out of the outer loop
       ! if (sum(valid_config) .ne. -2) then
       !    exit
       ! end if

    end do !finish iterating over N_nodes-1 nodes

  end function valid_config

  !function to test if the ring intersects any of the obstacles
  pure function obstacle_test(in_ring,in_obstacle_distribution)

    implicit none

    type(ring), intent(in) :: in_ring
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    integer :: i_node, i_obst
    integer :: obstacle_test(1:2)


    obstacle_test=(/-1,-1/)

    do i_obst=1,in_obstacle_distribution%N

       if (in_obstacle_distribution%placed(i_obst) .eq. 1) then
       
          do i_node=1,in_ring%N_nodes

             if (sum(abs(in_ring%x(:,i_node)-in_obstacle_distribution%x(:,i_obst))) .eq. 0) then
                obstacle_test=(/i_node,i_obst/)
                return
                !exit
             end if

          end do !end loop over nodes

       end if

    end do !end loop over obstacles

  end function obstacle_test

  
  !function to test if the ring intersects any of the obstacles
  function obstacle_test_fast(in_ring,in_obstacle_distribution)

    implicit none

    type(ring), intent(in) :: in_ring
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    logical :: obstacle_flag
    integer :: temp_obst
    integer :: i_node
    integer :: obstacle_test_fast(1:2)


    obstacle_test_fast=(/-1,-1/)


    obstacle_flag=.false.
    do i_node=1,in_ring%N_nodes

       call octree_search(obstacle_flag,temp_obst,in_ring%x(:,i_node),&
            in_obstacle_distribution%N,in_obstacle_distribution%x,&
            in_obstacle_distribution%obst_tree%root_node)

       if (obstacle_flag .eqv. .true.) then
          obstacle_test_fast=(/i_node,temp_obst/)
          return
       end if

    end do
       

  end function obstacle_test_fast


  !function to test single intersections
  function test_intersection_sing(in_ring,ignore_index,test_coord)

    implicit none

    type(ring), intent(in) :: in_ring
    type(coords), intent(in) :: test_coord
    integer, intent(in) :: ignore_index
    integer :: i_node, i_dim, count_match
    logical :: test_intersection_sing

    test_intersection_sing=.false.

    do i_node=1,in_ring%N_nodes
       if (i_node .ne. ignore_index) then
          count_match=0
          do i_dim=1,in_ring%dim
             if (test_coord%x(i_dim,1) .eq. in_ring%x(i_dim,i_node)) then
                count_match=count_match+1
             end if
          end do
          if (count_match .eq. in_ring%dim) then
             test_intersection_sing=.true.
             exit
          end if
       end if
    end do

  end function test_intersection_sing

  pure subroutine test_intersection_growth(intersect_flag,in_ring,N,x)

    implicit none

    logical, intent(out) :: intersect_flag

    type(ring), intent(in) :: in_ring
    integer, intent(in) :: N
    integer(ip), intent(in) :: x(1:in_ring%dim,1:N)
    
    integer :: i_test, i_node, i_dim, count_match

    intersect_flag=.false.
    
    do i_test=1,N
       do i_node=1,in_ring%N_nodes
          count_match=0
          do i_dim=1,in_ring%dim
             if (in_ring%x(i_dim,i_node) .eq.&
                  x(i_dim,i_test)) then
                count_match=count_match+1
             end if
          end do
          if (count_match .eq. in_ring%dim) then
             intersect_flag=.true.
             return
             !exit
          end if
       end do
    end do
    
  end subroutine test_intersection_growth

  !function to test multiple intersections for the motion of a subset of the ring
  pure subroutine test_intersection_mult(intersect_flag,in_ring,in_ring_subset)

    implicit none

    logical, intent(out) :: intersect_flag

    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_ring_subset

    integer :: i_node, i_dim, i_test, i_ignore, count_match
    integer :: prev_node_index

    intersect_flag=.false.

    !iterate over the nodes in the subset representing a motion of the ring
    do i_ignore=1,in_ring_subset%N

       !if an intersection has been detected break out of the loop
       ! if (intersect_flag .eqv. .true.) then
       !    exit
       ! else

       !if at the first node in the ring subset,
       !then iterate through the ring up to the first node and test intersections
       if (i_ignore .eq. 1) then

          if (in_ring_subset%node_indices(i_ignore) .gt. 1) then

             do i_node=1,in_ring_subset%node_indices(i_ignore)-1
                ! if (intersect_flag .eqv. .true.) then
                !    exit
                ! else
                do i_test=1,in_ring_subset%N
                   count_match=0
                   do i_dim=1,in_ring%dim
                      if (in_ring%x(i_dim,i_node) .eq.&
                           in_ring_subset%x(i_dim,i_test)) then
                         count_match=count_match+1
                      end if
                   end do
                   if (count_match .eq. in_ring%dim) then
                      intersect_flag=.true.
                      return
                      !exit
                   end if
                end do
                ! end if
             end do

          end if

          !if the first node in the ring subset is also the only node,
          !then iterate through the rest of ring and test the intersections
          if ((in_ring_subset%N .eq. 1) .and.&
               (in_ring_subset%node_indices(i_ignore) .lt. in_ring%N_nodes)) then

             do i_node=in_ring_subset%node_indices(i_ignore)+1,in_ring%N_nodes
                ! if (intersect_flag .eqv. .true.) then
                !    exit
                ! else
                do i_test=1,in_ring_subset%N
                   count_match=0
                   do i_dim=1,in_ring%dim
                      if (in_ring%x(i_dim,i_node) .eq.&
                           in_ring_subset%x(i_dim,i_test)) then
                         count_match=count_match+1
                      end if
                   end do
                   if (count_match .eq. in_ring%dim) then
                      intersect_flag=.true.
                      return
                      !exit
                   end if
                end do
                ! end if
             end do

          end if

          !if at the last node in the ring subset,
          !then iterate up the to the last node and test intersections
       elseif (i_ignore .eq. in_ring_subset%N) then

          if (in_ring_subset%node_indices(i_ignore)-prev_node_index .gt. 1) then

             do i_node=prev_node_index+1,in_ring_subset%node_indices(i_ignore)-1
                ! if (intersect_flag .eqv. .true.) then
                !    exit
                ! else
                do i_test=1,in_ring_subset%N
                   count_match=0
                   do i_dim=1,in_ring%dim
                      if (in_ring%x(i_dim,i_node) .eq.&
                           in_ring_subset%x(i_dim,i_test)) then
                         count_match=count_match+1
                      end if
                   end do
                   if (count_match .eq. in_ring%dim) then
                      intersect_flag=.true.
                      return
                      !exit
                   end if
                end do
                ! end if
             end do

          end if

          !if the last node in the ring subset is less than the number of nodes in the ring,
          !then iterate to the last node in the ring and test intersection
          if (in_ring_subset%node_indices(i_ignore) .lt. in_ring%N_nodes) then

             do i_node=in_ring_subset%node_indices(i_ignore)+1,in_ring%N_nodes
                ! if (intersect_flag .eqv. .true.) then
                !    exit
                ! else
                do i_test=1,in_ring_subset%N
                   count_match=0
                   do i_dim=1,in_ring%dim
                      if (in_ring%x(i_dim,i_node) .eq.&
                           in_ring_subset%x(i_dim,i_test)) then
                         count_match=count_match+1
                      end if
                   end do
                   if (count_match .eq. in_ring%dim) then
                      intersect_flag=.true.
                      return
                      !exit
                   end if
                end do
                ! end if
             end do

          end if

       else

          !if the current node and the previous node in the ring subset are separated
          !by more than one space, then iterate over the nodes between them and
          !test for intersections
          if (in_ring_subset%node_indices(i_ignore)-prev_node_index .gt. 1) then

             do i_node=prev_node_index+1,in_ring_subset%node_indices(i_ignore)-1
                ! if (intersect_flag .eqv. .true.) then
                !    exit
                ! else
                do i_test=1,in_ring_subset%N
                   count_match=0
                   do i_dim=1,in_ring%dim
                      if (in_ring%x(i_dim,i_node) .eq.&
                           in_ring_subset%x(i_dim,i_test)) then
                         count_match=count_match+1
                      end if
                   end do
                   if (count_match .eq. in_ring%dim) then
                      intersect_flag=.true.
                      return
                      !exit
                   end if
                end do
                ! end if
             end do

          end if

       end if

       prev_node_index=in_ring_subset%node_indices(i_ignore)

       !end if !end conditional used to break the loop

    end do !finish the iteration over the nodes in the ring subset

  end subroutine test_intersection_mult

  !function to test multiple intersections for the motion of a subset of the ring
  pure subroutine test_intersection_mult_deltaE(intersect_flag,deltaE,&
       in_ring,in_ring_subset,in_energy_params)

    implicit none

    logical, intent(out) :: intersect_flag
    real(rp), intent(out) :: deltaE
    
    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_ring_subset
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, i_test, i_ignore
    integer :: prev_node_index, sep

    real(rp) :: old_energy, new_energy

    intersect_flag=.false.

    deltaE=0.0d0

    !iterate over the nodes in the subset representing a motion of the ring
    do i_ignore=1,in_ring_subset%N

       !if an intersection has been detected break out of the loop
       if (intersect_flag .eqv. .true.) then
          exit
       else

          !if at the first node in the ring subset,
          !then iterate through the ring up to the first node and test intersections
          if (i_ignore .eq. 1) then

             if (in_ring_subset%node_indices(i_ignore) .gt. 1) then

                do i_node=1,in_ring_subset%node_indices(i_ignore)-1
                   ! if (intersect_flag .eqv. .true.) then
                   !    exit
                   ! else
                   do i_test=1,in_ring_subset%N
                      sep=sum(abs(in_ring%x(:,i_node)-in_ring_subset%x(:,i_test)))
                      if (sep .eq. 0) then
                         intersect_flag=.true.
                         return
                         !exit
                      else

                         if (in_energy_params%V_pairwise_type .ne. 0) then
                            new_energy=V_pairwise(in_ring%dim,sep,&
                                 in_ring%species(i_node),&
                                 in_ring%species(in_ring_subset%node_indices(i_test)),&
                                 in_ring%x(:,i_node),in_ring_subset%x(:,i_test),&
                                 in_energy_params)
                         end if

                         if ((sep .eq. 1) .and.&
                              (in_energy_params%V_nn_type .eq. 1)) then
                            new_energy=new_energy+in_energy_params%nn_params(1)
                         end if
                         
                      end if
                      
                      sep=sum(abs(in_ring%x(:,i_node)-&
                           in_ring%x(:,in_ring_subset%node_indices(i_test))))

                      if (in_energy_params%V_pairwise_type .ne. 0) then
                         old_energy=V_pairwise(in_ring%dim,sep,&
                              in_ring%species(i_node),&
                              in_ring%species(in_ring_subset%node_indices(i_test)),&
                              in_ring%x(:,i_node),in_ring%x(:,in_ring_subset%node_indices(i_test)),&
                              in_energy_params)
                      end if

                      if ((sep .eq. 1) .and.&
                           (in_energy_params%V_nn_type .eq. 1)) then
                         old_energy=old_energy+in_energy_params%nn_params(1)
                      end if

                      
                      deltaE=deltaE+new_energy
                      deltaE=deltaE-old_energy
                   end do
                   ! end if
                end do

             end if

             !if the first node in the ring subset is also the only node,
             !then iterate through the rest of ring and test the intersections
             if ((in_ring_subset%N .eq. 1) .and.&
                  (in_ring_subset%node_indices(i_ignore) .lt. in_ring%N_nodes)) then

                do i_node=in_ring_subset%node_indices(i_ignore)+1,in_ring%N_nodes
                   ! if (intersect_flag .eqv. .true.) then
                   !    exit
                   ! else
                   do i_test=1,in_ring_subset%N
                      sep=sum(abs(in_ring%x(:,i_node)-in_ring_subset%x(:,i_test)))
                      if (sep .eq. 0) then
                         intersect_flag=.true.
                         return
                         !exit
                      else

                         if (in_energy_params%V_pairwise_type .ne. 0) then
                            new_energy=V_pairwise(in_ring%dim,sep,&
                                 in_ring%species(i_node),&
                                 in_ring%species(in_ring_subset%node_indices(i_test)),&
                                 in_ring%x(:,i_node),in_ring_subset%x(:,i_test),&
                                 in_energy_params)
                         end if

                         if ((sep .eq. 1) .and.&
                              (in_energy_params%V_nn_type .eq. 1)) then
                            new_energy=new_energy+in_energy_params%nn_params(1)
                         end if
                      
                      end if
                      
                      sep=sum(abs(in_ring%x(:,i_node)-&
                           in_ring%x(:,in_ring_subset%node_indices(i_test))))

                      if (in_energy_params%V_pairwise_type .ne. 0) then
                         old_energy=V_pairwise(in_ring%dim,sep,&
                              in_ring%species(i_node),&
                              in_ring%species(in_ring_subset%node_indices(i_test)),&
                              in_ring%x(:,i_node),in_ring%x(:,in_ring_subset%node_indices(i_test)),&
                              in_energy_params)
                      end if

                      if ((sep .eq. 1) .and.&
                           (in_energy_params%V_nn_type .eq. 1)) then
                         old_energy=old_energy+in_energy_params%nn_params(1)
                      end if
                      
                      deltaE=deltaE+new_energy
                      deltaE=deltaE-old_energy
                   end do
                   ! end if
                end do

             end if

             !if at the last node in the ring subset,
             !then iterate up the to the last node and test intersections
          elseif (i_ignore .eq. in_ring_subset%N) then

             if (in_ring_subset%node_indices(i_ignore)-prev_node_index .gt. 1) then

                do i_node=prev_node_index+1,in_ring_subset%node_indices(i_ignore)-1
                   ! if (intersect_flag .eqv. .true.) then
                   !    exit
                   ! else
                   do i_test=1,in_ring_subset%N
                      sep=sum(abs(in_ring%x(:,i_node)-in_ring_subset%x(:,i_test)))
                      if (sep .eq. 0) then
                         intersect_flag=.true.
                         return
                         !exit
                      else

                         if (in_energy_params%V_pairwise_type .ne. 0) then
                            new_energy=V_pairwise(in_ring%dim,sep,&
                                 in_ring%species(i_node),&
                                 in_ring%species(in_ring_subset%node_indices(i_test)),&
                                 in_ring%x(:,i_node),in_ring_subset%x(:,i_test),&
                                 in_energy_params)
                         end if
                         
                         if ((sep .eq. 1) .and.&
                              (in_energy_params%V_nn_type .eq. 1)) then
                            new_energy=new_energy+in_energy_params%nn_params(1)
                         end if
                         
                      end if
                      
                      sep=sum(abs(in_ring%x(:,i_node)-&
                           in_ring%x(:,in_ring_subset%node_indices(i_test))))

                      if (in_energy_params%V_pairwise_type .ne. 0) then
                         old_energy=V_pairwise(in_ring%dim,sep,&
                              in_ring%species(i_node),&
                              in_ring%species(in_ring_subset%node_indices(i_test)),&
                              in_ring%x(:,i_node),in_ring%x(:,in_ring_subset%node_indices(i_test)),&
                              in_energy_params)
                      end if

                      if ((sep .eq. 1) .and.&
                           (in_energy_params%V_nn_type .eq. 1)) then
                         old_energy=old_energy+in_energy_params%nn_params(1)
                      end if
                      
                      deltaE=deltaE+new_energy
                      deltaE=deltaE-old_energy
                   end do
                   ! end if
                end do

             end if

             !if the last node in the ring subset is less than the number of nodes in the ring,
             !then iterate to the last node in the ring and test intersection
             if (in_ring_subset%node_indices(i_ignore) .lt. in_ring%N_nodes) then

                do i_node=in_ring_subset%node_indices(i_ignore)+1,in_ring%N_nodes
                   ! if (intersect_flag .eqv. .true.) then
                   !    exit
                   ! else
                   do i_test=1,in_ring_subset%N
                      sep=sum(abs(in_ring%x(:,i_node)-in_ring_subset%x(:,i_test)))
                      if (sep .eq. 0) then
                         intersect_flag=.true.
                         return
                         !exit
                      else

                         if (in_energy_params%V_pairwise_type .ne. 0) then
                            new_energy=V_pairwise(in_ring%dim,sep,&
                                 in_ring%species(i_node),&
                                 in_ring%species(in_ring_subset%node_indices(i_test)),&
                                 in_ring%x(:,i_node),in_ring_subset%x(:,i_test),&
                                 in_energy_params)
                         end if

                         if ((sep .eq. 1) .and.&
                              (in_energy_params%V_nn_type .eq. 1)) then
                            new_energy=new_energy+in_energy_params%nn_params(1)
                         end if
                         
                      end if
                      
                      sep=sum(abs(in_ring%x(:,i_node)-&
                           in_ring%x(:,in_ring_subset%node_indices(i_test))))

                      if (in_energy_params%V_pairwise_type .ne. 0) then
                         old_energy=V_pairwise(in_ring%dim,sep,&
                              in_ring%species(i_node),&
                              in_ring%species(in_ring_subset%node_indices(i_test)),&
                              in_ring%x(:,i_node),in_ring%x(:,in_ring_subset%node_indices(i_test)),&
                              in_energy_params)
                      end if
                         
                      if ((sep .eq. 1) .and.&
                           (in_energy_params%V_nn_type .eq. 1)) then
                         old_energy=old_energy+in_energy_params%nn_params(1)
                      end if
                      
                      deltaE=deltaE+new_energy
                      deltaE=deltaE-old_energy
                   end do
                   ! end if
                end do

             end if

          else

             !if the current node and the previous node in the ring subset are separated
             !by more than one space, then iterate over the nodes between them and
             !test for intersections
             if (in_ring_subset%node_indices(i_ignore)-prev_node_index .gt. 1) then

                do i_node=prev_node_index+1,in_ring_subset%node_indices(i_ignore)-1
                   ! if (intersect_flag .eqv. .true.) then
                   !    exit
                   ! else
                   do i_test=1,in_ring_subset%N
                      sep=sum(abs(in_ring%x(:,i_node)-in_ring_subset%x(:,i_test)))
                      if (sep .eq. 0) then
                         intersect_flag=.true.
                         return
                         !exit
                      else

                         if (in_energy_params%V_pairwise_type .ne. 0) then
                            new_energy=V_pairwise(in_ring%dim,sep,&
                                 in_ring%species(i_node),&
                                 in_ring%species(in_ring_subset%node_indices(i_test)),&
                                 in_ring%x(:,i_node),in_ring_subset%x(:,i_test),&
                                 in_energy_params)
                         end if

                         if ((sep .eq. 1) .and.&
                              (in_energy_params%V_nn_type .eq. 1)) then
                            new_energy=new_energy+in_energy_params%nn_params(1)
                         end if
                         
                      end if
                      
                      sep=sum(abs(in_ring%x(:,i_node)-&
                           in_ring%x(:,in_ring_subset%node_indices(i_test))))

                      if (in_energy_params%V_pairwise_type .ne. 0) then
                         old_energy=V_pairwise(in_ring%dim,sep,&
                              in_ring%species(i_node),&
                              in_ring%species(in_ring_subset%node_indices(i_test)),&
                              in_ring%x(:,i_node),in_ring%x(:,in_ring_subset%node_indices(i_test)),&
                              in_energy_params)
                      end if

                      if ((sep .eq. 1) .and.&
                           (in_energy_params%V_nn_type .eq. 1)) then
                         old_energy=old_energy+in_energy_params%nn_params(1)
                      end if
                      
                      deltaE=deltaE+new_energy
                      deltaE=deltaE-old_energy
                   end do
                   ! end if
                end do

             end if

          end if

          prev_node_index=in_ring_subset%node_indices(i_ignore)

       end if !end conditional used to break the loop

    end do !finish the iteration over the nodes in the ring subset

  end subroutine test_intersection_mult_deltaE


  !function to test the boundaries
  pure subroutine test_bounds(bounds_flag,in_ring,in_ring_subset)

    implicit none

    logical, intent(out) :: bounds_flag

    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_ring_subset

    integer :: i_node, i_bound
    real(rp) :: r_sqrd, lower_lim, upper_lim

    bounds_flag=.false.

    !test bounds for box in Cartesian space
    if (in_ring%bound_type .eq. 1) then
       do i_node=1,in_ring_subset%N
          do i_bound=1,in_ring%N_bounds
             if ((in_ring_subset%x(in_ring%bounds(0,i_bound),i_node) .lt.&
                  in_ring%bounds(1,i_bound)) .or.&
                  (in_ring_subset%x(in_ring%bounds(0,i_bound),i_node) .gt.&
                  in_ring%bounds(2,i_bound))) then
                bounds_flag=.true.
                return
                !exit
             end if
          end do
          ! if (bounds_flag .eqv. .true.) then
          !    exit
          ! end if
       end do
       !test bounds for the interior of a sphere
    elseif (in_ring%bound_type .eq. 2) then
       r_sqrd=in_ring%bounds(0,1)**2.0d0
       do i_node=1,in_ring_subset%N
          if (floor(sum((in_ring_subset%x(:,i_node)-in_ring%bounds(1:3,1))**2.0d0)) .gt.&
               r_sqrd) then
             bounds_flag=.true.
             return
             !exit
          end if
       end do
    elseif (in_ring%bound_type .eq. 3) then
       r_sqrd=in_ring%bounds(-1,1)**2.0d0
       lower_lim=in_ring%bounds(3,1)-in_ring%bounds(0,0)/2.0d0
       upper_lim=in_ring%bounds(3,1)+in_ring%bounds(0,0)/2.0d0
       do i_node=1,in_ring_subset%N
          if ((floor(sum((in_ring_subset%x(1:2,i_node)-in_ring%bounds(1:2,1))**2.0d0)) .gt.&
               r_sqrd) .or.&
               (in_ring_subset%x(3,i_node) .le. lower_lim) .or.&
               (in_ring_subset%x(3,i_node) .ge. upper_lim)) then
             bounds_flag=.true.
             return
             !exit
          end if
       end do
    end if

  end subroutine test_bounds

  !function to test the boundaries
  pure subroutine test_bounds_deltaE(bounds_flag,deltaE,&
       in_ring,in_ring_subset,in_energy_params)

    implicit none

    logical, intent(out) :: bounds_flag
    real(rp), intent(out) :: deltaE

    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_ring_subset
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, i_bound
    real(rp) :: r_sqrd
    real(rp) :: old_energy, new_energy

    bounds_flag=.false.
    deltaE=0.0d0

    !test bounds for box in Cartesian space
    if (in_ring%bound_type .eq. 1) then
       do i_node=1,in_ring_subset%N
          do i_bound=1,in_ring%N_bounds
             if ((in_ring_subset%x(in_ring%bounds(0,i_bound),i_node) .lt.&
                  in_ring%bounds(1,i_bound)) .or.&
                  (in_ring_subset%x(in_ring%bounds(0,i_bound),i_node) .gt.&
                  in_ring%bounds(2,i_bound))) then
                bounds_flag=.true.
                return
                !exit
             end if
          end do
          ! if (bounds_flag .eqv. .true.) then
          !    exit
          ! else
          old_energy=V_coord(in_ring%dim,&
               in_ring%species(in_ring_subset%node_indices(i_node)),&
               in_ring%x(:,in_ring_subset%node_indices(i_node)),&
               in_ring,in_energy_params)
          new_energy=V_coord(in_ring%dim,&
               in_ring%species(in_ring_subset%node_indices(i_node)),&
               in_ring_subset%x(:,i_node),&
               in_ring,in_energy_params)
          deltaE=deltaE+new_energy
          deltaE=deltaE-old_energy
          ! end if
       end do
       !test bounds for the interior of a sphere
    elseif (in_ring%bound_type .eq. 2) then
       r_sqrd=in_ring%bounds(0,1)**2.0d0
       do i_node=1,in_ring_subset%N
          if (floor(sum((in_ring_subset%x(:,i_node)-in_ring%bounds(1:3,1))**2.0d0)) .gt.&
               r_sqrd) then
             bounds_flag=.true.
             return
             !exit
          else
             old_energy=V_coord(in_ring%dim,&
                  in_ring%species(in_ring_subset%node_indices(i_node)),&
                  in_ring%x(:,in_ring_subset%node_indices(i_node)),&
                  in_ring,in_energy_params)
             new_energy=V_coord(in_ring%dim,&
                  in_ring%species(in_ring_subset%node_indices(i_node)),&
                  in_ring_subset%x(:,i_node),&
                  in_ring,in_energy_params)
             deltaE=deltaE+new_energy
             deltaE=deltaE-old_energy
          end if
       end do
    end if

  end subroutine test_bounds_deltaE

  !function to test obstacle-ring intersections when placing the obstacles
  pure function test_obstacle_place(in_obstacle_distribution,in_ring,i_obst)

    implicit none

    integer, intent(in) :: i_obst
    type(ring), intent(in) :: in_ring
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    logical :: test_obstacle_place

    integer :: i_node

    test_obstacle_place=.false.


    do i_node=1,in_ring%N_nodes

       if (sum(abs(in_ring%x(:,i_node)-in_obstacle_distribution%x(:,i_obst))) .eq. 0) then
          test_obstacle_place=.true.
          return
       end if

    end do

  end function test_obstacle_place

  !subroutine to test obstacle-ring intersections when moving the ring - can be pure if remove octree
  pure subroutine test_obstacle_move(obstacle_flag,in_obstacle_distribution,in_ring_subset)

    implicit none

    logical, intent(out) :: obstacle_flag

    type(ring_subset), intent(in) :: in_ring_subset
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    integer :: i_node, i_obst

    obstacle_flag=.false.

    do i_obst=1,in_obstacle_distribution%N

       if (in_obstacle_distribution%placed(i_obst) .eq. 1) then
       
          do i_node=1,in_ring_subset%N

             if (sum(abs(in_ring_subset%x(:,i_node)-in_obstacle_distribution%x(:,i_obst))) .eq. 0) then
                obstacle_flag=.true.
                return
             end if

          end do !end loop over nodes

       end if

    end do !end loop over obstacles


  end subroutine test_obstacle_move

  !subroutine to test obstacle-ring intersections when moving the ring - can be pure if remove octree
  subroutine test_obstacle_move_fast(obstacle_flag,in_obstacle_distribution,in_ring_subset)

    implicit none

    logical, intent(out) :: obstacle_flag

    type(ring_subset), intent(in) :: in_ring_subset
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    integer :: temp_obst
    integer :: i_node

    obstacle_flag=.false.

    do i_node=1,in_ring_subset%N

       call octree_search(obstacle_flag,temp_obst,in_ring_subset%x(:,i_node),&
            in_obstacle_distribution%N,in_obstacle_distribution%x,&
            in_obstacle_distribution%obst_tree%root_node)

       if (obstacle_flag .eqv. .true.) then
          return
       end if

    end do       

  end subroutine test_obstacle_move_fast

  !subroutine to test obstacle-ring intersections when moving the ring - can be made pure if octree is removed
  subroutine test_obstacle_move_deltaE(obstacle_flag,deltaE,&
       in_obstacle_distribution,in_ring_subset,in_energy_params)

    implicit none

    logical, intent(out) :: obstacle_flag
    real(rp), intent(out) :: deltaE

    type(ring_subset), intent(in) :: in_ring_subset
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params

    integer :: temp_point
    integer :: i_node

    deltaE=0.0d0
    
    obstacle_flag=.false.


    do i_node=1,in_ring_subset%N

       call octree_search(obstacle_flag,temp_point,in_ring_subset%x(:,i_node),&
            in_obstacle_distribution%N,in_obstacle_distribution%x,&
            in_obstacle_distribution%obst_tree%root_node)

       if (obstacle_flag .eqv. .true.) then
          return
       end if

    end do


  end subroutine test_obstacle_move_deltaE

  !subroutine to calculate the change in bending energy
  pure subroutine calc_bend_deltaE(deltaE,in_ring,&
       in_move,in_energy_params)

    implicit none

    real(rp), intent(out) :: deltaE

    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_move
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, i_trip, lo, hi
    integer :: triplet(1:3,1:2), dx_pair(1:in_ring%dim,2)

    deltaE=0.0d0

    !kink move
    if (in_move%N .eq. 1) then

       lo=1
       hi=in_move%N

    !reflect move
    else
       
       !move must wrap around the ring, crossing the origin
       if (((in_move%node_indices(1) .eq. 1) .and.&
            (in_move%node_indices(in_move%N) .eq. in_ring%N_nodes)) .and.&
            (in_move%node_indices(in_move%N) .lt. in_ring%N_nodes)) then

          do i_node=1,in_move%N-1
             if ((in_move%node_indices(i_node+1)-in_move%node_indices(i_node))&
                  .gt. 1) then
                lo=i_node+1
                hi=i_node
                exit
             end if
          end do

       !move is a continuous segment not crossing the origin
       else
          
          lo=1
          hi=in_move%N
          
       end if

    end if

    triplet(2,1)=in_move%node_indices(lo)-1
    triplet(2,2)=in_move%node_indices(hi)+1

    !find the triplets describing the coordinates about which the bending occurs
    do i_trip=1,2
       if (triplet(2,i_trip) .eq. 0) then
          triplet(2,i_trip)=in_ring%N_nodes
       elseif(triplet(2,i_trip) .eq. in_ring%N_nodes+1) then
          triplet(2,i_trip)=1
       end if
    end do

    do i_trip=1,2
       if (triplet(2,i_trip) .eq. 1) then
          triplet(1,i_trip)=in_ring%N_nodes
          triplet(3,i_trip)=triplet(2,i_trip)+1
       elseif (triplet(2,i_trip) .eq. in_ring%N_nodes) then
          triplet(1,i_trip)=triplet(2,i_trip)-1
          triplet(3,i_trip)=1
       else
          triplet(1,i_trip)=triplet(2,i_trip)-1
          triplet(3,i_trip)=triplet(2,i_trip)+1
       end if
    end do

    !calculate the old energy
    do i_trip=1,2
       dx_pair(:,1)=in_ring%x(:,triplet(2,i_trip))-in_ring%x(:,triplet(1,i_trip))
       dx_pair(:,2)=in_ring%x(:,triplet(3,i_trip))-in_ring%x(:,triplet(2,i_trip))
       deltaE=deltaE-V_bend(in_ring%dim,dx_pair,in_energy_params)
    end do

    !calculate the new energy
    do i_trip=1,2
       if (i_trip .eq. 1) then
          dx_pair(:,1)=in_ring%x(:,triplet(2,i_trip))-in_ring%x(:,triplet(2,i_trip))
          dx_pair(:,2)=in_move%x(:,lo)-in_ring%x(:,triplet(2,i_trip))
       else
          dx_pair(:,1)=in_ring%x(:,triplet(2,i_trip))-in_move%x(:,hi)
          dx_pair(:,2)=in_ring%x(:,triplet(3,i_trip))-in_ring%x(:,triplet(2,i_trip))
       end if
       deltaE=deltaE+V_bend(in_ring%dim,dx_pair,in_energy_params)
    end do
    
  end subroutine calc_bend_deltaE

  !subroutine to calculate the change in twisting energy
  subroutine calc_twist_deltaE(deltaE,in_ring,&
       in_move,in_energy_params)

    implicit none

    real(rp), intent(out) :: deltaE

    type(ring), intent(in) :: in_ring
    type(ring_subset), intent(in) :: in_move
    type(energy_params), intent(in) :: in_energy_params

    integer :: temp_twist, temp_pos
    integer :: i, i_node, i_side, i_quad, lo, hi
    integer :: quad(1:4,1:3,1:2)
    integer :: p_tens(-3:3,-3:3)
    integer(ip) :: dx_triplet(1:in_ring%dim,3), temp_coords(1:in_ring%dim,4)

    real(rp) :: temp_params(1:3,1:2)

    deltaE=0.0d0
    temp_params=0.0d0

    p_tens=perm_tensor()

    !kink move
    if (in_move%N .eq. 1) then

       lo=1
       hi=in_move%N

       !reflect move
    else

       !move must wrap around the ring, crossing the origin
       if ((in_move%node_indices(1) .eq. 1) .and.&
            (in_move%node_indices(in_move%N) .eq. in_ring%N_nodes)) then

          do i_node=1,in_move%N-1
             if ((in_move%node_indices(i_node+1)-in_move%node_indices(i_node))&
                  .gt. 1) then
                lo=i_node+1
                hi=i_node
                exit
             end if
          end do

       !move is a continuous segment not crossing the origin
       else

          lo=1
          hi=in_move%N

       end if

    end if

    ! print *,lo
    ! print *,hi
    ! print *,in_move%node_indices(lo)-1
    ! print *,in_move%node_indices(hi)+1

    quad(1,1,1)=in_move%node_indices(lo)-1
    quad(2,2,1)=in_move%node_indices(lo)-1
    quad(3,3,1)=in_move%node_indices(lo)-1
    quad(4,1,2)=in_move%node_indices(hi)+1
    quad(3,2,2)=in_move%node_indices(hi)+1
    quad(2,3,2)=in_move%node_indices(hi)+1

    !find the quadruplets describing the coordinates about which the bending occurs

    !low side
    do i_quad=1,3
       do i=1,4
          quad(i,i_quad,1)=quad(i_quad,i_quad,1)+i-i_quad
       end do
       !print *,quad(:,i_quad,1)
    end do

    !high side
    do i_quad=1,3
       do i=1,4
          quad(i,i_quad,2)=quad(4-i_quad+1,i_quad,2)+i-(4-i_quad+1)
       end do
       !print *,quad(:,i_quad,2)
    end do

    do i_side=1,2
       do i_quad=1,3
          do i=1,4
             if (quad(i,i_quad,i_side) .lt. 1) then
                quad(i,i_quad,i_side)=quad(i,i_quad,i_side)+in_ring%N_nodes
             elseif (quad(i,i_quad,i_side) .gt. in_ring%N_nodes) then
                quad(i,i_quad,i_side)=quad(i,i_quad,i_side)-in_ring%N_nodes
             end if
          end do
          !print *,quad(:,i_quad,i_side)
       end do
    end do
    
    !calculate the old energy
    do i_side=1,2
       do i_quad=1,3
          temp_twist=0
          do i=1,3
             dx_triplet(:,i)=in_ring%x(:,quad(i+1,i_quad,i_side))-in_ring%x(:,quad(i,i_quad,i_side))
             !print *,dx_triplet(:,i)
          end do
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_old=",temp_twist

          temp_params(i_quad,i_side)=calc_twist_param(in_ring%N_nodes,in_ring%species,&
               in_energy_params,quad(:,i_quad,i_side))

          ! print *,quad(:,i_quad,i_side)
          ! print *,in_ring%species(quad(:,i_quad,i_side))
          ! print *,temp_params(i_quad,i_side)
          
          deltaE=deltaE-(in_energy_params%twist_params(1)+&
               temp_params(i_quad,i_side))*temp_twist
          
       end do
    end do


    !calculate the new energy
    if (in_move%N .eq. 1) then

       !low side crossing - #1
       temp_twist=0

       temp_coords(:,1)=in_ring%x(:,quad(1,1,1))
       temp_coords(:,3)=in_ring%x(:,quad(3,1,1))
       temp_coords(:,4)=in_ring%x(:,quad(4,1,1))
       temp_coords(:,2)=in_move%x(:,1)

       !calculate the array of displacement vectors
       !print *,"triplet"
       do i=1,3
          dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
          !print *,dx_triplet(:,i)
       end do

       !calculate the energy new energy based on the displacement vectors
       call dihedral(dx_triplet,p_tens,temp_twist)
       !print *,"temp_twist_new=",temp_twist

       deltaE=deltaE+(in_energy_params%twist_params(1)+&
            temp_params(1,1))*temp_twist   

       !low side crossing - #2
       temp_twist=0

       temp_coords(:,1)=in_ring%x(:,quad(1,2,1))
       temp_coords(:,2)=in_ring%x(:,quad(2,2,1))
       temp_coords(:,4)=in_ring%x(:,quad(4,2,1))
       temp_coords(:,3)=in_move%x(:,1)

       !calculate the array of displacement vectors
       !print *,"triplet"
       do i=1,3
          dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
          !print *,dx_triplet(:,i)
       end do

       !calculate the energy new energy based on the displacement vectors
       call dihedral(dx_triplet,p_tens,temp_twist)
       !print *,"temp_twist_new=",temp_twist

       deltaE=deltaE+(in_energy_params%twist_params(1)+&
            temp_params(2,1))*temp_twist
       
       !low side no crossing
       do i_quad=3,3

          temp_twist=0

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,i)=in_ring%x(:,quad(i,i_quad,1))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=lo+(i-i_quad-1)
             if (in_move%node_indices(in_move%N) .eq. in_ring%N_nodes) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             temp_coords(:,i)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,1))*temp_twist          

       end do

       !high side no crossing
       do i_quad=3,3

          temp_twist=0
          !print *,"i_quad=",i_quad

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,4-i+1)=in_ring%x(:,quad(4-i+1,i_quad,2))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=hi-(i-i_quad-1)
             !print *,"i=",i
             !print *,"temp_pos=",temp_pos
             if (in_move%node_indices(1) .eq. 1) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             !print *,"  temp_pos=",temp_pos
             temp_coords(:,4-i+1)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,2))*temp_twist            

       end do

       
    elseif (in_move%N .eq. 2) then

       !low side crossing
       temp_twist=0

       temp_coords(:,1)=in_ring%x(:,quad(1,1,1))
       temp_coords(:,4)=in_ring%x(:,quad(4,1,1))
       temp_coords(:,2:3)=in_move%x(:,:)

       !calculate the array of displacement vectors
       !print *,"triplet"
       do i=1,3
          dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
          !print *,dx_triplet(:,i)
       end do

       !calculate the energy new energy based on the displacement vectors
       call dihedral(dx_triplet,p_tens,temp_twist)
       !print *,"temp_twist_new=",temp_twist

       deltaE=deltaE+(in_energy_params%twist_params(1)+&
            temp_params(1,1))*temp_twist

       !low side no crossing
       do i_quad=2,3

          temp_twist=0

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,i)=in_ring%x(:,quad(i,i_quad,1))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=lo+(i-i_quad-1)
             if (in_move%node_indices(in_move%N) .eq. in_ring%N_nodes) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             temp_coords(:,i)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,1))*temp_twist

       end do

       !high side no crossing
       do i_quad=2,3

          temp_twist=0
          !print *,"i_quad=",i_quad

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,4-i+1)=in_ring%x(:,quad(4-i+1,i_quad,2))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=hi-(i-i_quad-1)
             !print *,"i=",i
             !print *,"temp_pos=",temp_pos
             if (in_move%node_indices(1) .eq. 1) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             !print *,"  temp_pos=",temp_pos
             temp_coords(:,4-i+1)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,2))*temp_twist

       end do

       
    else

       !low side
       do i_quad=1,3
          
          temp_twist=0

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,i)=in_ring%x(:,quad(i,i_quad,1))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=lo+(i-i_quad-1)
             if (in_move%node_indices(in_move%N) .eq. in_ring%N_nodes) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             temp_coords(:,i)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,1))*temp_twist
          
       end do

       !high side
       do i_quad=1,3

          temp_twist=0
          !print *,"i_quad=",i_quad

          !add monomers from the non-moving section
          do i=1,i_quad
             temp_coords(:,4-i+1)=in_ring%x(:,quad(4-i+1,i_quad,2))
          end do

          !add monomers from the moving section
          do i=i_quad+1,4
             temp_pos=hi-(i-i_quad-1)
             !print *,"i=",i
             !print *,"temp_pos=",temp_pos
             if (in_move%node_indices(1) .eq. 1) then
                if (temp_pos .lt. 1) then
                   temp_pos=temp_pos+in_move%N
                elseif (temp_pos .gt. in_move%N) then
                   temp_pos=temp_pos-in_move%N
                end if
             end if
             !print *,"  temp_pos=",temp_pos
             temp_coords(:,4-i+1)=in_move%x(:,temp_pos)
          end do

          !calculate the array of displacement vectors
          do i=1,3
             dx_triplet(:,i)=temp_coords(:,i+1)-temp_coords(:,i)
             !print *,dx_triplet(:,i)
          end do

          !calculate the energy new energy based on the displacement vectors
          call dihedral(dx_triplet,p_tens,temp_twist)
          !print *,"temp_twist_new=",temp_twist

          deltaE=deltaE+(in_energy_params%twist_params(1)+&
               temp_params(i_quad,2))*temp_twist

       end do
       
    end if

    !print *,in_energy_params%twist_params(1)
    !print *,"deltaE_total=",deltaE

  end subroutine calc_twist_deltaE


  !function to test the boundaries
  pure subroutine test_bounds_ring(bounds_flag,in_ring)

    implicit none

    logical, intent(out) :: bounds_flag

    type(ring), intent(in) :: in_ring

    integer :: i_node, i_bound
    real(rp) :: r_sqrd, lower_lim, upper_lim

    bounds_flag=.false.

    !test bounds for box in Cartesian space
    if (in_ring%bound_type .eq. 1) then
       do i_node=1,in_ring%N_nodes
          do i_bound=1,in_ring%N_bounds
             if ((in_ring%x(in_ring%bounds(0,i_bound),i_node) .lt.&
                  in_ring%bounds(1,i_bound)) .or.&
                  (in_ring%x(in_ring%bounds(0,i_bound),i_node) .gt.&
                  in_ring%bounds(2,i_bound))) then
                bounds_flag=.true.
                return
                !exit
             end if
          end do
          ! if (bounds_flag .eqv. .true.) then
          !    exit
          ! end if
       end do
       !test bounds for the interior of a sphere
    elseif (in_ring%bound_type .eq. 2) then
       r_sqrd=in_ring%bounds(0,1)**2.0d0
       do i_node=1,in_ring%N_nodes
          if (floor(sum((in_ring%x(:,i_node)-in_ring%bounds(1:3,1))**2.0d0)) .gt.&
               r_sqrd) then
             bounds_flag=.true.
             return
             !exit
          end if
       end do
    elseif (in_ring%bound_type .eq. 3) then
       r_sqrd=in_ring%bounds(-1,1)**2.0d0
       lower_lim=in_ring%bounds(3,1)-in_ring%bounds(0,0)/2.0d0
       upper_lim=in_ring%bounds(3,1)+in_ring%bounds(0,0)/2.0d0
       do i_node=1,in_ring%N_nodes
          if ((floor(sum((in_ring%x(1:2,i_node)-in_ring%bounds(1:2,1))**2.0d0)) .gt.&
               r_sqrd) .or.&
               (in_ring%x(3,i_node) .le. lower_lim) .or.&
               (in_ring%x(3,i_node) .ge. upper_lim)) then
             bounds_flag=.true.
             return
             !exit
          end if
       end do
    end if

  end subroutine test_bounds_ring

  !subroutine to randomly place the initial ring
  subroutine rand_init_ring(in_ring,in_obstacle_distribution)

    implicit none

    type(ring), intent(inout) :: in_ring
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    logical :: valid_location, bounds_flag
    
    integer :: max_iter, count
    integer :: obst_test(1:2)
    integer(ip) :: x0(1:in_ring%dim,in_ring%N_nodes)

    valid_location=.false.
    max_iter=1000
    count=0

    x0=in_ring%x

    !write(unit_log,*)'called rand_init_ring'

    do while ((valid_location .eqv. .false.) .and.&
         (count .lt. max_iter))

       !write(unit_log,*)'entered loop'
       
       call propose_shift(in_ring,x0)

       !write(unit_log,*)in_ring%x

       bounds_flag=.true.

       call test_bounds_ring(bounds_flag,in_ring)

       if (bounds_flag .eqv. .true.) then
          count=count+1
          cycle
       end if

       obst_test = obstacle_test_fast(in_ring,in_obstacle_distribution)

       if ((sum(obst_test) .eq. -2) .and.&
            (bounds_flag .eqv. .false.)) then
          valid_location=.true.
       else
          count=count+1
          cycle
       end if
       
    end do

    if (valid_location .eqv. .false.) then
       !write(unit_log,*)'no new location'
       in_ring%x=x0
    end if
    
  end subroutine rand_init_ring

end module ring_config
