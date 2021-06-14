module ring_polymer

  !this module contains the subroutines and functions used to generate
  !configurations of the polymer ring
  
  use params
  use ring_objects
  use read_write
  use ring_procedures
  use searches_and_sorts
  use ring_config
  use energy_terms

  implicit none

contains
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! KINK MOVE ROUTINES AND FUNCTIONS !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !function to count the number of kink moves
  function count_kink_moves(in_ring)

    implicit none
    
    type(ring), intent(in) :: in_ring
    type(node_subsets) :: count_kink_moves
    integer :: i_node, i_fixed, prev_node_index
    integer :: alignment(1:in_ring%N_nodes)
    integer, allocatable :: new_entry(:)

    allocate(new_entry(1:3))
    new_entry=0
    alignment=0

    count_kink_moves=new_node_subsets(0,3)

    alignment=sum((cshift(in_ring%x,1,2)-in_ring%x)*&
         (in_ring%x-cshift(in_ring%x,-1,2)),1)
    
    !if there are no fixed nodes
    if (in_ring%fixed .eqv. .false.) then
       
       !loop over all nodes
       do i_node=1,in_ring%N_nodes

          if (alignment(i_node) .eq. 0) then
             
             if (i_node .eq. 1) then
                new_entry=(/in_ring%N_nodes,1,2/)
             elseif (i_node .eq. in_ring%N_nodes) then
                new_entry=(/i_node-1,i_node,1/)
             else
                new_entry=(/i_node-1,i_node,i_node+1/)
             end if

             call add_to_subset_list(count_kink_moves%nodes,&
                  count_kink_moves%N,count_kink_moves%N_max,&
                  count_kink_moves%M,new_entry)
             
          end if
          
       end do !end loop over all nodes

    !if there are fixed nodes
    else

       do i_fixed=1,in_ring%N_fixed

          if (i_fixed .eq. 1) then !at first fixed node

             !option 1
             if (in_ring%fixed_nodes(i_fixed) .gt. 1) then

                do i_node=1,in_ring%fixed_nodes(i_fixed)-1

                   if (alignment(i_node) .eq. 0) then

                      if (i_node .eq. 1) then
                         new_entry=(/in_ring%N_nodes,1,2/)
                      else
                         new_entry=(/i_node-1,i_node,i_node+1/)
                      end if

                      call add_to_subset_list(count_kink_moves%nodes,&
                           count_kink_moves%N,count_kink_moves%N_max,&
                           count_kink_moves%M,new_entry)

                   end if

                end do

             end if

             !option 2
             if ((in_ring%N_fixed .eq. 1) .and.&
                  (in_ring%fixed_nodes(i_fixed) .lt. in_ring%N_nodes)) then

                do i_node=in_ring%fixed_nodes(i_fixed)+1,in_ring%N_nodes

                   if (alignment(i_node) .eq. 0) then

                      if (i_node .eq. in_ring%N_nodes) then
                         new_entry=(/i_node-1,i_node,1/)
                      else
                         new_entry=(/i_node-1,i_node,i_node+1/)
                      end if

                      call add_to_subset_list(count_kink_moves%nodes,&
                           count_kink_moves%N,count_kink_moves%N_max,&
                           count_kink_moves%M,new_entry)

                   end if
                   
                end do

             end if

          elseif (i_fixed .eq. in_ring%N_fixed) then !at last fixed node

             !option 1
             if (in_ring%fixed_nodes(i_fixed)-prev_node_index .gt. 1) then

                do i_node=prev_node_index+1,in_ring%fixed_nodes(i_fixed)-1

                   if (alignment(i_node) .eq. 0) then

                      new_entry=(/i_node-1,i_node,i_node+1/)

                      call add_to_subset_list(count_kink_moves%nodes,&
                           count_kink_moves%N,count_kink_moves%N_max,&
                           count_kink_moves%M,new_entry)

                   end if
                   
                end do
                
             end if

             !option 2
             if (in_ring%fixed_nodes(i_fixed) .lt. in_ring%N_nodes) then

                do i_node=in_ring%fixed_nodes(i_fixed)+1,in_ring%N_nodes

                   if (alignment(i_node) .eq. 0) then

                      if (i_node .eq. in_ring%N_nodes) then
                         new_entry=(/i_node-1,i_node,1/)
                      else
                         new_entry=(/i_node-1,i_node,i_node+1/)
                      end if

                      call add_to_subset_list(count_kink_moves%nodes,&
                           count_kink_moves%N,count_kink_moves%N_max,&
                           count_kink_moves%M,new_entry)

                   end if

                end do

             end if
             
          else !between two fixed nodes

             if (in_ring%fixed_nodes(i_fixed)-prev_node_index .gt. 1) then

                do i_node=prev_node_index+1,in_ring%fixed_nodes(i_fixed)+1
                   
                   if (alignment(i_node) .eq. 0) then

                      new_entry=(/i_node-1,i_node,i_node+1/)

                      call add_to_subset_list(count_kink_moves%nodes,&
                           count_kink_moves%N,count_kink_moves%N_max,&
                           count_kink_moves%M,new_entry)

                   end if
                   
                end do
                
             end if
             
          end if !end conditionals texting interval

          prev_node_index=in_ring%fixed_nodes(i_fixed)
          
       end do !end loop over fixed nodes

    end if !end conditional for fixed nodes

    if (allocated(new_entry)) deallocate(new_entry)
    
  end function count_kink_moves
  

  !function to generate a kink move ring subset
  function kink_move(in_ring,triplet)

    implicit none
    type(ring), intent(in) :: in_ring
    integer, intent(in) :: triplet(1:3)
    type(ring_subset) :: kink_move

    integer :: i_dim, mismatched_dims(1:2), count

    kink_move=new_ring_subset(in_ring%dim,1)
    kink_move%node_indices(1)=triplet(2)
    kink_move%x(:,1)=in_ring%x(:,triplet(2))

    !write(unit_log,*)"kink_move_triplet=",triplet


    !if 2D swap the two dimensions in the middle node
    if (in_ring%dim .eq. 2) then
       
       do i_dim=1,2

          if (in_ring%x(i_dim,triplet(1)) .eq. in_ring%x(i_dim,triplet(2))) then
             kink_move%x(i_dim,1)=in_ring%x(i_dim,triplet(3))
          else
             kink_move%x(i_dim,1)=in_ring%x(i_dim,triplet(1))
          end if

       end do

    !if 3D find the mismatched dimensions
    elseif (in_ring%dim .eq. 3) then

       count=1

       do i_dim=1,3
          if (in_ring%x(i_dim,triplet(1)) .ne. in_ring%x(i_dim,triplet(3))) then
             mismatched_dims(count)=i_dim
             count=count+1
          end if
       end do

       !write(unit_log,*)"kink_move_mismatch=",mismatched_dims

       !swap the mismatched dimensions in the middle node
       do i_dim=1,2

          if (in_ring%x(mismatched_dims(i_dim),triplet(1)) .eq.&
               in_ring%x(mismatched_dims(i_dim),triplet(2))) then
             kink_move%x(mismatched_dims(i_dim),1)=&
                  in_ring%x(mismatched_dims(i_dim),triplet(3))
          else
             kink_move%x(mismatched_dims(i_dim),1)=&
                  in_ring%x(mismatched_dims(i_dim),triplet(1))
          end if

       end do
       
    end if
    
  end function kink_move
  
  !subroutine to repeat kink moves in prior direction
  subroutine repeat_kink_moves(fin_node,dir,&
       in_ring,node_moved,max_dist)

    implicit none

    integer, intent(out) :: fin_node, dir

    type(ring), intent(inout) :: in_ring
    
    integer, intent(in) :: node_moved, max_dist
    
    type(ring_subset) :: new_kink_move
    
    logical :: valid_fwd, valid_bwd, move_success

    integer :: i_node, j_node, alignment, curr_node, dist
    integer :: triplet(1:3)
    integer(ip) :: dx(1:in_ring%dim,1:2)

    valid_fwd=.false.
    valid_bwd=.false.

    fin_node=-1

    
    ! test if we can repeat kink moves in the backwards direction
    j_node=count_transform_alt(-2,node_moved,in_ring%N_nodes)
    i_node=count_transform_alt(-1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=node_moved

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_bwd=.true.
    end if

    !test if we can repeat kink_moves in the forwards direction
    j_node=node_moved
    i_node=count_transform_alt(1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=count_transform_alt(2,node_moved,in_ring%N_nodes)

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_fwd=.true.
    end if

    !choose forwards or backwards moves
    if ((valid_fwd .eqv. .true.) .and.&
         (valid_bwd .eqv. .true.)) then
       !write(unit_log,*)"both dir valid"
       if (int_rand(0,1) .eq. 1) then
          valid_fwd=.false.
          !write(unit_log,*)"choose backward"
       else
          valid_bwd=.false.
          !write(unit_log,*)"choose forward"
       end if
    end if

    !repeat kink moves in the backward direction
    if (valid_bwd .eqv. .true.) then

       dir=-1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"bwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move(move_success,in_ring,new_kink_move)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if
          
       end do
       
    end if

    !repeat kink moves in the forward direction
    if (valid_fwd .eqv. .true.) then

       dir=1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"fwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move(move_success,in_ring,new_kink_move)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if

       end do

    end if
    
  end subroutine repeat_kink_moves
  
  !subroutine to repeat kink moves in the presence of obstacles
  subroutine repeat_kink_moves_obstacles(fin_node,dir,&
       in_ring,node_moved,max_dist,&
       in_obstacle_distribution)

    implicit none

    integer, intent(out) :: fin_node, dir

    type(ring), intent(inout) :: in_ring
    
    integer, intent(in) :: node_moved, max_dist
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution

    type(ring_subset) :: new_kink_move
    
    logical :: valid_fwd, valid_bwd, move_success

    integer :: i_node, j_node, alignment, curr_node, dist
    integer :: triplet(1:3)
    integer(ip) :: dx(1:in_ring%dim,1:2)

    valid_fwd=.false.
    valid_bwd=.false.

    fin_node=-1

    
    ! test if we can repeat kink moves in the backwards direction
    j_node=count_transform_alt(-2,node_moved,in_ring%N_nodes)
    i_node=count_transform_alt(-1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=node_moved

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_bwd=.true.
    end if

    !test if we can repeat kink_moves in the forwards direction
    j_node=node_moved
    i_node=count_transform_alt(1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=count_transform_alt(2,node_moved,in_ring%N_nodes)

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_fwd=.true.
    end if

    !choose forwards or backwards moves
    if ((valid_fwd .eqv. .true.) .and.&
         (valid_bwd .eqv. .true.)) then
       !write(unit_log,*)"both dir valid"
       if (int_rand(0,1) .eq. 1) then
          valid_fwd=.false.
          !write(unit_log,*)"choose backward"
       else
          valid_bwd=.false.
          !write(unit_log,*)"choose forward"
       end if
    end if

    !repeat kink moves in the backward direction
    if (valid_bwd .eqv. .true.) then

       dir=-1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"bwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_obstacles(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if
          
       end do
       
    end if

    !repeat kink moves in the forward direction
    if (valid_fwd .eqv. .true.) then

       dir=1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"fwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_obstacles(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if

       end do

    end if
    
  end subroutine repeat_kink_moves_obstacles

  !subroutine to repeat kink moves with MH
  subroutine repeat_kink_moves_MH(fin_node,dir,&
       in_ring,node_moved,max_dist,&
       in_energy_params)

    implicit none

    integer, intent(out) :: fin_node, dir

    type(ring), intent(inout) :: in_ring
    
    integer, intent(in) :: node_moved, max_dist
    type(energy_params), intent(in) :: in_energy_params

    type(ring_subset) :: new_kink_move
    
    logical :: valid_fwd, valid_bwd, move_success

    integer :: i_node, j_node, alignment, curr_node, dist
    integer :: triplet(1:3)
    integer(ip) :: dx(1:in_ring%dim,1:2)

    valid_fwd=.false.
    valid_bwd=.false.

    fin_node=-1

    
    ! test if we can repeat kink moves in the backwards direction
    j_node=count_transform_alt(-2,node_moved,in_ring%N_nodes)
    i_node=count_transform_alt(-1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=node_moved

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_bwd=.true.
    end if

    !test if we can repeat kink_moves in the forwards direction
    j_node=node_moved
    i_node=count_transform_alt(1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=count_transform_alt(2,node_moved,in_ring%N_nodes)

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_fwd=.true.
    end if

    !choose forwards or backwards moves
    if ((valid_fwd .eqv. .true.) .and.&
         (valid_bwd .eqv. .true.)) then
       !write(unit_log,*)"both dir valid"
       if (int_rand(0,1) .eq. 1) then
          valid_fwd=.false.
          !write(unit_log,*)"choose backward"
       else
          valid_bwd=.false.
          !write(unit_log,*)"choose forward"
       end if
    end if

    !repeat kink moves in the backward direction
    if (valid_bwd .eqv. .true.) then

       dir=-1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"bwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_MH(move_success,in_ring,new_kink_move,&
                  in_energy_params)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if
          
       end do
       
    end if

    !repeat kink moves in the forward direction
    if (valid_fwd .eqv. .true.) then

       dir=1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"fwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_MH(move_success,in_ring,new_kink_move,&
                  in_energy_params)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if

       end do

    end if
    
  end subroutine repeat_kink_moves_MH

  !subroutine to repeat kink moves with MH and in the presence of obstacles
  subroutine repeat_kink_moves_obstacles_MH(fin_node,dir,&
       in_ring,node_moved,max_dist,&
       in_obstacle_distribution,in_energy_params)

    implicit none

    integer, intent(out) :: fin_node, dir

    type(ring), intent(inout) :: in_ring
    
    integer, intent(in) :: node_moved, max_dist
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params

    type(ring_subset) :: new_kink_move
    
    logical :: valid_fwd, valid_bwd, move_success

    integer :: i_node, j_node, alignment, curr_node, dist
    integer :: triplet(1:3)
    integer(ip) :: dx(1:in_ring%dim,1:2)

    valid_fwd=.false.
    valid_bwd=.false.

    fin_node=-1

    
    ! test if we can repeat kink moves in the backwards direction
    j_node=count_transform_alt(-2,node_moved,in_ring%N_nodes)
    i_node=count_transform_alt(-1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=node_moved

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_bwd=.true.
    end if

    !test if we can repeat kink_moves in the forwards direction
    j_node=node_moved
    i_node=count_transform_alt(1,node_moved,in_ring%N_nodes)

    dx(:,1)=in_ring%x(:,i_node)-in_ring%x(:,j_node)
    
    j_node=count_transform_alt(2,node_moved,in_ring%N_nodes)

    dx(:,2)=in_ring%x(:,j_node)-in_ring%x(:,i_node)

    alignment=sum(dx(:,1)*dx(:,2))

    if (alignment .eq. 0) then
       valid_fwd=.true.
    end if

    !choose forwards or backwards moves
    if ((valid_fwd .eqv. .true.) .and.&
         (valid_bwd .eqv. .true.)) then
       !write(unit_log,*)"both dir valid"
       if (int_rand(0,1) .eq. 1) then
          valid_fwd=.false.
          !write(unit_log,*)"choose backward"
       else
          valid_bwd=.false.
          !write(unit_log,*)"choose forward"
       end if
    end if

    !repeat kink moves in the backward direction
    if (valid_bwd .eqv. .true.) then

       dir=-1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"bwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_obstacles_MH(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution,in_energy_params)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if
          
       end do
       
    end if

    !repeat kink moves in the forward direction
    if (valid_fwd .eqv. .true.) then

       dir=1

       dist=0
       move_success=.true.
       curr_node=count_transform_alt(dir,node_moved,in_ring%N_nodes)
       !write(unit_log,*)"fwd repeats:",node_moved,dir,curr_node

       do while ((move_success .eqv. .true.) .and.&
            (dist .lt. max_dist))

          j_node=count_transform_alt(-1,curr_node,in_ring%N_nodes)
          i_node=count_transform_alt(1,curr_node,in_ring%N_nodes)

          dx(:,1)=in_ring%x(:,curr_node)-in_ring%x(:,j_node)
          dx(:,2)=in_ring%x(:,i_node)-in_ring%x(:,curr_node)

          alignment=sum(dx(:,1)*dx(:,2))

          if (alignment .eq. 0) then
             triplet=(/j_node,curr_node,i_node/)
             !write(unit_log,*)triplet
             new_kink_move=kink_move(in_ring,triplet)
             call perform_move_obstacles_MH(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution,in_energy_params)
          else
             move_success=.false.
          end if

          if (move_success .eqv. .true.) then
             dist=dist+1
             fin_node=curr_node
             !write(unit_log,*)"fin_node=",fin_node
             curr_node=count_transform_alt(dir,curr_node,in_ring%N_nodes)
          end if

       end do

    end if
    
  end subroutine repeat_kink_moves_obstacles_MH
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! REFLECT MOVE ROUTINES AND FUNCTIONS !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !function to count the number of reflect moves
  function count_reflect_moves(in_ring)

    implicit none

    type(ring), intent(in) :: in_ring
    type(direction_pairs) :: count_reflect_moves
    type(direction_pairs) :: temp_reflect_moves, temp_reflect_moves_fixed

    logical :: prev_pair
    integer :: pair_temp(1:2), dist_dir(1:2)
    integer :: count, count_trans, matches
    integer :: i_node, i_dim


    logical :: prev_node_trigger, search_flag
    integer :: index_preorigin_cross, lo_1, hi_1, lo_2, hi_2, prev_node_N
    integer, allocatable :: temp_swap(:,:)
    
    temp_reflect_moves=new_direction_pairs(0,2)

    prev_node_trigger=.false.
    search_flag=.false.
    
    !iterate over the nodes in the ring
    do i_node=1,in_ring%N_nodes

       !initialize with no origin crossing
       index_preorigin_cross=-1

       !assign the index from the end of the previous node iterations
       if (prev_node_trigger .eqv. .true.) then
          prev_node_N=temp_reflect_moves%N
          search_flag=.true.
       end if

       !start the count by offsetting by 3 nodes, this is the minimum distance, for a node subset
       !to be shifted off the axis of the reflect move
       count=3

       !iterate to loop back around the ring
       !do while (count .lt. in_ring%N_nodes-2)
       do while (count .lt. min(in_ring%N_nodes/4,500))

          !transform the count moved along the ring to account for crossing the origin of the ring
          count_trans=count_transform(count,i_node,in_ring%N_nodes)

          !calculate the number of matches along each dimension from the node indexed by i_node to
          !the node indexed by the transformed count variable
          matches=0
          do i_dim=1,in_ring%dim
             if (in_ring%x(i_dim,i_node) .eq. in_ring%x(i_dim,count_trans)) then
                matches=matches+1
             end if
          end do

          !if the number of matched dimensions is equal to the number of dimensions minus one, then proceed
          if (matches .eq. in_ring%dim-1) then

             pair_temp(1:2)=(/i_node,count_trans/)

             !if the number of moves found is greater than zero, then test if the move is equivalent to a past move
             if (temp_reflect_moves%N .gt. 0) then
                
                !if i_node is past the first node and the search flag is triggered, then perform a search
                if ((i_node .ne. 1) .and. (search_flag .eqv. .true.)) then

                   !assume the pair is repeated
                   prev_pair=.true.

                   !search the along the first nodes in all of the previous pairs for the second node
                   !of the proposed reflect pair
                   call binary_search(lo_1,hi_1,temp_reflect_moves%nodes(1,1:prev_node_N),&
                        prev_node_N,pair_temp(2))
                   
                   !if there is an interval matching the second node of the proposed reflect pair, then proceed
                   if ((lo_1 .ne. -1) .and. (hi_1 .ne. -1)) then

                      !search the interval of the second nodes for the first node of the proposed reflect pair
                      call binary_search(lo_2,hi_2,temp_reflect_moves%nodes(2,lo_1:hi_1),&
                           hi_1-lo_1+1,pair_temp(1))

                      !if there is not a match to the second node of the proposed reflect pair, set the previous
                      !pair flag to false
                      if ((lo_2 .eq. -1) .and. (hi_2 .eq. -1)) then
                         prev_pair=.false.
                      end if
                      
                   else !if there is not a matching interval for the second node of the proposed reflect pair,
                      !set the previous pair flag to false
                      prev_pair=.false.
                   end if

                else !if the initial conditional was false, set the previous pair flag to false
                   prev_pair=.false.
                end if

                !if the proposed pair was not a previous reflect pair, then
                if (prev_pair .eqv. .false.) then

                   !calculate the direction pair
                   dist_dir=min_distance(in_ring%N_nodes,pair_temp(1),pair_temp(2))
                   
                   !test if origin has been crossed
                   if ((count_trans .lt. i_node) .and. (index_preorigin_cross .eq. -1)) then
                      index_preorigin_cross=temp_reflect_moves%N
                   end if

                   !add the reflect pair and the direction pair to the direction pairs list
                   call add_to_direction_pairs_list(temp_reflect_moves%nodes,&
                        temp_reflect_moves%dirs,temp_reflect_moves%N,&
                        temp_reflect_moves%N_max,pair_temp,dist_dir)
                   
                end if !end conditional testing for a previous pair
                
             else !if there are no direction pairs so far, start the direction pairs list

                !calculate the direction pair
                dist_dir=min_distance(in_ring%N_nodes,pair_temp(1),pair_temp(2))

                !initialize the direction pairs list with the first pair found
                temp_reflect_moves=new_direction_pairs(1,2)
                temp_reflect_moves%nodes(:,1)=pair_temp
                temp_reflect_moves%dirs(:,1)=dist_dir

                !trigger the previous node flag
                prev_node_trigger=.true.
                
             end if !end the conditional for number of direction pairs
             
          end if !end the conditional for the number of matched dimensions

          count=count+1 !increase the displacement along the ring
          
       end do !finish the iterator increasing the count

       !swap portions of the array here to maintain a sorted array following an origin crossing

       !Here is a brief explanation:
       !Suppose i_node = N_nodes - 10, we may find the direction pairs [N_nodes-10, N_nodes-5] and
       !pass around the ring, crossing the origin, to find another direction pair, say
       ![N_node-10, 5], if we add this directly to the direction pairs list, then the elements
       !along the second index of the pairs of nodes will be unsorted within each interval defined
       !by a single value in the first index of the pairs of nodes, this will cause issues for the
       !binary searches used to avoid redundant moves

       
       if ((index_preorigin_cross .ne. -1) .and. (index_preorigin_cross .gt. prev_node_N)) then
          allocate(temp_swap(1:2,1:temp_reflect_moves%N-index_preorigin_cross))
          !swap the nodes
          temp_swap=temp_reflect_moves%nodes(:,prev_node_N+1:index_preorigin_cross)
          temp_reflect_moves%nodes(:,prev_node_N+1:&
               temp_reflect_moves%N-index_preorigin_cross+prev_node_N)=&
               temp_reflect_moves%nodes(:,index_preorigin_cross+1:temp_reflect_moves%N)
          temp_reflect_moves%nodes(:,temp_reflect_moves%N-index_preorigin_cross+prev_node_N+1:&
               temp_reflect_moves%N)=temp_swap
          !swap the dirs
          temp_swap=temp_reflect_moves%dirs(:,prev_node_N+1:index_preorigin_cross)
          temp_reflect_moves%dirs(:,prev_node_N+1:&
               temp_reflect_moves%N-index_preorigin_cross+prev_node_N)=&
               temp_reflect_moves%dirs(:,index_preorigin_cross+1:temp_reflect_moves%N)
          temp_reflect_moves%dirs(:,temp_reflect_moves%N-index_preorigin_cross+prev_node_N+1:&
               temp_reflect_moves%N)=temp_swap
          if (allocated(temp_swap)) deallocate(temp_swap)
       end if

       !if the first i_node with a direction pair has been triggered, reduce the reflect pairs

       !the motivation for this is to reduce the memory overhead
       !edit by BRG - removing this may speed up the enumeration of reflect moves
       
       if ((prev_node_trigger .eqv. .true.) .and. (favor_memory .eqv. .true.)) then
          temp_reflect_moves=red_reflect_moves(in_ring,temp_reflect_moves)
       end if

    end do !end loop over nodes
    
    !remove contiguous lines and test for fixed nodes
    if (in_ring%fixed .eqv. .false.) then

       !if there are no fixed nodes just remove the contiguous lines
       count_reflect_moves=red_reflect_moves(in_ring,temp_reflect_moves)
       
    else

       !remove the reflect pairs with a fixed node between them
       temp_reflect_moves_fixed=red_reflect_moves_fixed(in_ring,temp_reflect_moves)

       !if the remaining number of reflect pairs is greater than 0, then
       !remove the contiguous lines
       if (temp_reflect_moves_fixed%N .gt. 0) then

          count_reflect_moves=red_reflect_moves(in_ring,temp_reflect_moves_fixed)

       end if
              
    end if

    if (allocated(temp_reflect_moves%nodes)) deallocate(temp_reflect_moves%nodes)
    if (allocated(temp_reflect_moves%dirs)) deallocate(temp_reflect_moves%dirs)
    if (allocated(temp_reflect_moves_fixed%nodes)) deallocate(temp_reflect_moves_fixed%nodes)
    if (allocated(temp_reflect_moves_fixed%dirs)) deallocate(temp_reflect_moves_fixed%dirs)
    
  end function count_reflect_moves

  !function to reduce the reflect moves by removing redundant reflect moves
  function red_reflect_moves(in_ring,in_reflect_moves)

    implicit none

    type(ring), intent(in) :: in_ring
    type(direction_pairs), intent(in) :: in_reflect_moves
    type(direction_pairs) :: red_reflect_moves

    integer :: pair_temp(1:2), dist_dir_old(1:2), dist_dir_new(1:2)
    integer :: j_new, k_new, j_old, k_old, k_diff
    integer :: i_pair

    k_diff=1
    
    !iterate over the list of pairs
    do i_pair=1,in_reflect_moves%N

       !if starting with the first pair, add it
       if (i_pair .eq. 1) then

          red_reflect_moves=new_direction_pairs(1,2)
          red_reflect_moves%nodes(:,1)=in_reflect_moves%nodes(:,1)
          red_reflect_moves%dirs(:,1)=in_reflect_moves%dirs(:,1)
          k_diff=1

       else !if not the first pair, proceed with testing for redundancy

          j_new=in_reflect_moves%nodes(1,i_pair) !first node of new pair
          k_new=in_reflect_moves%nodes(2,i_pair) !second node of new pair
          pair_temp=(/j_new,k_new/)

          dist_dir_new=in_reflect_moves%dirs(:,i_pair) !distance of the new pair

          j_old=red_reflect_moves%nodes(1,red_reflect_moves%N)  !first node of the old pair
          k_old=red_reflect_moves%nodes(2,red_reflect_moves%N) !second node of the old pair

          !test if the first nodes in the pair are the same
          if (j_new .eq. j_old) then

             if (k_old .eq. k_new-k_diff) then !if the second node of the new pair is directly beside the second
                ! node of the old pair, then proceed

                dist_dir_old=red_reflect_moves%dirs(:,red_reflect_moves%N) !distance of the old pair

                !if the new distance is less than the old distance, then replace the old with the new reflect pair,
                !this will reduce the number of nodes whose new positions will need to be calculated by the rotation matrix
                if (dist_dir_new(1) .lt. dist_dir_old(1)) then
                   red_reflect_moves%nodes(:,red_reflect_moves%N)=pair_temp
                   red_reflect_moves%dirs(:,red_reflect_moves%N)=dist_dir_new
                else !if the new distance is not less than the old distance, increase the second node difference
                   k_diff=k_diff+1
                end if

             else ! if the second node of the new pair is not directly beside the second node of the old pair,
                ! then add the new reflect pair to the reduced list and reset the second node difference
                k_diff=1
                call add_to_direction_pairs_list(red_reflect_moves%nodes,&
                     red_reflect_moves%dirs,red_reflect_moves%N,&
                     red_reflect_moves%N_max,pair_temp,dist_dir_new)
             end if

          else !if the first nodes are not the same, then add the new reflect pair to the reduced list and
             ! reset the second node difference
             k_diff=1
             call add_to_direction_pairs_list(red_reflect_moves%nodes,&
                  red_reflect_moves%dirs,red_reflect_moves%N,&
                  red_reflect_moves%N_max,pair_temp,dist_dir_new)
          end if

       end if !end conditional testing for first pair

    end do !finish iterating over the reflect pairs

  end function red_reflect_moves

  !function to reduce the reflect moves by removing the reflect moves containing a fixed node
  function red_reflect_moves_fixed(in_ring,in_reflect_moves)

    implicit none

    type(ring), intent(in) :: in_ring
    type(direction_pairs), intent(in) :: in_reflect_moves
    type(direction_pairs) :: red_reflect_moves_fixed

    integer :: temp_dist_dir(1:2)
    integer :: i_pair
    logical :: between_flag, initial_pair

    initial_pair=.false.

    !iterate over the list of pairs
    do i_pair=1,in_reflect_moves%N

       !test if the initial pair has been found
       if (initial_pair .eqv. .false.) then

          between_flag=.false.

          !test fixed between for default direction
          between_flag=fixed_between(in_ring%N_fixed,in_ring%fixed_nodes,&
               in_reflect_moves%nodes(:,i_pair),in_reflect_moves%dirs(:,i_pair))

          !set initial pair if fixed is not between
          if (between_flag .eqv. .false.) then
             
             red_reflect_moves_fixed=new_direction_pairs(1,2)
             red_reflect_moves_fixed%nodes(:,1)=in_reflect_moves%nodes(:,i_pair)
             red_reflect_moves_fixed%dirs(:,1)=in_reflect_moves%dirs(:,i_pair)

             initial_pair=.true.
             
          else !if a fixed node was between, test the reversed direction

             between_flag=.false.

             !create a temporary distance with the reversed direction
             temp_dist_dir=(/in_ring%N_nodes-in_reflect_moves%dirs(1,i_pair),&
                  -1*in_reflect_moves%dirs(2,i_pair)/)

             !test fixed between for reversed direction
             between_flag=fixed_between(in_ring%N_fixed,in_ring%fixed_nodes,&
                  in_reflect_moves%nodes(:,i_pair),temp_dist_dir)

             !set initial pair as the reversed if fixed is not between
             if (between_flag .eqv. .false.) then
                
                red_reflect_moves_fixed=new_direction_pairs(1,2)
                red_reflect_moves_fixed%nodes(:,1)=in_reflect_moves%nodes(:,i_pair)
                red_reflect_moves_fixed%dirs(:,1)=temp_dist_dir

                initial_pair=.true.

             end if
             
          end if !end conditional swapping direction

       else !if the initial pair was found, test the later pairs

          between_flag=.false.

          !test fixed between for default direction with the new candidate pair
          between_flag=fixed_between(in_ring%N_fixed,in_ring%fixed_nodes,&
               in_reflect_moves%nodes(:,i_pair),in_reflect_moves%dirs(:,i_pair))

          !if the new candidate pair does not have a fixed node between, add it
          if (between_flag .eqv. .false.) then

             !add the new pair to the list
             call add_to_direction_pairs_list(red_reflect_moves_fixed%nodes,&
                  red_reflect_moves_fixed%dirs,red_reflect_moves_fixed%N,&
                  red_reflect_moves_fixed%N_max,in_reflect_moves%nodes(:,i_pair),&
                  in_reflect_moves%dirs(:,i_pair))

          else !if a fixed node was between, test the reversed direction

             between_flag=.false.

             !create a temporary distance with the reversed direction
             temp_dist_dir=(/in_ring%N_nodes-in_reflect_moves%dirs(1,i_pair),&
                  -1*in_reflect_moves%dirs(2,i_pair)/)

             !test fixed between for reversed direction
             between_flag=fixed_between(in_ring%N_fixed,in_ring%fixed_nodes,&
                  in_reflect_moves%nodes(:,i_pair),temp_dist_dir)

             !if the new reversed candidate pair does not have a fixed node between, add it
             if (between_flag .eqv. .false.) then

                !add the reversed pair to the list
                call add_to_direction_pairs_list(red_reflect_moves_fixed%nodes,&
                     red_reflect_moves_fixed%dirs,red_reflect_moves_fixed%N,&
                     red_reflect_moves_fixed%N_max,in_reflect_moves%nodes(:,i_pair),&
                     temp_dist_dir)

             end if

          end if !end conditional swapping the direction
          
       end if !end conditional for the initial pair being found

    end do !finish iterating over the set of reflect pairs

    !check if there were no reflect pairs without fixed nodes in between
    if (initial_pair .eqv. .false.) then
       red_reflect_moves_fixed=new_direction_pairs(1,2)
       red_reflect_moves_fixed%N=0
    end if

  end function red_reflect_moves_fixed

  !function to generate a reflect move ring subset
  function reflect_move(in_ring,pair,dist_dir)

    implicit none
    
    type(ring), intent(in) :: in_ring
    integer, intent(in) :: pair(1:2), dist_dir(1:2)
    type(ring_subset) :: reflect_move

    integer :: i_node, count
    real(rp) :: rot_matrix(1:in_ring%dim,in_ring%dim), x_mid(1:in_ring%dim)
    real(rp) :: x_temp(1:in_ring%dim,1:dist_dir(1)-1)


    reflect_move=new_ring_subset(in_ring%dim,dist_dir(1)-1)

    !find the midpoint of the pair and generate a rotation matrix
    x_mid=(in_ring%x(:,pair(1))+in_ring%x(:,pair(2)))/2.0d0
    rot_matrix=rot_mat(in_ring%dim,in_ring%x(:,pair(1)),in_ring%x(:,pair(2)),1)

    !increasing orientation of pair
    if (pair(2) .gt. pair(1)) then
       !forward direction
       if (dist_dir(2) .eq. 1) then
          count=1
          do i_node=pair(1)+1,pair(2)-1
             reflect_move%node_indices(count)=i_node
             x_temp(:,count)=in_ring%x(:,i_node)-x_mid
             count=count+1
          end do
       !reverse direction
       else
          count=1
          if (pair(1) .gt. 1) then
             do i_node=1,pair(1)-1
                reflect_move%node_indices(count)=i_node
                x_temp(:,count)=in_ring%x(:,i_node)-x_mid
                count=count+1
             end do
          end if
          if (pair(2) .lt. in_ring%N_nodes) then
             do i_node=pair(2)+1,in_ring%N_nodes
                reflect_move%node_indices(count)=i_node
                x_temp(:,count)=in_ring%x(:,i_node)-x_mid
                count=count+1
             end do
          end if
       end if
    !decreasing orientation of pair
    else
       !forward direction
       if (dist_dir(2) .eq. 1) then
          count=1
          if (pair(2) .gt. 1) then
             do i_node=1,pair(2)-1
                reflect_move%node_indices(count)=i_node
                x_temp(:,count)=in_ring%x(:,i_node)-x_mid
                count=count+1
             end do
          end if
          if (pair(1) .lt. in_ring%N_nodes) then
             do i_node=pair(1)+1,in_ring%N_nodes
                reflect_move%node_indices(count)=i_node
                x_temp(:,count)=in_ring%x(:,i_node)-x_mid
                count=count+1
             end do
          end if
       !reverse direction
       else
          count=1
          do i_node=pair(2)+1,pair(1)-1
             reflect_move%node_indices(count)=i_node
             x_temp(:,count)=in_ring%x(:,i_node)-x_mid
             count=count+1
          end do
       end if
    end if
    
    !apply the rotational matrix to the subset of nodes between the reflect pair
    x_temp=matmul(rot_matrix,x_temp)
    ! call dgemm('n','n',&
    !      in_ring%dim,reflect_move%N,in_ring%dim,&
    !      1.0d0,&
    !      rot_matrix,in_ring%dim,&
    !      dble(x_temp),in_ring%dim,&
    !      0.0d0,&
    !      x_temp,in_ring%dim)
    
    do i_node=1,reflect_move%N
       x_temp(:,i_node)=x_mid+x_temp(:,i_node)
    end do
    
    !round the result and convert it to integers
    reflect_move%x=int(nint(x_temp))
    
  end function reflect_move
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! OBSTACLE ROUTINES AND FUNCTIONS !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !subroutine to attempt to place an obstacle
  subroutine place_obstacle(success,in_obstacle_distribution,in_ring,in_deform_params)

    implicit none

    logical, intent(out) :: success

    type(obstacle_distribution), intent(inout) :: in_obstacle_distribution

    type(ring), intent(in) :: in_ring
    type(deform_params), intent(in) :: in_deform_params

    logical :: intersect_flag

    integer :: lo, hi
    integer :: site_count
    integer :: N_rejected, count
    integer :: i_obst
    integer :: rejected_place(1:in_obstacle_distribution%N-in_obstacle_distribution%N_placed)

    success=.false.
    N_rejected=0
    count=0

    !iterate while a placement is unsuccessful, the number rejected are less than the number unplaced,
    !and the iteration limit has not exceeded the maximum
    do while ((success .eqv. .false.) .and.&
         (N_rejected .lt. in_obstacle_distribution%N-in_obstacle_distribution%N_placed) .and.&
         (count .lt. in_deform_params%max_iter))

       !sample a random integer
       i_obst=int_rand(1,in_obstacle_distribution%N)

       !test if the proposed obstacle is an obstacle that has already been rejected, the cycle
       if (N_rejected .gt. 0) then
          call binary_search(lo,hi,rejected_place(1:N_rejected),N_rejected,i_obst)
          if ((lo .ne. -1) .and. (hi .ne. -1)) then
             count=count+1
             cycle
          end if
       end if

       !test if the proposed obstacle is an obstacle that has already been placed, then cycle
       if (in_obstacle_distribution%N_placed .gt. 0) then
          if (in_obstacle_distribution%placed(i_obst) .eq. 1) then
             count=count+1
             cycle
          end if
       end if

       !test if the proposed obstacle intersects the ring
       intersect_flag=test_obstacle_place(in_obstacle_distribution,in_ring,i_obst)

       !if the obstacle intersects the ring, then reject the move, increase the number rejected,
       !add it to the rejected list, sort the list, and increase the count
       if (intersect_flag .eqv. .true.) then
          N_rejected=N_rejected+1
          rejected_place(N_rejected)=i_obst
          call shell_sort(rejected_place(1:N_rejected),N_rejected)
          count=count+1
       else !if the obstacle does not intersect the ring, break the while loop, increase the number
          !placed, add the accepted move the placed list, and sort the placed list
          success=.true.
          in_obstacle_distribution%N_placed=in_obstacle_distribution%N_placed+1
          in_obstacle_distribution%placed(i_obst)=1
          call binary_search(lo,hi,&
               in_obstacle_distribution%mapping(1:in_obstacle_distribution%N),&
               in_obstacle_distribution%N,&
               in_obstacle_distribution%mapping(i_obst))
          
          if (in_obstacle_distribution%obst_type .eq. 1) then
             site_count=1
          elseif (in_obstacle_distribution%obst_type .eq. 2) then
             site_count=1+2*in_ring%dim
          elseif (in_obstacle_distribution%obst_type .eq. 3) then
             site_count=(2**in_ring%dim)*(1+2*in_ring%dim)
          else
             site_count=-1
          end if
          
          if (sum(in_obstacle_distribution%mapping(lo:hi)) .eq. site_count) then
             in_obstacle_distribution%obst_base%placed(in_obstacle_distribution%mapping(i_obst))=1
             in_obstacle_distribution%obst_base%N_placed=in_obstacle_distribution%obst_base%N_placed+1
          end if
       end if
       
    end do !finish iterating to sample the obstacle placements
    
  end subroutine place_obstacle

  !subroutine to attempt to place an obstacle
  subroutine place_obstacles(success,in_obstacle_distribution,in_ring,in_deform_params)

    implicit none

    logical, intent(out) :: success

    type(obstacle_distribution), intent(inout) :: in_obstacle_distribution

    type(ring), intent(in) :: in_ring
    type(deform_params), intent(in) :: in_deform_params

    logical :: intersect_flag

    integer :: lo, hi
    integer :: site_count
    integer :: N_rejected, N_accepted, N_obst_init, count
    integer :: i_obst
    integer :: rejected_place(1:in_obstacle_distribution%N-in_obstacle_distribution%N_placed)

    success=.false.
    N_rejected=0
    N_accepted=0
    count=0

    N_obst_init=in_obstacle_distribution%N_placed

    !iterate while a placement is unsuccessful, the number rejected are less than the number unplaced,
    !and the iteration limit has not exceeded the maximum
    do while ((count .lt. in_deform_params%max_iter) .and.&
         (N_rejected+N_accepted .lt. in_obstacle_distribution%N-N_obst_init) .and.&
         (N_accepted .lt. in_deform_params%obst_mult) .and.&
         (success .eqv. .false.))

       !sample a random integer
       i_obst=int_rand(1,in_obstacle_distribution%N)

       !test if the proposed obstacle is an obstacle that has already been rejected, then cycle
       if (N_rejected .gt. 0) then
          call binary_search(lo,hi,rejected_place(1:N_rejected),N_rejected,i_obst)
          if ((lo .ne. -1) .and. (hi .ne. -1)) then
             count=count+1
             cycle
          end if
       end if

       !test if the proposed obstacle is an obstacle that has already been placed, then cycle
       if (in_obstacle_distribution%N_placed .gt. 0) then
          if (in_obstacle_distribution%placed(i_obst) .eq. 1) then
             count=count+1
             cycle
          end if
       end if

       !test if the proposed obstacle intersects the ring
       intersect_flag=test_obstacle_place(in_obstacle_distribution,in_ring,i_obst)

       !if the obstacle intersects the ring, then reject the move, increase the number rejected,
       !add it to the rejected list, sort the list, and increase the count
       if (intersect_flag .eqv. .true.) then
          N_rejected=N_rejected+1
          rejected_place(N_rejected)=i_obst
          call shell_sort(rejected_place(1:N_rejected),N_rejected)
          count=count+1
       else !if the obstacle does not intersect the ring, break the while loop, increase the number
          !placed, add the accepted move the placed list, and sort the placed list
          in_obstacle_distribution%N_placed=in_obstacle_distribution%N_placed+1
          in_obstacle_distribution%placed(i_obst)=1
          call binary_search(lo,hi,&
               in_obstacle_distribution%mapping(1:in_obstacle_distribution%N),&
               in_obstacle_distribution%N,&
               in_obstacle_distribution%mapping(i_obst))
          
          if (in_obstacle_distribution%obst_type .eq. 1) then
             site_count=1
          elseif (in_obstacle_distribution%obst_type .eq. 2) then
             site_count=1+2*in_ring%dim
          elseif (in_obstacle_distribution%obst_type .eq. 3) then
             site_count=(2**in_ring%dim)*(1+2*in_ring%dim)
          else
             site_count=-1
          end if
          
          if (sum(in_obstacle_distribution%mapping(lo:hi)) .eq. site_count) then
             in_obstacle_distribution%obst_base%placed(in_obstacle_distribution%mapping(i_obst))=1
             in_obstacle_distribution%obst_base%N_placed=in_obstacle_distribution%obst_base%N_placed+1
          end if
       end if
       
    end do !finish iterating to sample the obstacle placements

    if ((N_accepted .eq. in_deform_params%obst_mult) .or.&
         (N_accepted .eq. in_obstacle_distribution%N-N_obst_init)) then
       success=.true.
    end if
    
  end subroutine place_obstacles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! GROW ROUTINES !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !subroutine to perform a grow based on the subset of a ring
  !returns if the grow is successful and only updates the ring if the grow is successful
  subroutine perform_grow(success,in_ring,in_grow)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring
    
    type(ring_subset), intent(in) :: in_grow
    
    logical :: intersect, bound


    intersect=.true.
    bound=.true.

    !test the boundaries for the grow
    call test_bounds(bound,in_ring,in_grow)

    !if the grow remains within the boundaries, then test the self-intersections
    if (bound .eqv. .false.) then
       call test_intersection_growth(intersect,in_ring,in_grow%N,in_grow%x(:,1:in_grow%N))
    end if

    !if the grow remains within the boundaries, and there are no self-intersetions,
    !then change the ring coordinates to the grow and indicate that the grow was
    !successful
    if ((intersect .eqv. .false.) .and. (bound .eqv. .false.)) then
       success=.true.
       call insert_node_pairs(in_ring,in_grow%N,in_grow%x,in_grow%node_indices)
    else !if the grow was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_grow

  !subroutine to perform a grow based on the subset of a ring within the obstacle distribution
  !returns if the grow is successful and only updates the ring if the grow is successful
  subroutine perform_grow_obstacles(success,in_ring,in_grow,in_obstacle_distribution)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring
    
    type(ring_subset), intent(in) :: in_grow
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    
    logical :: intersect, bound, obstacle
    

    intersect=.true.
    bound=.true.
    obstacle=.true.

    !test the boundaries for the grow
    call test_bounds(bound,in_ring,in_grow)

    !if the grow remains within the boundaries, then test the obstacle intersections
    if (bound .eqv. .false.) then
       call test_obstacle_move_fast(obstacle,in_obstacle_distribution,in_grow)
    end if

    !if the grow remains within the boundaries, and there are no obstacle intersections,
    !then test self-intersections
    if ((bound .eqv. .false.) .and.&
         (obstacle .eqv. .false.)) then
       !call test_intersection_mult(intersect,in_ring,in_grow)
       call test_intersection_growth(intersect,in_ring,in_grow%N,in_grow%x(:,1:in_grow%N))
    end if

    !if the grow remains within the boundaries, there are no self-intersections, and
    !there are no obstacle intersections, then change the ring coordinates to the grow
    !and indicate that the grow was successful
    if ((intersect .eqv. .false.) .and.&
         (obstacle .eqv. .false.).and.&
         (bound .eqv. .false.)) then
       success=.true.
       !call insert_node_pair(in_ring,in_grow%x,in_grow%node_indices)
       call insert_node_pairs(in_ring,in_grow%N,in_grow%x,in_grow%node_indices)
    else !if the grow was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_grow_obstacles

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! MOVE ROUTINES !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !subroutine to perform a move based on the subset of a ring
  !returns if the move is successful and only updates the ring if the move is successful
  pure subroutine perform_move(success,in_ring,in_move)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring
    
    type(ring_subset), intent(in) :: in_move
    
    logical :: intersect, bound
    integer :: i_node


    intersect=.true.
    bound=.true.

    !test the boundaries for the move
    call test_bounds(bound,in_ring,in_move)

    !if the move remains within the boundaries, then test the self-intersections
    if (bound .eqv. .false.) then
       call test_intersection_mult(intersect,in_ring,in_move)
    end if

    !if the move remains within the boundaries, and there are no self-intersetions,
    !then change the ring coordinates to the move and indicate that the move was
    !successful
    if ((intersect .eqv. .false.) .and. (bound .eqv. .false.)) then
       success=.true.
       do i_node=1,in_move%N
          in_ring%x(:,in_move%node_indices(i_node))=in_move%x(:,i_node)
       end do
    else !if the move was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_move

  !subroutine to perform a move based on the subset of a ring within the obstacle distribution
  !returns if the move is successful and only updates the ring if the move is successful
  !can be made pure if test_obstacles_move does not use octree
  subroutine perform_move_obstacles(success,in_ring,in_move,in_obstacle_distribution)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring
    
    type(ring_subset), intent(in) :: in_move
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    
    logical :: intersect, bound, obstacle
    integer :: i_node
    

    intersect=.true.
    bound=.true.
    obstacle=.true.

    !test the boundaries for the move
    call test_bounds(bound,in_ring,in_move)

    !if the move remains within the boundaries, then test the obstacle intersections
    if (bound .eqv. .false.) then
       call test_obstacle_move_fast(obstacle,in_obstacle_distribution,in_move)
    end if

    !if the move remains within the boundaries, and there are no obstacle intersections,
    !then test self-intersections
    if ((bound .eqv. .false.) .and.&
         (obstacle .eqv. .false.)) then
       call test_intersection_mult(intersect,in_ring,in_move)
    end if

    !if the move remains within the boundaries, there are no self-intersections, and
    !there are no obstacle intersections, then change the ring coordinates to the move
    !and indicate that the move was successful
    if ((intersect .eqv. .false.) .and.&
         (obstacle .eqv. .false.).and.&
         (bound .eqv. .false.)) then
       success=.true.
       do i_node=1,in_move%N
          in_ring%x(:,in_move%node_indices(i_node))=in_move%x(:,i_node)
       end do
    else !if the move was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_move_obstacles

  !subroutine to perform a move based on the subset of a ring
  !returns if the move is successful and only updates the ring if the move is successful
  subroutine perform_move_MH(success,in_ring,in_move,in_energy_params)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring
    
    type(ring_subset), intent(in) :: in_move
    type(energy_params), intent(in) :: in_energy_params
    
    logical :: intersect, bound
    
    integer :: i_node

    real(rp) :: temp_deltaE, total_deltaE, u

    total_deltaE=0.0d0

    intersect=.true.
    bound=.true.

    !test the boundaries for the move
    if (in_energy_params%V_coord_type .ne. 0) then
       call test_bounds_deltaE(bound,temp_deltaE,&
            in_ring,in_move,in_energy_params)
    else
       temp_deltaE=0.0d0
       call test_bounds(bound,in_ring,in_move)
    end if

    total_deltaE=total_deltaE+temp_deltaE
    
    !if the move remains within the boundaries, then test the self-intersections
    if (bound .eqv. .false.) then
       if ((in_energy_params%V_nn_type .ne. 0) .or.&
            (in_energy_params%V_pairwise_type .ne. 0)) then
          call test_intersection_mult_deltaE(intersect,temp_deltaE,&
               in_ring,in_move,in_energy_params)
       else
          temp_deltaE=0.0d0
          call test_intersection_mult(intersect,in_ring,in_move)
       end if
    end if

    total_deltaE=total_deltaE+temp_deltaE

    !if the move remains within the boundaries, and there are no self-intersetions,
    !then change the ring coordinates to the move and indicate that the move was
    !successful
    if ((intersect .eqv. .false.) .and. (bound .eqv. .false.)) then

       !calculate energy change due to bending
       if (in_energy_params%V_bend_type .ne. 0) then
          call calc_bend_deltaE(temp_deltaE,in_ring,&
               in_move,in_energy_params)
          total_deltaE=total_deltaE+temp_deltaE
       end if

       !calculate energy change due to twisting
       if (in_energy_params%V_twist_type .ne. 0) then
          call calc_twist_deltaE(temp_deltaE,in_ring,&
               in_move,in_energy_params)
          total_deltaE=total_deltaE+temp_deltaE
       end if
       
       call random_number(u)
       if (u .le. min(1.0d0,exp(-in_energy_params%beta*total_deltaE))) then
          success=.true.
          do i_node=1,in_move%N
             in_ring%x(:,in_move%node_indices(i_node))=in_move%x(:,i_node)
          end do
       else
          success=.false.
       end if
    else !if the move was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_move_MH

  !subroutine to perform a move based on the subset of a ring within the obstacle distribution
  !returns if the move is successful and only updates the ring if the move is successful
  subroutine perform_move_obstacles_MH(success,in_ring,&
       in_move,in_obstacle_distribution,in_energy_params)

    implicit none
    
    logical, intent(out) :: success
    
    type(ring), intent(inout) :: in_ring

    type(ring_subset), intent(in) :: in_move
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params
    
    logical :: intersect, bound, obstacle
    integer :: i_node

    real(rp) :: temp_deltaE, total_deltaE, u

    total_deltaE=0.0d0

    intersect=.true.
    bound=.true.
    obstacle=.true.

    !test the boundaries for the move
    if (in_energy_params%V_coord_type .ne. 0) then
       call test_bounds_deltaE(bound,temp_deltaE,&
            in_ring,in_move,in_energy_params)
    else
       temp_deltaE=0.0d0
       call test_bounds(bound,in_ring,in_move)
    end if

    total_deltaE=total_deltaE+temp_deltaE

    !if the move remains within the boundaries, then test the obstacle intersections
    if (bound .eqv. .false.) then
       if (in_energy_params%V_obst_type .ne. 0) then
          call test_obstacle_move_deltaE(obstacle,temp_deltaE,&
               in_obstacle_distribution,in_move,in_energy_params)
       else
          temp_deltaE=0.0d0
          call test_obstacle_move_fast(obstacle,in_obstacle_distribution,in_move)
       end if
    end if

    total_deltaE=total_deltaE+temp_deltaE
    
    !if the move remains within the boundaries, then test the self-intersections
    if ((bound .eqv. .false.) .and.&
         (obstacle .eqv. .false.)) then
       if ((in_energy_params%V_nn_type .ne. 0) .or.&
            (in_energy_params%V_pairwise_type .ne. 0)) then
          call test_intersection_mult_deltaE(intersect,temp_deltaE,&
               in_ring,in_move,in_energy_params)
       else
          temp_deltaE=0.0d0
          call test_intersection_mult(intersect,in_ring,in_move)
       end if
    end if

    total_deltaE=total_deltaE+temp_deltaE
    
    !if the move remains within the boundaries, there are no self-intersections, and
    !there are no obstacle intersections, then change the ring coordinates to the move
    !and indicate that the move was successful
    if ((intersect .eqv. .false.) .and.&
         (obstacle .eqv. .false.).and.&
         (bound .eqv. .false.)) then

       !calculate the energy change due to bending
       if (in_energy_params%V_bend_type .ne. 0) then
          call calc_bend_deltaE(temp_deltaE,in_ring,&
               in_move,in_energy_params)
          total_deltaE=total_deltaE+temp_deltaE
       end if

       !calculate energy change due to twisting
       if (in_energy_params%V_twist_type .ne. 0) then
          call calc_twist_deltaE(temp_deltaE,in_ring,&
               in_move,in_energy_params)
          total_deltaE=total_deltaE+temp_deltaE
       end if
       
       call random_number(u)
       if (u .le. min(1.0d0,exp(-in_energy_params%beta*total_deltaE))) then
          success=.true.
          do i_node=1,in_move%N
             in_ring%x(:,in_move%node_indices(i_node))=in_move%x(:,i_node)
          end do
       else
          success=.false.
       end if
    else !if the move was unsuccessful, indicate a failure
       success=.false.
    end if
    
  end subroutine perform_move_obstacles_MH
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MAIN SUBROUTINES !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !main routine used to generate ring configurations
  subroutine move_ring(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(grow_params), intent(in) :: in_grow_params
    type(deform_params), intent(in) :: in_deform_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(deform_params) :: restart_deform_params
    type(node_subsets) :: kink_moves
    type(direction_pairs) :: reflect_moves
    type(ring_subset) :: new_kink_move, new_reflect_move

    logical :: move_success
    logical :: length_error, intersect_error
    logical :: kink_cycle
    logical :: node_1_between, node_2_between, ref_cycle
    
    integer :: unit_rep
    integer :: write_update_freq
    integer :: i_iter, i_move, j_move, count, move_type
    integer :: N_moves, N_rejected, N_accepted, N_remaining, multiplicity
    integer :: rand_move_type
    integer :: initial_length, test_length
    integer :: fin_node, dir, i_node, j_node, mid_node
    integer :: config_test(1:2)
    integer :: ref_moves_accepted(1:2*in_deform_params%ref_mult)
    integer :: kink_moves_accepted(1:2*in_deform_params%kink_mult)
    integer :: incompatible_kink_nodes(1:2*in_deform_params%kink_mult)
    integer, allocatable :: moves_rejected(:), moves_remaining(:)

    character(len=120) :: format_string


    sim_completed=.false.
    
    write_update_freq=1

    move_type=0
    multiplicity=1
    
    if (in_deform_params%num_iter .le. 100) then
       write_update_freq=1
    elseif (in_deform_params%num_iter .le. 500) then
       write_update_freq=5
    elseif (in_deform_params%num_iter .le. 1000) then
       write_update_freq=10
    elseif (in_deform_params%num_iter .le. 5000) then
       write_update_freq=50
    elseif (in_deform_params%num_iter .le. 10000) then
       write_update_freq=100
    elseif (in_deform_params%num_iter .le. 50000) then
       write_update_freq=500
    elseif (in_deform_params%num_iter .le. 100000) then
       write_update_freq=1000
    elseif (in_deform_params%num_iter .le. 1000000) then
       write_update_freq=10000
    end if

    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    initial_length=ring_length(in_ring)

    write(unit_log,"(a,i5.5,a)")'[move_ring] --- i_rep=',i_rep,' --- BEGIN'
    
    !main loop
    !write(unit_log,*)' [move_ring] entering main loop'
    do i_iter=1,in_deform_params%num_iter

       !using the frequency information, select the type of move to be enumerated and
       !then sampled from
       if (in_deform_params%ergodic_flag .eq. 0) then
          if (i_iter .gt. in_deform_params%ref_wait) then
             !restrict move types based on frequency
             if (i_iter .le. in_deform_params%ref_init+in_deform_params%ref_wait) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             elseif (mod(i_iter,in_deform_params%ref_freq) .eq. 0) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             else
                !write(unit_log,*)' [move_ring_MH] only kink moves'
                kink_moves=count_kink_moves(in_ring)
                N_moves=kink_moves%N
                multiplicity=in_deform_params%kink_mult
                move_type=0
             end if
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             multiplicity=in_deform_params%kink_mult
             move_type=0
          end if
       elseif (in_deform_params%ergodic_flag .eq. 1) then
          rand_move_type=int_rand(1,in_deform_params%ref_freq)
          if (rand_move_type .eq. in_deform_params%ref_freq) then
             !write(unit_log,*)' [move_ring_MH] only reflect moves'
             multiplicity=int_rand(0,2*in_deform_params%ref_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             reflect_moves=count_reflect_moves(in_ring)
             N_moves=reflect_moves%N
             move_type=1
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             multiplicity=int_rand(0,2*in_deform_params%kink_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             move_type=0
          end if
       end if

       allocate(moves_rejected(1:N_moves),moves_remaining(1:N_moves))
       moves_rejected=0
       moves_remaining=(/(i_move,i_move=1,N_moves,1)/)

       !repeat multiple kink moves without renumerating the kinks
       if (move_type .eq. 0) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          kink_moves_accepted=0
          incompatible_kink_nodes=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)
             

             !write(unit_log,*)'     [move_ring] calculating kink move'
             new_kink_move=kink_move(in_ring,kink_moves%nodes(:,moves_remaining(i_move)))
             !write(unit_log,*)'     [move_ring] executing kink move'
             call perform_move(move_success,in_ring,new_kink_move)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=moves_remaining(i_move)
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                kink_moves_accepted(N_accepted)=moves_remaining(i_move)
                !incompatible_kink_nodes(N_accepted)=kink_moves%nodes(2,moves_remaining(i_move))

                fin_node=-1
                mid_node=kink_moves%nodes(2,moves_remaining(i_move))
                
                if (in_deform_params%kink_repeat .gt. 0) then
                   call repeat_kink_moves(fin_node,dir,&
                        in_ring,&
                        mid_node,&
                        in_deform_params%kink_repeat)
                end if

                
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)

                kink_cycle=.true.
                do while (kink_cycle .eqv. .true.)
                   kink_cycle=.false.
                   do j_move=1,N_remaining
                      if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                           mid_node))then
                         !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                         N_rejected=N_rejected+1
                         count=count+1
                         moves_rejected(N_rejected)=moves_remaining(j_move)
                         call reduce_remaining_moves(N_moves,j_move,&
                              moves_remaining,N_remaining)
                         kink_cycle=.true.
                         exit
                      end if
                   end do
                end do


                if (fin_node .ne. -1) then
                   if (dir .eq. -1) then
                      if (fin_node .lt. mid_node) then
                        ! write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node+in_ring%N_nodes-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   elseif (dir .eq. 1) then
                      if (fin_node .lt. mid_node) then
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-(fin_node+in_ring%N_nodes)-1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node+1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   end if
                end if

                move_success=.false.
             end if

          end do !end loop to find valid move

       !repeat multiple reflect moves without renumerating the reflects
       elseif (move_type .eq. 1) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          ref_moves_accepted=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)


             !test if reflect move is a repeat and if it is compatible with reflect moves
             if (N_accepted .gt. 0) then
                
                ref_cycle=.false.
                do j_move=1,N_accepted

                   node_1_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(1,i_move))

                   node_2_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(2,i_move))

                   !test for agreement with the betweenness
                   if (node_1_between .neqv. node_2_between) then
                      ref_cycle=.true.
                      exit
                   end if

                end do

                !cycle the while loop if incompatible
                if (ref_cycle .eqv. .true.) then

                   N_rejected=N_rejected+1
                   count=count+1
                   moves_rejected(N_rejected)=i_move
                   call reduce_remaining_moves(N_moves,i_move,&
                        moves_remaining,N_remaining)
                   cycle

                end if
                   
             end if

             
             !write(unit_log,*)'     [move_ring] calculating reflect move'
             new_reflect_move=reflect_move(in_ring,reflect_moves%nodes(:,i_move),&
                  reflect_moves%dirs(:,i_move))
             !write(unit_log,*)'     [move_ring] executing reflect move'
             call perform_move(move_success,in_ring,new_reflect_move)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                ref_moves_accepted(N_accepted)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
                move_success=.false.
             end if

          end do !end loop to find valid move

       end if

       if (allocated(kink_moves%nodes)) deallocate(kink_moves%nodes)
       if (allocated(new_kink_move%x)) deallocate(new_kink_move%x)
       if (allocated(new_kink_move%node_indices)) deallocate(new_kink_move%node_indices)
       if (allocated(reflect_moves%nodes)) deallocate(reflect_moves%nodes)
       if (allocated(reflect_moves%dirs)) deallocate(reflect_moves%dirs)
       if (allocated(new_reflect_move%x)) deallocate(new_reflect_move%x)
       if (allocated(new_reflect_move%node_indices)) deallocate(new_reflect_move%node_indices)
       if (allocated(moves_rejected)) deallocate(moves_rejected)
       if (allocated(moves_remaining)) deallocate(moves_remaining)


       !run additional error tests
       if ((mod(i_iter,in_deform_params%err_check_freq) .eq. 0) .or. &
            (i_iter .eq. in_deform_params%num_iter)) then

          length_error=.true.
          intersect_error=.true.
          
          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if (test_length .ne. initial_length) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring] --- i_rep=',i_rep,&
                  ' --- ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [move_ring] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring] --- i_rep=',i_rep,&
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

             !create a new set of deform params
             restart_deform_params=in_deform_params
             restart_deform_params%num_iter=restart_deform_params%num_iter-i_iter

             !adjust the deform params for the iterations completed so far
             if (i_iter .gt. restart_deform_params%ref_wait) then
                restart_deform_params%ref_wait=0
                if (i_iter .gt. restart_deform_params%ref_init+&
                     restart_deform_params%ref_wait) then
                   restart_deform_params%ref_init=0
                else
                   restart_deform_params%ref_init=restart_deform_params%ref_init+&
                        restart_deform_params%ref_wait-i_iter
                end if
             else
                restart_deform_params%ref_wait=restart_deform_params%ref_wait-i_iter
             end if

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  restart_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (i_iter .eq. in_deform_params%num_iter) then
                sim_completed=.true.
             end if

          end if
          
       end if !end conditional for error checking
       
       if (mod(i_iter,write_update_freq) .eq. 0) then
          format_string="(a,i5.5,a,i8,a,i8,a)"
          write(unit_log,format_string)' [move_ring] --- i_rep=',i_rep,&
               ' --- progress: ',&
               i_iter,' of ',in_deform_params%num_iter,' iterations completed'
       end if
       
    end do !end main loop

    write(unit_log,"(a,i5.5,a)")'[move_ring] --- i_rep=',i_rep,' --- END'
    
  end subroutine move_ring

  !main routine used to generate ring configurations, when there are obstacles present or
  !obstacles that need to placed during the motion of the ring
  subroutine move_ring_obstacles(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    type(obstacle_distribution), intent(inout) :: in_obstacle_distribution

    type(grow_params), intent(in) :: in_grow_params
    type(deform_params), intent(in) :: in_deform_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(deform_params) :: restart_deform_params
    type(node_subsets) :: kink_moves
    type(direction_pairs) :: reflect_moves
    type(ring_subset) :: new_kink_move, new_reflect_move

    logical :: move_success, obst_success
    logical :: length_error, intersect_error, obstacle_error
    logical :: kink_cycle
    logical :: node_1_between, node_2_between, ref_cycle
    
    integer :: unit_rep
    integer :: write_update_freq
    integer :: i_iter, i_move, j_move, count, move_type
    integer :: N_moves, N_rejected, N_accepted, N_remaining, multiplicity
    integer :: rand_move_type
    integer :: initial_length, test_length
    integer :: fin_node, dir, i_node, j_node, mid_node
    integer :: config_test(1:2), node_obst_test(1:2)
    integer :: ref_moves_accepted(1:2*in_deform_params%ref_mult)
    integer :: kink_moves_accepted(1:2*in_deform_params%kink_mult)
    integer :: incompatible_kink_nodes(1:2*in_deform_params%kink_mult)
    integer, allocatable :: moves_rejected(:), moves_remaining(:)

    character(len=120) :: format_string


    sim_completed=.false.

    write_update_freq=1

    move_type=0
    multiplicity=1
    
    if (in_deform_params%num_iter .le. 100) then
       write_update_freq=1
    elseif (in_deform_params%num_iter .le. 500) then
       write_update_freq=5
    elseif (in_deform_params%num_iter .le. 1000) then
       write_update_freq=10
    elseif (in_deform_params%num_iter .le. 5000) then
       write_update_freq=50
    elseif (in_deform_params%num_iter .le. 10000) then
       write_update_freq=100
    elseif (in_deform_params%num_iter .le. 50000) then
       write_update_freq=500
    elseif (in_deform_params%num_iter .le. 100000) then
       write_update_freq=1000
    elseif (in_deform_params%num_iter .le. 1000000) then
       write_update_freq=10000
    end if

    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    initial_length=ring_length(in_ring)

    write(unit_log,"(a,i5.5,a)")'[move_ring_obstacles] --- i_rep=',i_rep,' --- BEGIN'
    
    !write(unit_log,*)' [move_ring_obstacles] entering main loop'
    !main loop
    do i_iter=1,in_deform_params%num_iter

       !if the obstacle frequency is matched, then attempt to place an obstacle from
       !the list of unplaced obstacles
       if (in_deform_params%obst_freq .ne. 0) then
          
          if ((mod(i_iter,in_deform_params%obst_freq) .eq. 0) .and.&
               (in_obstacle_distribution%N_placed .lt. in_obstacle_distribution%N)) then
             obst_success=.false.
             call place_obstacles(obst_success,in_obstacle_distribution,&
                  in_ring,in_deform_params)
             if (obst_success .eqv. .true.) then
                write(unit_log,"(a,i5.5,a,i6.6,a)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                     ' --- obstacle success, ',in_obstacle_distribution%N_placed,&
                     ' obstacles now placed'
             else
                write(unit_log,"(a,i5.5,a,i6.6,a)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                     ' --- obstacle failure, ',in_obstacle_distribution%N_placed,&
                     ' obstacles already placed'
             end if
          end if

          if ((in_obstacle_distribution%N_placed .eq. in_obstacle_distribution%N) .and.&
               (i_iter .lt. in_deform_params%num_iter)) then
             cycle
          end if
          
       end if

       !using the frequency information, select the type of move to be enumerated and
       !then sampled from
       if (in_deform_params%ergodic_flag .eq. 0) then
          if (i_iter .gt. in_deform_params%ref_wait) then
             !restrict move types based on frequency
             if (i_iter .le. in_deform_params%ref_init+in_deform_params%ref_wait) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             elseif (mod(i_iter,in_deform_params%ref_freq) .eq. 0) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             else
                !write(unit_log,*)' [move_ring_MH] only kink moves'
                kink_moves=count_kink_moves(in_ring)
                N_moves=kink_moves%N
                multiplicity=in_deform_params%kink_mult
                move_type=0
             end if
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             multiplicity=in_deform_params%kink_mult
             move_type=0
          end if
       elseif (in_deform_params%ergodic_flag .eq. 1) then
          rand_move_type=int_rand(1,in_deform_params%ref_freq)
          if (rand_move_type .eq. in_deform_params%ref_freq) then
             !write(unit_log,*)' [move_ring_MH] only reflect moves'
             multiplicity=int_rand(0,2*in_deform_params%ref_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             reflect_moves=count_reflect_moves(in_ring)
             N_moves=reflect_moves%N
             move_type=1
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             multiplicity=int_rand(0,2*in_deform_params%kink_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             move_type=0
          end if
       end if


       allocate(moves_rejected(1:N_moves),moves_remaining(1:N_moves))
       moves_rejected=0
       moves_remaining=(/(i_move,i_move=1,N_moves,1)/)

       !repeat multiple kink moves without renumerating the kinks
       if (move_type .eq. 0) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          kink_moves_accepted=0
          incompatible_kink_nodes=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)
             

             !write(unit_log,*)'     [move_ring] calculating kink move'
             new_kink_move=kink_move(in_ring,kink_moves%nodes(:,moves_remaining(i_move)))
             !write(unit_log,*)'     [move_ring] executing kink move'
             call perform_move_obstacles(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=moves_remaining(i_move)
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                kink_moves_accepted(N_accepted)=moves_remaining(i_move)
                !incompatible_kink_nodes(N_accepted)=kink_moves%nodes(2,moves_remaining(i_move))

                fin_node=-1
                mid_node=kink_moves%nodes(2,moves_remaining(i_move))
                
                if (in_deform_params%kink_repeat .gt. 0) then
                   call repeat_kink_moves_obstacles(fin_node,dir,&
                        in_ring,&
                        mid_node,&
                        in_deform_params%kink_repeat,&
                        in_obstacle_distribution)
                end if

                
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)

                kink_cycle=.true.
                do while (kink_cycle .eqv. .true.)
                   kink_cycle=.false.
                   do j_move=1,N_remaining
                      if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                           mid_node))then
                         !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                         N_rejected=N_rejected+1
                         count=count+1
                         moves_rejected(N_rejected)=moves_remaining(j_move)
                         call reduce_remaining_moves(N_moves,j_move,&
                              moves_remaining,N_remaining)
                         kink_cycle=.true.
                         exit
                      end if
                   end do
                end do


                if (fin_node .ne. -1) then
                   if (dir .eq. -1) then
                      if (fin_node .lt. mid_node) then
                        ! write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node+in_ring%N_nodes-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   elseif (dir .eq. 1) then
                      if (fin_node .lt. mid_node) then
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-(fin_node+in_ring%N_nodes)-1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node+1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   end if
                end if

                move_success=.false.
             end if

          end do !end loop to find valid move

       !repeat multiple reflect moves without renumerating the reflects
       elseif (move_type .eq. 1) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          ref_moves_accepted=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)


             !test if reflect move is a repeat and if it is compatible with reflect moves
             if (N_accepted .gt. 0) then
                
                ref_cycle=.false.
                do j_move=1,N_accepted

                   node_1_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(1,i_move))

                   node_2_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(2,i_move))

                   !test for agreement with the betweenness
                   if (node_1_between .neqv. node_2_between) then
                      ref_cycle=.true.
                      exit
                   end if

                end do

                !cycle the while loop if incompatible
                if (ref_cycle .eqv. .true.) then

                   N_rejected=N_rejected+1
                   count=count+1
                   moves_rejected(N_rejected)=i_move
                   call reduce_remaining_moves(N_moves,i_move,&
                        moves_remaining,N_remaining)
                   cycle

                end if
                   
             end if

             
             !write(unit_log,*)'     [move_ring] calculating reflect move'
             new_reflect_move=reflect_move(in_ring,reflect_moves%nodes(:,i_move),&
                  reflect_moves%dirs(:,i_move))
             !write(unit_log,*)'     [move_ring] executing reflect move'
             call perform_move_obstacles(move_success,in_ring,new_reflect_move,&
                  in_obstacle_distribution)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                ref_moves_accepted(N_accepted)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
                move_success=.false.
             end if

          end do !end loop to find valid move

       end if

       if (allocated(kink_moves%nodes)) deallocate(kink_moves%nodes)
       if (allocated(new_kink_move%x)) deallocate(new_kink_move%x)
       if (allocated(new_kink_move%node_indices)) deallocate(new_kink_move%node_indices)
       if (allocated(reflect_moves%nodes)) deallocate(reflect_moves%nodes)
       if (allocated(reflect_moves%dirs)) deallocate(reflect_moves%dirs)
       if (allocated(new_reflect_move%x)) deallocate(new_reflect_move%x)
       if (allocated(new_reflect_move%node_indices)) deallocate(new_reflect_move%node_indices)
       if (allocated(moves_rejected)) deallocate(moves_rejected)
       if (allocated(moves_remaining)) deallocate(moves_remaining)

       !run additional error tests
       if ((mod(i_iter,in_deform_params%err_check_freq) .eq. 0) .or. &
            (i_iter .eq. in_deform_params%num_iter)) then

          length_error=.true.
          intersect_error=.true.
          obstacle_error=.true.
          
          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if (test_length .ne. initial_length) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- interesection at nodes',&
                  config_test(1),' and ',config_test(2)
             exit
          else
             intersect_error=.false.
          end if

          !test for obstacle intersections
          node_obst_test=obstacle_test_fast(in_ring,in_obstacle_distribution)
          !node_obst_test=obstacle_test(in_ring,in_obstacle_distribution)
          
          if (sum(node_obst_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring-obstacle intersection'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring_obstacles] --- i_rep=',i_rep,&
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

             !create a new set of deform params
             restart_deform_params=in_deform_params
             restart_deform_params%num_iter=restart_deform_params%num_iter-i_iter

             !adjust the deform params for the iterations completed so far
             if (i_iter .gt. restart_deform_params%ref_wait) then
                restart_deform_params%ref_wait=0
                if (i_iter .gt. restart_deform_params%ref_init+&
                     restart_deform_params%ref_wait) then
                   restart_deform_params%ref_init=0
                else
                   restart_deform_params%ref_init=restart_deform_params%ref_init+&
                        restart_deform_params%ref_wait-i_iter
                end if
             else
                restart_deform_params%ref_wait=restart_deform_params%ref_wait-i_iter
             end if

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  restart_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (i_iter .eq. in_deform_params%num_iter) then
                sim_completed=.true.
             end if
             
          end if
          
       end if !end conditional for error checking
       
       if (mod(i_iter,write_update_freq) .eq. 0) then
          format_string="(a,i5.5,a,i8,a,i8,a)"
          write(unit_log,format_string)' [move_ring_obstacles] --- i_rep=',i_rep,&
               ' --- progress: ',&
               i_iter,' of ',in_deform_params%num_iter,' iterations completed'
       end if
       
    end do !end main loop

    write(unit_log,"(a,i5.5,a)")'[move_ring_obstacles] --- i_rep=',i_rep,' --- END'
    
  end subroutine move_ring_obstacles

  !main routine used to generate ring configurations
  subroutine move_ring_MH(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_energy_params,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    
    type(obstacle_distribution), intent(in) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params

    type(grow_params), intent(in) :: in_grow_params
    type(deform_params), intent(in) :: in_deform_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(deform_params) :: restart_deform_params
    type(node_subsets) :: kink_moves
    type(direction_pairs) :: reflect_moves
    type(ring_subset) :: new_kink_move, new_reflect_move

    logical :: move_success
    logical :: length_error, intersect_error
    logical :: kink_cycle
    logical :: node_1_between, node_2_between, ref_cycle
    
    integer :: unit_rep
    integer :: write_update_freq
    integer :: i_iter, i_move, j_move, count, move_type
    integer :: N_moves, N_rejected, N_accepted, N_remaining, multiplicity
    integer :: rand_move_type
    integer :: initial_length, test_length
    integer :: fin_node, dir, i_node, j_node, mid_node
    integer :: nodes_moved
    integer :: config_test(1:2)
    integer :: ref_moves_accepted(1:2*in_deform_params%ref_mult)
    integer :: kink_moves_accepted(1:2*in_deform_params%kink_mult)
    integer :: incompatible_kink_nodes(1:2*in_deform_params%kink_mult)
    integer, allocatable :: moves_rejected(:), moves_remaining(:)

    real(rp) :: temp_energy

    character(len=120) :: format_string
    character(len=120) :: energy_file


    sim_completed=.false.
    
    write_update_freq=1

    move_type=0
    multiplicity=1
    
    if (in_deform_params%num_iter .le. 100) then
       write_update_freq=1
    elseif (in_deform_params%num_iter .le. 500) then
       write_update_freq=5
    elseif (in_deform_params%num_iter .le. 1000) then
       write_update_freq=10
    elseif (in_deform_params%num_iter .le. 5000) then
       write_update_freq=50
    elseif (in_deform_params%num_iter .le. 10000) then
       write_update_freq=100
    elseif (in_deform_params%num_iter .le. 50000) then
       write_update_freq=500
    elseif (in_deform_params%num_iter .le. 100000) then
       write_update_freq=1000
    elseif (in_deform_params%num_iter .le. 1000000) then
       write_update_freq=10000
    end if

    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    initial_length=ring_length(in_ring)

    nodes_moved=0

    if (track_energy .eqv. .true.) then
       call total_energy(temp_energy,in_ring,in_energy_params)
       write(unit_log,"(a,i5.5,a,es20.10)")'[move_ring_MH] --- i_rep=',i_rep,' --- BEGIN, energy = ',temp_energy
       energy_file=adjustl(restart_file(:index(restart_file,".")-9))
       energy_file=trim(energy_file)//"_energy.dat"
       open(UNIT=unit_rep,FILE=trim(energy_file),ACTION='write',STATUS='unknown',POSITION='append')
       write(unit_rep,"(i9,a,es20.10)")nodes_moved,',',temp_energy
       close(unit_rep)
    else
       write(unit_log,"(a,i5.5,a)")'[move_ring_MH] --- i_rep=',i_rep,' --- BEGIN'
    end if
    
    !main loop
    !write(unit_log,*)' [move_ring_MH] entering main loop'
    do i_iter=1,in_deform_params%num_iter

       !using the frequency information, select the type of move to be enumerated and
       !then sampled from
       if (in_deform_params%ergodic_flag .eq. 0) then
          if (i_iter .gt. in_deform_params%ref_wait) then
             !restrict move types based on frequency
             if (i_iter .le. in_deform_params%ref_init+in_deform_params%ref_wait) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             elseif (mod(i_iter,in_deform_params%ref_freq) .eq. 0) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             else
                !write(unit_log,*)' [move_ring_MH] only kink moves'
                kink_moves=count_kink_moves(in_ring)
                N_moves=kink_moves%N
                multiplicity=in_deform_params%kink_mult
                move_type=0
             end if
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             multiplicity=in_deform_params%kink_mult
             move_type=0
          end if
       elseif (in_deform_params%ergodic_flag .eq. 1) then
          rand_move_type=int_rand(1,in_deform_params%ref_freq)
          if (rand_move_type .eq. in_deform_params%ref_freq) then
             !write(unit_log,*)' [move_ring_MH] only reflect moves'
             multiplicity=int_rand(0,2*in_deform_params%ref_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             reflect_moves=count_reflect_moves(in_ring)
             N_moves=reflect_moves%N
             move_type=1
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             multiplicity=int_rand(0,2*in_deform_params%kink_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             move_type=0
          end if
       end if

       allocate(moves_rejected(1:N_moves),moves_remaining(1:N_moves))
       moves_rejected=0
       moves_remaining=(/(i_move,i_move=1,N_moves,1)/)

       !repeat multiple kink moves without renumerating the kinks
       if (move_type .eq. 0) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          kink_moves_accepted=0
          incompatible_kink_nodes=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)
             

             !write(unit_log,*)'     [move_ring] calculating kink move'
             new_kink_move=kink_move(in_ring,kink_moves%nodes(:,moves_remaining(i_move)))
             !write(unit_log,*)'     [move_ring] executing kink move'
             call perform_move_MH(move_success,in_ring,new_kink_move,&
                  in_energy_params)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=moves_remaining(i_move)
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                kink_moves_accepted(N_accepted)=moves_remaining(i_move)
                !incompatible_kink_nodes(N_accepted)=kink_moves%nodes(2,moves_remaining(i_move))

                nodes_moved=nodes_moved+new_kink_move%N

                fin_node=-1
                mid_node=kink_moves%nodes(2,moves_remaining(i_move))
                
                if (in_deform_params%kink_repeat .gt. 0) then
                   call repeat_kink_moves_MH(fin_node,dir,&
                        in_ring,&
                        mid_node,&
                        in_deform_params%kink_repeat,&
                        in_energy_params)
                end if

                
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)

                kink_cycle=.true.
                do while (kink_cycle .eqv. .true.)
                   kink_cycle=.false.
                   do j_move=1,N_remaining
                      if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                           mid_node))then
                         !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                         N_rejected=N_rejected+1
                         count=count+1
                         moves_rejected(N_rejected)=moves_remaining(j_move)
                         call reduce_remaining_moves(N_moves,j_move,&
                              moves_remaining,N_remaining)
                         kink_cycle=.true.
                         exit
                      end if
                   end do
                end do


                if (fin_node .ne. -1) then
                   if (dir .eq. -1) then
                      if (fin_node .lt. mid_node) then
                        ! write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node+in_ring%N_nodes-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   elseif (dir .eq. 1) then
                      if (fin_node .lt. mid_node) then
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-(fin_node+in_ring%N_nodes)-1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node+1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   end if
                end if

                move_success=.false.
             end if

          end do !end loop to find valid move

       !repeat multiple reflect moves without renumerating the reflects
       elseif (move_type .eq. 1) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          ref_moves_accepted=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)


             !test if reflect move is a repeat and if it is compatible with reflect moves
             if (N_accepted .gt. 0) then
                
                ref_cycle=.false.
                do j_move=1,N_accepted

                   node_1_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(1,i_move))

                   node_2_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(2,i_move))

                   !test for agreement with the betweenness
                   if (node_1_between .neqv. node_2_between) then
                      ref_cycle=.true.
                      exit
                   end if

                end do

                !cycle the while loop if incompatible
                if (ref_cycle .eqv. .true.) then

                   N_rejected=N_rejected+1
                   count=count+1
                   moves_rejected(N_rejected)=i_move
                   call reduce_remaining_moves(N_moves,i_move,&
                        moves_remaining,N_remaining)
                   cycle

                end if
                   
             end if

             
             !write(unit_log,*)'     [move_ring] calculating reflect move'
             new_reflect_move=reflect_move(in_ring,reflect_moves%nodes(:,i_move),&
                  reflect_moves%dirs(:,i_move))
             !write(unit_log,*)'     [move_ring] executing reflect move'
             call perform_move_MH(move_success,in_ring,new_reflect_move,&
                  in_energy_params)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                ref_moves_accepted(N_accepted)=i_move
                nodes_moved=nodes_moved+new_reflect_move%N
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
                move_success=.false.
             end if

          end do !end loop to find valid move

       end if

       if (allocated(kink_moves%nodes)) deallocate(kink_moves%nodes)
       if (allocated(new_kink_move%x)) deallocate(new_kink_move%x)
       if (allocated(new_kink_move%node_indices)) deallocate(new_kink_move%node_indices)
       if (allocated(reflect_moves%nodes)) deallocate(reflect_moves%nodes)
       if (allocated(reflect_moves%dirs)) deallocate(reflect_moves%dirs)
       if (allocated(new_reflect_move%x)) deallocate(new_reflect_move%x)
       if (allocated(new_reflect_move%node_indices)) deallocate(new_reflect_move%node_indices)
       if (allocated(moves_rejected)) deallocate(moves_rejected)
       if (allocated(moves_remaining)) deallocate(moves_remaining)


       !run additional error tests
       if ((mod(i_iter,in_deform_params%err_check_freq) .eq. 0) .or. &
            (i_iter .eq. in_deform_params%num_iter)) then

          length_error=.true.
          intersect_error=.true.
          
          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if (test_length .ne. initial_length) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_MH] --- i_rep=',i_rep,&
                  ' --- ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_MH] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_MH] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [move_ring_MH] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_MH] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring_MH] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring_MH] --- i_rep=',i_rep,&
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

             !create a new set of deform params
             restart_deform_params=in_deform_params
             restart_deform_params%num_iter=restart_deform_params%num_iter-i_iter

             !adjust the deform params for the iterations completed so far
             if (i_iter .gt. restart_deform_params%ref_wait) then
                restart_deform_params%ref_wait=0
                if (i_iter .gt. restart_deform_params%ref_init+&
                     restart_deform_params%ref_wait) then
                   restart_deform_params%ref_init=0
                else
                   restart_deform_params%ref_init=restart_deform_params%ref_init+&
                        restart_deform_params%ref_wait-i_iter
                end if
             else
                restart_deform_params%ref_wait=restart_deform_params%ref_wait-i_iter
             end if

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  restart_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (track_energy .eqv. .true.) then
                format_string="(a,i5.5,a)"
                !write(unit_log,format_string)'  [move_ring_MH] --- i_rep=',i_rep,' --- calculate energy'
                call total_energy(temp_energy,in_ring,in_energy_params)
                open(UNIT=unit_rep,FILE=trim(energy_file),ACTION='write',STATUS='unknown',POSITION='append')
                write(unit_rep,"(i9,a,es20.10)")nodes_moved,',',temp_energy
                close(unit_rep)
             end if

             if (i_iter .eq. in_deform_params%num_iter) then
                sim_completed=.true.
             end if

          end if
          
       end if !end conditional for error checking
       
       if (mod(i_iter,write_update_freq) .eq. 0) then
          format_string="(a,i5.5,a,i8,a,i8,a)"
          write(unit_log,format_string)' [move_ring_MH] --- i_rep=',i_rep,&
               ' --- progress: ',&
               i_iter,' of ',in_deform_params%num_iter,' iterations completed'
       end if
       
    end do !end main loop

    if (track_energy .eqv. .true.) then
       ! call total_energy(temp_energy,in_ring,in_energy_params)
       write(unit_log,"(a,i5.5,a,es20.10,a,i12)")'[move_ring_MH] --- i_rep=',i_rep,&
            ' --- END, energy = ',temp_energy,', nodes_moved = ',nodes_moved
       ! open(UNIT=unit_rep,FILE=trim(energy_file),ACTION='write',STATUS='unknown',POSITION='append')
       ! write(unit_rep,"(i9,a,es20.10)")nodes_moved,',',temp_energy
       ! close(unit_rep)
    else
       write(unit_log,"(a,i5.5,a)")'[move_ring_MH] --- i_rep=',i_rep,' --- END'
    end if
    
  end subroutine move_ring_MH

  !main routine used to generate ring configurations, when there are obstacles present or
  !obstacles that need to placed during the motion of the ring
  subroutine move_ring_obstacles_MH(sim_completed,&
       in_ring,in_obstacle_distribution,&
       in_energy_params,&
       in_deform_params,in_grow_params,&
       restart_file,i_rep)

    implicit none

    logical, intent(out) :: sim_completed
    
    type(ring), intent(inout) :: in_ring
    type(obstacle_distribution), intent(inout) :: in_obstacle_distribution
    type(energy_params), intent(in) :: in_energy_params

    type(grow_params), intent(in) :: in_grow_params
    type(deform_params), intent(in) :: in_deform_params
    integer, intent(in) :: i_rep
    character(len=120), intent(in) :: restart_file

    type(deform_params) :: restart_deform_params
    type(node_subsets) :: kink_moves
    type(direction_pairs) :: reflect_moves
    type(ring_subset) :: new_kink_move, new_reflect_move

    logical :: move_success, obst_success
    logical :: length_error, intersect_error, obstacle_error
    logical :: kink_cycle
    logical :: node_1_between, node_2_between, ref_cycle
    
    integer :: unit_rep
    integer :: write_update_freq
    integer :: i_iter, i_move, j_move, count, move_type
    integer :: N_moves, N_rejected, N_accepted, N_remaining, multiplicity
    integer :: rand_move_type
    integer :: initial_length, test_length
    integer :: fin_node, dir, i_node, j_node, mid_node
    integer :: nodes_moved
    integer :: config_test(1:2), node_obst_test(1:2)
    integer :: ref_moves_accepted(1:2*in_deform_params%ref_mult)
    integer :: kink_moves_accepted(1:2*in_deform_params%kink_mult)
    integer :: incompatible_kink_nodes(1:2*in_deform_params%kink_mult)
    integer, allocatable :: moves_rejected(:), moves_remaining(:)

    real(rp) :: temp_energy

    character(len=120) :: format_string
    character(len=120) :: energy_file


    sim_completed=.false.

    write_update_freq=1

    move_type=0
    multiplicity=1
    
    if (in_deform_params%num_iter .le. 100) then
       write_update_freq=1
    elseif (in_deform_params%num_iter .le. 500) then
       write_update_freq=5
    elseif (in_deform_params%num_iter .le. 1000) then
       write_update_freq=10
    elseif (in_deform_params%num_iter .le. 5000) then
       write_update_freq=50
    elseif (in_deform_params%num_iter .le. 10000) then
       write_update_freq=100
    elseif (in_deform_params%num_iter .le. 50000) then
       write_update_freq=500
    elseif (in_deform_params%num_iter .le. 100000) then
       write_update_freq=1000
    elseif (in_deform_params%num_iter .le. 1000000) then
       write_update_freq=10000
    end if

    unit_rep=offset_units+(i_rep-in_deform_params%min_rep)

    initial_length=ring_length(in_ring)

    nodes_moved=0

    if (track_energy .eqv. .true.) then
       call total_energy(temp_energy,in_ring,in_energy_params)
       write(unit_log,"(a,i5.5,a,es20.10)")'[move_ring_obstacles_MH] --- i_rep=',i_rep,' --- BEGIN, energy = ',temp_energy
       energy_file=adjustl(restart_file(:index(restart_file,".")-9))
       energy_file=trim(energy_file)//"_energy.dat"
       open(UNIT=unit_rep,FILE=trim(energy_file),ACTION='write',STATUS='unknown',POSITION='append')
       write(unit_rep,"(i9,a,es20.10)")nodes_moved,',',temp_energy
       close(unit_rep)
    else
       write(unit_log,"(a,i5.5,a)")'[move_ring_obstacles_MH] --- i_rep=',i_rep,' --- BEGIN'
    end if
    
    !write(unit_log,*)' [move_ring_obstacles_MH] entering main loop'
    !main loop
    do i_iter=1,in_deform_params%num_iter

       !if the obstacle frequency is matched, then attempt to place an obstacle from
       !the list of unplaced obstacles
       if (in_deform_params%obst_freq .ne. 0) then
          
          if ((mod(i_iter,in_deform_params%obst_freq) .eq. 0) .and.&
               (in_obstacle_distribution%N_placed .lt. in_obstacle_distribution%N)) then
             obst_success=.false.
             call place_obstacles(obst_success,in_obstacle_distribution,&
                  in_ring,in_deform_params)
             if (obst_success .eqv. .true.) then
                write(unit_log,"(a,i5.5,a,i6.6,a)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                     ' --- obstacle success, ',in_obstacle_distribution%N_placed,&
                     ' obstacles now placed'
             else
                write(unit_log,"(a,i5.5,a,i6.6,a)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                     ' --- obstacle failure, ',in_obstacle_distribution%N_placed,&
                     ' obstacles already placed'
             end if
          end if

          if ((in_obstacle_distribution%N_placed .eq. in_obstacle_distribution%N) .and.&
               (i_iter .lt. in_deform_params%num_iter)) then
             cycle
          end if
          
       end if

       !using the frequency information, select the type of move to be enumerated and
       !then sampled from
       if (in_deform_params%ergodic_flag .eq. 0) then
          if (i_iter .gt. in_deform_params%ref_wait) then
             !restrict move types based on frequency
             if (i_iter .le. in_deform_params%ref_init+in_deform_params%ref_wait) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             elseif (mod(i_iter,in_deform_params%ref_freq) .eq. 0) then
                !write(unit_log,*)' [move_ring_MH] only reflect moves'
                reflect_moves=count_reflect_moves(in_ring)
                N_moves=reflect_moves%N
                multiplicity=in_deform_params%ref_mult
                move_type=1
             else
                !write(unit_log,*)' [move_ring_MH] only kink moves'
                kink_moves=count_kink_moves(in_ring)
                N_moves=kink_moves%N
                multiplicity=in_deform_params%kink_mult
                move_type=0
             end if
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             multiplicity=in_deform_params%kink_mult
             move_type=0
          end if
       elseif (in_deform_params%ergodic_flag .eq. 1) then
          rand_move_type=int_rand(1,in_deform_params%ref_freq)
          if (rand_move_type .eq. in_deform_params%ref_freq) then
             !write(unit_log,*)' [move_ring_MH] only reflect moves'
             multiplicity=int_rand(0,2*in_deform_params%ref_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             reflect_moves=count_reflect_moves(in_ring)
             N_moves=reflect_moves%N
             move_type=1
          else
             !write(unit_log,*)' [move_ring_MH] only kink moves'
             multiplicity=int_rand(0,2*in_deform_params%kink_mult)
             if (multiplicity .eq. 0) then
                cycle
             end if
             kink_moves=count_kink_moves(in_ring)
             N_moves=kink_moves%N
             move_type=0
          end if
       end if

       
       allocate(moves_rejected(1:N_moves),moves_remaining(1:N_moves))
       moves_rejected=0
       moves_remaining=(/(i_move,i_move=1,N_moves,1)/)

       !repeat multiple kink moves without renumerating the kinks
       if (move_type .eq. 0) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          kink_moves_accepted=0
          incompatible_kink_nodes=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)
             

             !write(unit_log,*)'     [move_ring] calculating kink move'
             new_kink_move=kink_move(in_ring,kink_moves%nodes(:,moves_remaining(i_move)))
             !write(unit_log,*)'     [move_ring] executing kink move'
             call perform_move_obstacles_MH(move_success,in_ring,new_kink_move,&
                  in_obstacle_distribution,in_energy_params)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=moves_remaining(i_move)
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                kink_moves_accepted(N_accepted)=moves_remaining(i_move)
                !incompatible_kink_nodes(N_accepted)=kink_moves%nodes(2,moves_remaining(i_move))

                nodes_moved=nodes_moved+new_kink_move%N
                
                fin_node=-1
                mid_node=kink_moves%nodes(2,moves_remaining(i_move))
                
                if (in_deform_params%kink_repeat .gt. 0) then
                   call repeat_kink_moves_obstacles_MH(fin_node,dir,&
                        in_ring,&
                        mid_node,&
                        in_deform_params%kink_repeat,&
                        in_obstacle_distribution,in_energy_params)
                end if

                
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)

                kink_cycle=.true.
                do while (kink_cycle .eqv. .true.)
                   kink_cycle=.false.
                   do j_move=1,N_remaining
                      if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                           mid_node) .or.&
                           (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                           mid_node))then
                         !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                         N_rejected=N_rejected+1
                         count=count+1
                         moves_rejected(N_rejected)=moves_remaining(j_move)
                         call reduce_remaining_moves(N_moves,j_move,&
                              moves_remaining,N_remaining)
                         kink_cycle=.true.
                         exit
                      end if
                   end do
                end do


                if (fin_node .ne. -1) then
                   if (dir .eq. -1) then
                      if (fin_node .lt. mid_node) then
                        ! write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node+in_ring%N_nodes-fin_node-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   elseif (dir .eq. 1) then
                      if (fin_node .lt. mid_node) then
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-(fin_node+in_ring%N_nodes)-1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      else
                         !write(unit_log,*)"mid_node=",mid_node,"fin_node=",fin_node,"dir=",dir
                         do i_node=0,mid_node-fin_node+1,-1
                            j_node=count_transform_alt(i_node,fin_node,in_ring%N_nodes)
                            !write(unit_log,*)"j_node=",j_node
                            kink_cycle=.true.
                            do while (kink_cycle .eqv. .true.)
                               kink_cycle=.false.
                               do j_move=1,N_remaining
                                  if ((kink_moves%nodes(1,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(2,moves_remaining(j_move)) .eq.&
                                       j_node) .or.&
                                       (kink_moves%nodes(3,moves_remaining(j_move)) .eq.&
                                       j_node))then
                                     !write(unit_log,*)kink_moves%nodes(:,moves_remaining(j_move))
                                     N_rejected=N_rejected+1
                                     count=count+1
                                     moves_rejected(N_rejected)=moves_remaining(j_move)
                                     call reduce_remaining_moves(N_moves,j_move,&
                                          moves_remaining,N_remaining)
                                     kink_cycle=.true.
                                     exit
                                  end if
                               end do
                            end do
                         end do
                      end if
                   end if
                end if

                move_success=.false.
             end if

          end do !end loop to find valid move

       !repeat multiple reflect moves without renumerating the reflects
       elseif (move_type .eq. 1) then

          move_success=.false.
          count=0
          N_rejected=0
          N_accepted=0
          N_remaining=N_moves

          ref_moves_accepted=0
          
          
          !loop to find a valid move
          !write(unit_log,*)'  [move_ring] entering sampling loop'
          do while ((count .lt. in_deform_params%max_iter) .and.&
               (N_rejected+N_accepted .lt. N_moves) .and.&
               (N_accepted .lt. multiplicity) .and.&
               (move_success .eqv. .false.))

             !write(unit_log,*)'   [move_ring] generating new sample move'
             !generate new move
             i_move=int_rand(1,N_remaining)


             !test if reflect move is a repeat and if it is compatible with reflect moves
             if (N_accepted .gt. 0) then
                
                ref_cycle=.false.
                do j_move=1,N_accepted

                   node_1_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(1,i_move))

                   node_2_between=between(reflect_moves%nodes(:,ref_moves_accepted(j_move)),&
                        reflect_moves%dirs(:,ref_moves_accepted(j_move)),&
                        reflect_moves%nodes(2,i_move))

                   !test for agreement with the betweenness
                   if (node_1_between .neqv. node_2_between) then
                      ref_cycle=.true.
                      exit
                   end if

                end do

                !cycle the while loop if incompatible
                if (ref_cycle .eqv. .true.) then

                   N_rejected=N_rejected+1
                   count=count+1
                   moves_rejected(N_rejected)=i_move
                   call reduce_remaining_moves(N_moves,i_move,&
                        moves_remaining,N_remaining)
                   cycle

                end if
                   
             end if

             
             !write(unit_log,*)'     [move_ring] calculating reflect move'
             new_reflect_move=reflect_move(in_ring,reflect_moves%nodes(:,i_move),&
                  reflect_moves%dirs(:,i_move))
             !write(unit_log,*)'     [move_ring] executing reflect move'
             call perform_move_obstacles_MH(move_success,in_ring,new_reflect_move,&
                  in_obstacle_distribution,in_energy_params)

             !add rejected move to reject list and increment count
             if (move_success .eqv. .false.) then
                !write(unit_log,*)'   [move_ring] rejecting sample move'
                N_rejected=N_rejected+1
                count=count+1
                moves_rejected(N_rejected)=i_move
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
             elseif (N_accepted .lt. multiplicity) then
                N_accepted=N_accepted+1
                ref_moves_accepted(N_accepted)=i_move
                nodes_moved=nodes_moved+new_reflect_move%N
                call reduce_remaining_moves(N_moves,i_move,&
                     moves_remaining,N_remaining)
                move_success=.false.
             end if

          end do !end loop to find valid move

       end if

       if (allocated(kink_moves%nodes)) deallocate(kink_moves%nodes)
       if (allocated(new_kink_move%x)) deallocate(new_kink_move%x)
       if (allocated(new_kink_move%node_indices)) deallocate(new_kink_move%node_indices)
       if (allocated(reflect_moves%nodes)) deallocate(reflect_moves%nodes)
       if (allocated(reflect_moves%dirs)) deallocate(reflect_moves%dirs)
       if (allocated(new_reflect_move%x)) deallocate(new_reflect_move%x)
       if (allocated(new_reflect_move%node_indices)) deallocate(new_reflect_move%node_indices)
       if (allocated(moves_rejected)) deallocate(moves_rejected)
       if (allocated(moves_remaining)) deallocate(moves_remaining)

       !run additional error tests
       if ((mod(i_iter,in_deform_params%err_check_freq) .eq. 0) .or. &
            (i_iter .eq. in_deform_params%num_iter)) then

          length_error=.true.
          intersect_error=.true.
          obstacle_error=.true.
          
          !calculate the length
          test_length=ring_length(in_ring)

          !test if length is the same
          if (test_length .ne. initial_length) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- initial ring length = ',initial_length
             write(unit_log,"(a,i5.5,a,i6)")'  [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- current ring length = ',test_length
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring length inconsistency'
             exit
          else
             length_error=.false.
          end if

          !test for intersections
          config_test=valid_config(in_ring)

          if (sum(config_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- exiting routine due to self-intersecting ring configuration'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- interesection at nodes',&
                  config_test(1),' and ',config_test(2)
             exit
          else
             intersect_error=.false.
          end if

          !test for obstacle intersections
          node_obst_test=obstacle_test_fast(in_ring,in_obstacle_distribution)
          !node_obst_test=obstacle_test(in_ring,in_obstacle_distribution)
          
          if (sum(node_obst_test) .ne. -2) then
             write(unit_log,"(a,i5.5,a,i8)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' ---ERROR at i_iter = ',i_iter
             write(unit_log,"(a,i5.5,a)")' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
                  ' --- exiting routine due to ring-obstacle intersection'
             write(unit_log,"(a,i5.5,a,i6,a,i6)")'  [move_ring_obstacles_MH] --- i_rep=',i_rep,&
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

             !create a new set of deform params
             restart_deform_params=in_deform_params
             restart_deform_params%num_iter=restart_deform_params%num_iter-i_iter

             !adjust the deform params for the iterations completed so far
             if (i_iter .gt. restart_deform_params%ref_wait) then
                restart_deform_params%ref_wait=0
                if (i_iter .gt. restart_deform_params%ref_init+&
                     restart_deform_params%ref_wait) then
                   restart_deform_params%ref_init=0
                else
                   restart_deform_params%ref_init=restart_deform_params%ref_init+&
                        restart_deform_params%ref_wait-i_iter
                end if
             else
                restart_deform_params%ref_wait=restart_deform_params%ref_wait-i_iter
             end if

             !write the restart file
             call write_restart_file(in_ring,in_obstacle_distribution,&
                  restart_deform_params,in_grow_params,unit_rep,.false.)

             close(unit_rep)

             if (track_energy .eqv. .true.) then
                format_string="(a,i5.5,a)"
                !write(unit_log,format_string)'  [move_ring_MH] --- i_rep=',i_rep,' --- calculate energy'
                call total_energy(temp_energy,in_ring,in_energy_params)
                open(UNIT=unit_rep,FILE=trim(energy_file),ACTION='write',STATUS='unknown',POSITION='append')
                write(unit_rep,"(i9,a,es20.10)")nodes_moved,',',temp_energy
                close(unit_rep)
             end if

             if (i_iter .eq. in_deform_params%num_iter) then
                sim_completed=.true.
             end if
             
          end if
          
       end if !end conditional for error checking
       
       if (mod(i_iter,write_update_freq) .eq. 0) then
          format_string="(a,i5.5,a,i8,a,i8,a)"
          write(unit_log,format_string)' [move_ring_obstacles_MH] --- i_rep=',i_rep,&
               ' --- progress: ',&
               i_iter,' of ',in_deform_params%num_iter,' iterations completed'
       end if
       
    end do !end main loop

    if (track_energy .eqv. .true.) then
       call total_energy(temp_energy,in_ring,in_energy_params)
       write(unit_log,"(a,i5.5,a,es20.10,a,i12)")'[move_ring_obstacles_MH] --- i_rep=',i_rep,&
            ' --- END, energy = ',temp_energy,', nodes_moved = ',nodes_moved
    else
       write(unit_log,"(a,i5.5,a)")'[move_ring_obstacles_MH] --- i_rep=',i_rep,' --- END'
    end if
    
  end subroutine move_ring_obstacles_MH
  
end module ring_polymer
