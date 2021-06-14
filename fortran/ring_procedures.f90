module ring_procedures

  !this module contains the subroutines and functions used as helper procedures
  !in the main routines to move the ring and evaluate the energies

  use params
  use ring_objects

  implicit none

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! HELPER ROUTINES AND FUNCTIONS !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !function to generate rotation matrices for the lattice
  function rot_mat(dim,x1,x2,type_flag)

    implicit none
    
    integer, intent(in) :: dim, type_flag
    integer(ip), intent(in) :: x1(1:dim), x2(1:dim)
    integer :: theta, c, s
    real(rp) :: rot_mat(1:dim,1:dim), dx(1:dim)
    
    dx=x1-x2
    dx=dx/sqrt(sum(dx**2.0d0))

    if (dim .eq. 2) then
       rot_mat(1,1)=2*dx(1)*dx(1)-1
       rot_mat(2,1)=2*dx(1)*dx(2)
       rot_mat(1,2)=2*dx(2)*dx(1)
       rot_mat(2,2)=2*dx(2)*dx(2)-1
    else
       
       if (type_flag .eq. 0) then
          theta=2
       else
          theta=int_rand(1,3)
       end if

       if (theta .eq. 1) then
          c=0
          s=1
       elseif (theta .eq. 2) then
          c=-1
          s=0
       elseif (theta .eq. 3) then
          c=0
          s=-1
       end if

       rot_mat(1,1)=dx(1)*dx(1)*(1-c)+c
       rot_mat(2,1)=dx(1)*dx(2)*(1-c)+dx(3)*s
       rot_mat(3,1)=dx(1)*dx(3)*(1-c)-dx(2)*s
       rot_mat(1,2)=dx(2)*dx(1)*(1-c)-dx(3)*s
       rot_mat(2,2)=dx(2)*dx(2)*(1-c)+c
       rot_mat(3,2)=dx(2)*dx(3)*(1-c)+dx(1)*s
       rot_mat(1,3)=dx(3)*dx(1)*(1-c)+dx(2)*s
       rot_mat(2,3)=dx(3)*dx(2)*(1-c)-dx(1)*s
       rot_mat(3,3)=dx(3)*dx(3)*(1-c)+c
       
    end if
    
  end function rot_mat

  !function to generate an integer on the interval [lb,ub]
  function int_rand(lb,ub)
    
    implicit none
    
    integer, intent(in) :: lb, ub
    real(rp) :: u
    integer :: int_rand

    call random_number(u)
    int_rand=lb+floor((ub+1-lb)*u)

  end function int_rand

  subroutine generate_permutation(a,n)

    implicit none

    integer, intent(in) :: n
    integer, intent(inout) :: a(1:n)

    integer :: i, u, temp

    do i=n,2,-1
       u=int_rand(1,i)
       if (u .ne. i) then
          temp=a(u)
          a(u)=a(i)
          a(i)=temp
       end if
    end do

  end subroutine generate_permutation

  !function to calculate the minimum distance and directionality between two nodes on the ring
  pure function min_distance(N_nodes,j_node,k_node)

    implicit none
    
    integer, intent(in) :: N_nodes, j_node, k_node
    integer :: half, min_distance(1:2)

    half=N_nodes/2

    if (k_node .ge. j_node) then
       if (k_node-j_node .ge. half) then
          min_distance(1)=j_node+n_nodes-k_node
          min_distance(2)=-1
       else
          min_distance(1)=k_node-j_node
          min_distance(2)=1
       end if
    else
       if (j_node-k_node .ge. half) then
          min_distance(1)=k_node+n_nodes-j_node
          min_distance(2)=1
       else
          min_distance(1)=j_node-k_node
          min_distance(2)=-1
       end if
    end if
    
  end function min_distance

  !function to test if node is between direction pair
  pure function between(pair,direction,test_node)

    implicit none

    integer, intent(in) :: test_node
    integer, intent(in) :: pair(1:2), direction(1:2)
    logical :: between

    between=.false.
    
    if (pair(1) .lt. pair(2)) then

       if (direction(2) .eq. 1) then
          if ((test_node .gt. pair(1)) .and. (test_node .lt. pair(2))) then
             between=.true.
          end if
       else
          if ((test_node .gt. pair(2)) .or. (test_node .lt. pair(1))) then
             between=.true.
          end if
       end if
       
    else

       if (direction(2) .eq. 1) then
          if ((test_node .gt. pair(1)) .or. (test_node .lt. pair(2))) then
             between=.true.
          end if
       else
          if ((test_node .gt. pair(2)) .and. (test_node .lt. pair(1))) then
             between=.true.
          end if
       end if

    end if

  end function between

  !function to test for fixed nodes between direction pair
  pure function fixed_between(N_fixed,fixed_nodes,pair,direction)

    implicit none

    integer, intent(in) :: N_fixed
    integer, intent(in) :: fixed_nodes(1:N_fixed), pair(1:2), direction(1:2)
    logical :: fixed_between

    integer :: i_node
    
    fixed_between=.false.

    do i_node=1,N_fixed
       fixed_between=between(pair,direction,fixed_nodes(i_node))
       if (fixed_between .eqv. .true.) then
          exit
       end if
    end do
       
  end function fixed_between

  !function to transform a displacement along the ring to the proper index
  pure function count_transform(count,i_node,N_nodes)

    implicit none
    
    integer, intent(in) :: count, i_node, N_nodes
    integer :: count_transform

    count_transform=count+i_node
    
    if (count_transform .gt. N_nodes) then
       count_transform=count_transform-N_nodes
    elseif (count_transform .lt. 1) then
       count_transform=count_transform+N_nodes
    end if
    
  end function count_transform

    !function to transform a displacement along the ring to the proper index
  pure function count_transform_alt(count,i_node,N_nodes)

    implicit none
    
    integer, intent(in) :: count, i_node, N_nodes
    integer :: count_transform_alt

    count_transform_alt=count+i_node
    
    if (count_transform_alt .gt. N_nodes) then
       count_transform_alt=count_transform_alt-N_nodes
    elseif (count_transform_alt .lt. 1) then
       count_transform_alt=count_transform_alt+N_nodes
    end if
    
  end function count_transform_alt

  pure subroutine update_species(ref_node,N,N_final,final_species,new_species)

   implicit none

   integer, intent(in) :: ref_node, N, N_final
   integer, intent(in) :: final_species(1:N_final)
   integer, intent(out) :: new_species(1:N)

   integer :: i_node, temp_node

   temp_node=ref_node-N/2+1
   if (temp_node .lt. 1) then
      temp_node=temp_node+N_final
   end if
   
   do i_node=1,N
      new_species(i_node)=final_species(count_transform(i_node,temp_node,N_final))
   end do
   
   
 end subroutine update_species

 subroutine propose_shift(in_ring,x0)

   implicit none

   type(ring), intent(inout) :: in_ring
   integer(ip), intent(in) :: x0(1:in_ring%dim,in_ring%N_nodes)

   integer :: i_node, i_dim, rand_shift

   integer(ip) :: x_base(1:in_ring%dim,in_ring%N_nodes), x_com(1:in_ring%dim)

   x_com=sum(x0,dim=2)/in_ring%N_nodes

   !write(unit_log,*)'called propose_shift'
   
   do i_node=1,in_ring%N_nodes
      x_base(:,i_node)=x0(:,i_node)-x_com
   end do

   if (in_ring%bound_type .eq. 2) then
      do i_dim=1,in_ring%dim
         rand_shift=int_rand(in_ring%bounds(i_dim,1)-in_ring%bounds(0,1),&
              in_ring%bounds(i_dim,1)+in_ring%bounds(0,1))
         !write(unit_log,*)rand_shift
         x_base(i_dim,:)=x_base(i_dim,:)+rand_shift
              
      end do
   end if

   in_ring%x=x_base

 end subroutine propose_shift
 

 !function to calculate permutation tensor
 pure function perm_tensor()

   implicit none

   integer :: perm_tensor(-3:3,-3:3)

   perm_tensor=0

   perm_tensor(-3,-2)=-1
   perm_tensor(-3,-1)=2
   perm_tensor(-2,-1)=-3

   perm_tensor(-3,2)=1
   perm_tensor(-3,1)=-2
   perm_tensor(-2,1)=3

   perm_tensor(-2,3)=-1
   perm_tensor(-1,3)=2
   perm_tensor(-1,2)=-3

   perm_tensor(2,3)=1
   perm_tensor(1,3)=-2
   perm_tensor(1,2)=3

   perm_tensor=perm_tensor-transpose(perm_tensor)
   
 end function perm_tensor

 !function to calculate cross product
 pure function cross_product(p_tens,v1,v2)

   implicit none

   integer, intent(in) :: p_tens(-3:3,-3:3), v1, v2

   integer :: cross_product

   cross_product=p_tens(v1,v2)
   
 end function cross_product

 !function to convert to basis vector
 pure function e_basis(v)

   implicit none

   integer, intent(in) :: v(1:3)

   integer :: i_dim
   integer :: e_basis

   e_basis=0

   do i_dim=1,3
      if (v(i_dim) .eq. 1) then
         e_basis=i_dim
      elseif(v(i_dim) .eq. -1) then
         e_basis=-i_dim
      end if
   end do

 end function e_basis

 !subroutine to calculate twist of segment
 pure subroutine dihedral(dx,p_tens,twist)

   implicit none

   integer(ip), intent(in) :: dx(1:3,1:3)
   integer, intent(in) :: p_tens(-3:3,-3:3)
   integer, intent(out) :: twist

   integer :: v1, v2, e_arr(1:3)
   integer :: i_dx

   do i_dx=1,3
      e_arr(i_dx)=e_basis(dx(:,i_dx))
   end do

   v1=cross_product(p_tens,e_arr(1),e_arr(2))
   v2=cross_product(p_tens,e_arr(2),e_arr(3))

   !print *,"v1=",v1
   !print *,"v2=",v2

   twist=cross_product(p_tens,v1,v2)

   !print *,"twist=",twist

   if (twist .gt. 0) then
      twist=1
   elseif (twist .lt. 0) then
      twist=-1
   end if

 end subroutine dihedral

 !subroutine to calculate total twist
 pure subroutine total_twist(in_ring,tot_twist)

   implicit none

   type(ring), intent(in) :: in_ring

   integer, intent(out) :: tot_twist

   integer :: p_tens(-3:3,-3:3)
   integer(ip) :: dx(in_ring%dim,in_ring%N_nodes), temp_dx(1:in_ring%dim,1:3)

   integer :: temp_twist
   integer :: i_dx

   tot_twist=0

   p_tens=perm_tensor()

   dx=in_ring%x-cshift(in_ring%x,shift=-1,dim=2)

   do i_dx=1,in_ring%N_nodes

      if (i_dx .eq. 1) then
         temp_dx(:,1)=in_ring%x(:,in_ring%N_nodes)
         temp_dx(:,2)=in_ring%x(:,i_dx)
         temp_dx(:,3)=in_ring%x(:,i_dx+1)
      elseif (i_dx .eq. in_ring%N_nodes) then
         temp_dx(:,1)=in_ring%x(:,i_dx-1)
         temp_dx(:,2)=in_ring%x(:,i_dx)
         temp_dx(:,3)=in_ring%x(:,1)
      else
         temp_dx(:,1)=in_ring%x(:,i_dx-1)
         temp_dx(:,2)=in_ring%x(:,i_dx)
         temp_dx(:,3)=in_ring%x(:,i_dx+1)
      end if

      call dihedral(temp_dx,p_tens,temp_twist)

      tot_twist=tot_twist+temp_twist

   end do

 end subroutine total_twist

 pure subroutine reduce_remaining_moves(N_moves,i_move,moves_remaining,N_remaining)

   implicit none

   integer, intent(in) :: N_moves, i_move
   integer, intent(inout) :: moves_remaining(1:N_moves), N_remaining

   if (i_move .eq. 1) then
      moves_remaining(1:N_remaining-1)=moves_remaining(2:N_remaining)
   ! elseif (i_move .eq. N_remaining) then
   !    moves_remaining(1:N_remaining-1)=moves_remaining(1:N_remaining-2)
   elseif (i_move .lt. N_remaining) then
      moves_remaining(i_move:N_remaining-1)=moves_remaining(i_move+1:N_remaining)
   end if

   N_remaining=N_remaining-1
   
 end subroutine reduce_remaining_moves

end module ring_procedures
