module octree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  this implementation of octrees in fortran is based on the following source code  !!!!
!!!!!!!!!!!!!!!!!!!!!!!  https://github.com/dongli/fortran-octree  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use params
  
  implicit none

  type tree_type
     type(node_type), pointer :: root_node
  end type tree_type

  type node_type
     logical :: leaf_flag
     integer :: depth
     integer :: bbox(1:3,1:2)
     integer :: N, N_active
     integer :: points(8), active_points(8)
     type(node_type), pointer :: parent
     type(node_type), pointer :: children(:)
  end type node_type
  
contains

  !initializing the bounding box within which the octree is defined
  pure function init_bbox(x)

    implicit none
    integer(ip), intent(in) :: x(:,:)

    integer :: init_bbox(1:3,1:2)
    integer :: max_power_2(1:3)
    integer :: i_dim, temp

    init_bbox(:,1)=minval(x,2)
    init_bbox(:,2)=maxval(x,2)

    max_power_2=2
    
    do i_dim=1,3
       if (init_bbox(i_dim,2) .gt. -1) then
          temp=2
          do while (temp .lt. init_bbox(i_dim,2)+1)
             temp=2*temp
          end do
          max_power_2(i_dim)=temp
       end if
       if (init_bbox(i_dim,1) .lt. 0) then
          temp=2
          do while (temp .lt. abs(init_bbox(i_dim,1)))
             temp=2*temp
          end do
          if (temp .gt. max_power_2(i_dim)) then
             max_power_2(i_dim)=temp
          end if
       end if
    end do

    if (maxval(init_bbox(:,1)) .ge. 0) then
       init_bbox(:,1)=0
       init_bbox(:,2)=maxval(max_power_2)-1
    end if

    if (minval(init_bbox(:,2)) .le. -1) then
       init_bbox(:,2)=-1
       init_bbox(:,1)=-maxval(max_power_2)
    end if

    if ((maxval(init_bbox(:,1)) .le. -1) .and.&
         (minval(init_bbox(:,2)) .ge. 0)) then
       init_bbox(:,1)=-maxval(max_power_2)
       init_bbox(:,2)=maxval(max_power_2)-1
    end if
    
  end function init_bbox

  !recursive subdividing of the parent bounding box
  pure function subdivide_bbox(bbox)

    implicit none

    integer, intent(in) :: bbox(3,2)

    integer :: subdivide_bbox(3,2,8)
    integer :: i, j, k, i_subbox
    integer :: half(3)

    i_subbox=1

    half=(bbox(:,2)-bbox(:,1)+1)/2

    do i=-1,1,2
       do j=-1,1,2
          do k=-1,1,2
             if (i .lt. 0) then
                subdivide_bbox(1,1,i_subbox)=bbox(1,1)
                subdivide_bbox(1,2,i_subbox)=bbox(1,1)+half(1)-1
             else
                subdivide_bbox(1,1,i_subbox)=bbox(1,2)-half(1)+1
                subdivide_bbox(1,2,i_subbox)=bbox(1,2)
             end if
             if (j .lt. 0) then
                subdivide_bbox(2,1,i_subbox)=bbox(2,1)
                subdivide_bbox(2,2,i_subbox)=bbox(2,1)+half(2)-1
             else
                subdivide_bbox(2,1,i_subbox)=bbox(2,2)-half(2)+1
                subdivide_bbox(2,2,i_subbox)=bbox(2,2)
             end if
             if (k .lt. 0) then
                subdivide_bbox(3,1,i_subbox)=bbox(3,1)
                subdivide_bbox(3,2,i_subbox)=bbox(3,1)+half(3)-1
             else
                subdivide_bbox(3,1,i_subbox)=bbox(3,2)-half(3)+1
                subdivide_bbox(3,2,i_subbox)=bbox(3,2)
             end if

             i_subbox=i_subbox+1
          end do
       end do
    end do
    
  end function

  !initializing the octree
  subroutine octree_init(tree,bbox)

    implicit none

    type(tree_type), intent(inout) :: tree
    integer, intent(in) :: bbox(3,2)

    allocate(tree%root_node)
    if (.not. associated(tree%root_node)) then
       !print *,'allocating root node'
       allocate(tree%root_node)
    end if
    tree%root_node%leaf_flag=.true.
    !print *,'after allocate test'
    !print *,tree%root_node%N
    call reset_node(tree%root_node)
    !print *,'after reset call test'
    tree%root_node%depth=1
    tree%root_node%bbox=bbox
    
  end subroutine octree_init

  !removing the octree
  subroutine octree_final(tree)

    implicit none

    type(tree_type), intent(inout) :: tree

    call clean_node(tree%root_node)
    !write(unit_log,*)'before deallocate in octree_final'
    deallocate(tree%root_node)
    
  end subroutine octree_final

  !resetting the state of a node in the octree
  subroutine reset_node(node)

    implicit none
    
    type(node_type), intent(inout) :: node

    node%N=0
    !print *,'before parent'
    if (associated(node%parent)) then
       nullify(node%parent)
    end if
    !print *,'before children'
    !if (associated(node%children)) then
    ! print *,'node depth=',node%depth
    ! print *,'size of children=',size(node%children)
    if (node%leaf_flag .eqv. .false.) then
    !if (size(node%children) .eq. 8) then
       !print *,'before deallocate'
       deallocate(node%children)
       !print *,'after deallocate'
    end if
    !print *,'before nullify children'
    node%leaf_flag=.true.
    nullify(node%children)
    !print *,'after nullify children'
    
  end subroutine reset_node

  !recursively cleaning a node
  recursive subroutine clean_node(node)

    implicit none

    type(node_type), intent(inout) :: node
    integer :: i_node

    if (node%leaf_flag .eqv. .false.) then
       do i_node=1,8
          call clean_node(node%children(i_node))
       end do
       !write(unit_log,*)'before deallocate children'
       deallocate(node%children)
    end if
    
  end subroutine clean_node

  !creating child nodes from a parent node
  subroutine subdivide_node(node)

    type(node_type), intent(inout), target :: node
    integer :: sub_bbox(3,2,8)
    integer :: i_child

    sub_bbox=subdivide_bbox(node%bbox)

    node%leaf_flag=.false.
    allocate(node%children(1:8))

    !print *,'node depth=',node%depth
    do i_child=1,8
       !print *,'i_child=',i_child
       !print *,'before reset node'
       node%children(i_child)%leaf_flag=.true.
       call reset_node(node%children(i_child))
       !print *,'after reset node'
       node%children(i_child)%depth=node%depth+1
       node%children(i_child)%bbox=sub_bbox(:,:,i_child)
       node%children(i_child)%parent=>node
    end do
    
  end subroutine subdivide_node
  
  !building an octree from coordinate data
  recursive subroutine octree_build(N,id,x,active,node)

    implicit none

    integer, intent(in) :: N
    integer, intent(in) :: id(1:N), active(1:N)
    integer(ip), intent(in) :: x(1:3,1:N)
    type(node_type), intent(inout) :: node

    logical :: point_repeat
    integer :: N_contained, bbox_satisfy
    integer :: i_child, i_point, i_dim, j_point
    integer :: id_contained(1:N), active_contained(1:N)
    integer(ip) :: x_contained(1:3,1:N)

    !print *,"N=",N
    if (N .eq. 0) then
       node%N=0
       node%N_active=0
       node%points=0
       node%active_points=0
       return
    elseif (N .le. 8) then
       node%N=N
       !print *,"node%N=",node%N
       node%N_active=sum(active)
       node%points=0
       node%active_points=0
       node%points(1:N)=id
       node%active_points(1:N)=active
       return
    else
       !print *,'begin subdividing at depth=',node%depth
       call subdivide_node(node)
       !print *,'before loop over children'

       !loop over the children
       do i_child=1,8
          
          N_contained=0
          !loop over the points to find those contained in each child node
          do i_point=1,N
             
             bbox_satisfy=0
             do i_dim=1,3
                if ((x(i_dim,i_point) .ge. node%children(i_child)%bbox(i_dim,1)) .and.&
                     (x(i_dim,i_point) .le. node%children(i_child)%bbox(i_dim,2))) then
                   bbox_satisfy=bbox_satisfy+1
                end if
             end do

             !if the points are within the boundary box continue
             if (bbox_satisfy .eq. 3) then

                
                if (N_contained .gt. 0) then !look for repeats
                   
                   point_repeat=.false.
                   do j_point=1,N_contained
                      if (sum(abs(x_contained(:,j_point)-x(:,i_point))) .eq. 0) then
                         point_repeat=.true.
                         exit
                      end if
                   end do

                   !if the point is not a repeat, then add the point
                   if (point_repeat .eqv. .false.) then
                      N_contained=N_contained+1
                      id_contained(N_contained)=id(i_point)
                      active_contained(N_contained)=active(i_point)
                      x_contained(:,N_contained)=x(:,i_point)
                   end if

                else !if this is the first point
                   N_contained=N_contained+1
                   id_contained(N_contained)=id(i_point)
                   active_contained(N_contained)=active(i_point)
                   x_contained(:,N_contained)=x(:,i_point)
                end if
                
             end if !end conditional for points within boundary box
             
          end do !end looping over points
          
          !print *,'before recurring build'

          !print *,"i_child=",i_child
          !print *,"N_contained=",N_contained
          
          !recursively build octree for the child
          call octree_build(N_contained,&
               id_contained(1:N_contained),&
               x_contained(:,1:N_contained),&
               active_contained(1:N_contained),&
               node%children(i_child))

          ! if (N_contained .gt. 0) then
          !    call octree_build(N_contained,&
          !         id_contained(1:N_contained),&
          !         x_contained(:,1:N_contained),&
          !         active_contained(1:N_contained),&
          !         node%children(i_child))
          ! else
          !    node%children(i_child)%N=0
          !    node%children(i_child)%N_active=0
          ! end if
          
       end do !end loop over children
       !print *,'end subdividing'
    end if
    
  end subroutine octree_build

  !searching an octree for a coordinate
  recursive subroutine octree_search(found_flag,point_found,x_test,N,x,in_node)

    implicit none

    logical, intent(out) :: found_flag
    integer, intent(out) :: point_found
    integer, intent(in) :: N
    integer(ip), intent(in) :: x_test(3)
    integer(ip), intent(in) :: x(3,N)
    type(node_type), intent(in), target :: in_node

    integer :: i_point, i_child, i_dim, bbox_satisfy
    type(node_type), pointer :: node

    found_flag=.false.
    
    node=>in_node

    if (node%leaf_flag .eqv. .false.) then
       do i_child=1,8
          bbox_satisfy=0
          do i_dim=1,3
             if ((x_test(i_dim) .ge. node%children(i_child)%bbox(i_dim,1)) .and.&
                  (x_test(i_dim) .le. node%children(i_child)%bbox(i_dim,2))) then
                bbox_satisfy=bbox_satisfy+1
             end if
          end do
          if (bbox_satisfy .eq. 3) then
             call octree_search(found_flag,point_found,&
                  x_test,N,x,node%children(i_child))
          end if
       end do
    else
       if (node%N_active .eq. 0) then
          return
       else
          !print *,node%active_points
          ! write(unit_log,*)'x_test = ',x_test
          ! write(unit_log,"(A)")'-------leaf---------'
          ! write(unit_log,"(A,6i4)")'bbox: ',node%bbox
          ! write(unit_log,"(A,i4)")'depth: ',node%depth
          ! write(unit_log,"(A,i6)")'N: ',node%N
          ! write(unit_log,*)node%active_points
          ! do i_point=1,node%N
          !    !if (node%active_points(i_point) .eq. 1) then
          !    write(unit_log,*)x(:,node%points(i_point))
          !    !end if
          ! end do
          do i_point=1,node%N
             if (node%active_points(i_point) .eq. 1) then
                if (sum(abs(x_test-x(:,node%points(i_point)))) .eq. 0) then
                   point_found=node%points(i_point)
                   found_flag=.true.
                   return
                end if
             end if
          end do
       end if
    end if
    
  end subroutine octree_search

  !copying an octree
  subroutine copy_tree(out_tree,in_tree)

    type(tree_type), intent(inout) :: out_tree
    type(tree_type), intent(in) :: in_tree

    ! if (associated(out_tree%root_node) .eqv. .true.) then
    !    write(unit_log,*)'before octree_final in copy tree'
    !    write(unit_log,*)out_tree%root_node%leaf_flag
    !    call octree_final(out_tree)
    !    write(unit_log,*)'after octree_final in copy tree'
    ! end if
    !write(unit_log,*)'before octree_init'
    call octree_init(out_tree,in_tree%root_node%bbox)
    call copy_node(out_tree%root_node,in_tree%root_node)

  end subroutine copy_tree

  !copying a node
  recursive subroutine copy_node(out_node,in_node)

    type(node_type), intent(inout) :: out_node
    type(node_type), intent(in) :: in_node

    integer :: i_child
    

    if (in_node%leaf_flag .eqv. .true.) then
       out_node%N=in_node%N
       out_node%N_active=in_node%N_active
       out_node%points=in_node%points
       out_node%active_points=in_node%active_points
       return
    else
       call subdivide_node(out_node)
       do i_child=1,8
          call copy_node(out_node%children(i_child),in_node%children(i_child))
       end do
    end if

  end subroutine copy_node
  
  !printing the contents of a node
  subroutine print_node(node,N,x)

    type(node_type), intent(in) :: node
    integer, intent(in) :: N
    integer(ip), intent(in) :: x(1:3,1:N)

    integer :: i_point

    if (node%N .gt. 0) then
       write(unit_log,"(A)")'-------leaf---------'
       write(unit_log,"(A,6i4)")'bbox: ',node%bbox
       write(unit_log,"(A,i4)")'depth: ',node%depth
       write(unit_log,"(A,i6)")'N: ',node%N
       write(unit_log,*)node%active_points
       do i_point=1,node%N
          !if (node%active_points(i_point) .eq. 1) then
          write(unit_log,*)x(:,node%points(i_point))
          !end if
       end do
       !write(unit_log,*)'Leaf: ',node%leaf_flag
    end if

  end subroutine print_node

  !printing an entire octree
  recursive subroutine print_tree(node,N,x)

    type(node_type), intent(in) :: node
    integer, intent(in) :: N
    integer(ip), intent(in) :: x(1:3,1:N)
    
    integer :: i_child

    if (node%leaf_flag .eqv. .false.) then
       write(unit_log,"(A)")'-------branch---------'
       write(unit_log,"(A,6i4)")'bbox: ',node%bbox
       write(unit_log,"(A,i4)")'depth: ',node%depth
       do i_child=1,8
          call print_tree(node%children(i_child),N,x)
       end do
    else
       call print_node(node,N,x)
    end if

  end subroutine print_tree
  
end module octree
