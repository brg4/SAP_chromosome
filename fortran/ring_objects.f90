module ring_objects

  !this module contains the derived types used to model the ring polymer and
  !subroutines and functions for managing those derived types
  
  use params
  use octree
  
  implicit none

  !derived type used to hold the parameters describing the deformation
  type deform_params
     integer :: ergodic_flag
     integer :: ref_init, ref_wait, ref_freq, ref_mult
     integer :: kink_mult, kink_repeat
     integer :: num_iter, max_iter, err_check_freq
     integer :: obst_freq, obst_mult
     integer :: min_rep, max_rep, N_reps
     character(len=120) :: energy_file
  end type deform_params

  !derived type used to hold the parameters describing the growth
  type grow_params
     integer :: grow_flag
     integer :: growth_target, initial_growth, max_growth, rand_init, ref_node
     integer :: max_iter, err_check_freq
     integer :: min_iter_move, max_wait
     integer :: N_species_final
     integer, allocatable :: species_final(:)
     real(rp) :: alpha(1:2), eta(1:2)
     character(len=120) :: energy_file, species_file
  end type grow_params

  !derived type describing the ring and its constraints
  type ring
     logical :: fixed
     integer :: dim, bound_type
     integer :: N_nodes, N_fixed, N_bounds, N_species
     integer(ip), allocatable :: x(:,:)
     integer, allocatable :: fixed_nodes(:), species(:), bounds(:,:)
  end type ring

  !derived type used to hold coordinate information
  type coords
     integer :: dim, N
     integer(ip), allocatable :: x(:,:)
  end type coords

  !extension of the coords derived type used to describe a subset of the ring
  type, extends(coords) :: ring_subset
     integer, allocatable :: node_indices(:)
  end type ring_subset

  type, extends(coords) :: obstacle_base
     integer :: N_placed
     integer, allocatable :: placed(:), obst_types(:)
  end type obstacle_base

  !extension of the coords derived type used to describe a distribution of obstacles
  type, extends(coords) :: obstacle_distribution
     integer :: N_placed, obst_type
     integer, allocatable :: placed(:), mapping(:)
     type(obstacle_base) :: obst_base
     type(tree_type) :: obst_tree
  end type obstacle_distribution

  !derived type used to hold subsets of nodes
  !this is used for the enumeration of kink and reflect moves
  type node_subsets
     integer :: N, M, N_max
     integer, allocatable :: nodes(:,:)
  end type node_subsets

  !extension of the node_subsets derived type used to describe candidate reflect moves
  type, extends(node_subsets) :: direction_pairs
     integer, allocatable :: dirs(:,:)
  end type direction_pairs


contains

  !function to allocate new ring
  function new_ring(dim,N_nodes,N_fixed,fixed,N_species,bound_type,N_bounds)
    implicit none
    integer, intent(in) :: dim, bound_type
    integer, intent(in) :: N_nodes, N_fixed, N_species, N_bounds
    logical, intent(in) :: fixed
    type(ring) :: new_ring
    new_ring%dim=dim
    new_ring%N_nodes=N_nodes
    new_ring%fixed=fixed
    new_ring%N_fixed=N_fixed
    new_ring%N_species=N_species
    new_ring%bound_type=bound_type
    new_ring%N_bounds=N_bounds
    if(allocated(new_ring%x)) deallocate(new_ring%x)
    if(allocated(new_ring%fixed_nodes)) deallocate(new_ring%fixed_nodes)
    if(allocated(new_ring%species)) deallocate(new_ring%species)
    if(allocated(new_ring%bounds)) deallocate(new_ring%bounds)
    allocate(new_ring%x(1:dim,1:N_nodes))
    allocate(new_ring%species(1:N_nodes))
    allocate(new_ring%fixed_nodes(1:N_fixed))
    if (bound_type .eq. 0) then
       allocate(new_ring%bounds(1:1,1:1))
    elseif (bound_type .eq. 1) then
       allocate(new_ring%bounds(0:2,1:N_bounds))
    elseif (bound_type .eq. 2) then
       allocate(new_ring%bounds(0:dim,1:1))
    end if
    new_ring%x=0
    new_ring%fixed_nodes=0
    new_ring%species=0
    new_ring%bounds=0
  end function new_ring

  !function to allocate new coords
  function new_coords(dim,N)
    implicit none
    integer, intent(in) :: dim, N
    type(coords) :: new_coords
    new_coords%dim=dim
    new_coords%N=N
    if(allocated(new_coords%x)) deallocate(new_coords%x)
    allocate(new_coords%x(1:dim,1:N))
    new_coords%x=0
  end function new_coords

  !function to allocate new ring subset
  function new_ring_subset(dim,N)
    implicit none
    integer, intent(in) :: dim, N
    type(ring_subset) :: new_ring_subset
    new_ring_subset%dim=dim
    new_ring_subset%N=N
    if(allocated(new_ring_subset%x)) deallocate(new_ring_subset%x)
    if(allocated(new_ring_subset%node_indices)) deallocate(new_ring_subset%node_indices)
    allocate(new_ring_subset%x(1:dim,1:N))
    allocate(new_ring_subset%node_indices(1:N))
    new_ring_subset%x=0
    new_ring_subset%node_indices=0
  end function new_ring_subset

  !function to allocate new node subsets
  function new_node_subsets(N,M)
    implicit none
    type(node_subsets) :: new_node_subsets
    integer, intent(in) :: N, M
    new_node_subsets%N=N
    new_node_subsets%N_max=N
    new_node_subsets%M=M
    if(allocated(new_node_subsets%nodes)) deallocate(new_node_subsets%nodes)
    allocate(new_node_subsets%nodes(1:M,1:N))
    new_node_subsets%nodes=0
  end function new_node_subsets
  
  !function to allocate new direction pairs
  function new_direction_pairs(N,M)
    implicit none
    type(direction_pairs) :: new_direction_pairs
    integer, intent(in) :: N, M
    new_direction_pairs%N=N
    new_direction_pairs%N_max=N
    new_direction_pairs%M=M
    if(allocated(new_direction_pairs%nodes)) deallocate(new_direction_pairs%nodes)
    if(allocated(new_direction_pairs%dirs)) deallocate(new_direction_pairs%dirs)
    allocate(new_direction_pairs%nodes(1:M,1:N))
    allocate(new_direction_pairs%dirs(1:2,1:N))
    new_direction_pairs%nodes=0
    new_direction_pairs%dirs=0
  end function new_direction_pairs

  !function to allocate new obstacle distribution
  function new_obstacle_distribution(dim,N,obst_type)
    implicit none
    type(obstacle_distribution) :: new_obstacle_distribution
    integer, intent(in) :: dim, N, obst_type
    new_obstacle_distribution%dim=dim
    new_obstacle_distribution%N=N
    new_obstacle_distribution%obst_type=obst_type
    new_obstacle_distribution%N_placed=0
    if(allocated(new_obstacle_distribution%placed)) deallocate(new_obstacle_distribution%placed)
    if(allocated(new_obstacle_distribution%mapping)) deallocate(new_obstacle_distribution%mapping)
    if(allocated(new_obstacle_distribution%x)) deallocate(new_obstacle_distribution%x)
    allocate(new_obstacle_distribution%placed(1:N))
    allocate(new_obstacle_distribution%mapping(1:N))
    allocate(new_obstacle_distribution%x(1:dim,1:N))
    new_obstacle_distribution%placed=0
    new_obstacle_distribution%mapping=0
    new_obstacle_distribution%x=0
  end function new_obstacle_distribution

  !function to initialize the deformation parameters
  function new_deform_params(energy_file,&
       ergodic_flag,num_iter,&
       ref_init,ref_wait,ref_freq,ref_mult,&
       kink_mult,kink_repeat,&
       obst_freq,obst_mult,&
       max_iter,err_check_freq,&
       min_rep,max_rep)
    implicit none
    type(deform_params) :: new_deform_params
    integer, intent(in) :: num_iter, ergodic_flag, max_iter, err_check_freq
    integer, intent(in) :: ref_init , ref_wait, ref_freq, ref_mult, kink_mult, kink_repeat
    integer, intent(in) :: obst_freq, obst_mult
    integer, intent(in) :: min_rep, max_rep
    character(len=120), intent(in) :: energy_file
    new_deform_params%energy_file=energy_file
    new_deform_params%ergodic_flag=ergodic_flag
    new_deform_params%num_iter=num_iter
    new_deform_params%ref_init=ref_init
    new_deform_params%ref_wait=ref_wait
    new_deform_params%ref_freq=ref_freq
    new_deform_params%ref_mult=ref_mult
    new_deform_params%kink_mult=kink_mult
    new_deform_params%kink_repeat=kink_repeat
    new_deform_params%obst_freq=obst_freq
    new_deform_params%obst_mult=obst_mult
    new_deform_params%max_iter=max_iter
    new_deform_params%err_check_freq=err_check_freq
    new_deform_params%min_rep=min_rep
    new_deform_params%max_rep=max_rep
    new_deform_params%N_reps=max_rep-min_rep+1
  end function new_deform_params

  !function to initialize the deformation parameters
  function new_grow_params(energy_file,species_file,&
       grow_flag,growth_target,max_growth,initial_growth,rand_init,ref_node,&
       alpha,eta,&
       min_iter_move,max_wait,&
       err_check_freq,max_iter)
    implicit none
    type(grow_params) :: new_grow_params
    integer, intent(in) :: grow_flag
    integer, intent(in) :: growth_target, initial_growth, max_growth, rand_init, ref_node
    integer, intent(in) :: min_iter_move, max_wait
    integer, intent(in) :: err_check_freq, max_iter
    real(rp), intent(in) :: alpha(1:2), eta(1:2)
    character(len=120), intent(in) :: energy_file, species_file
    new_grow_params%energy_file=energy_file
    new_grow_params%species_file=species_file
    new_grow_params%grow_flag=grow_flag
    new_grow_params%growth_target=growth_target
    new_grow_params%max_growth=max_growth
    new_grow_params%initial_growth=initial_growth
    new_grow_params%rand_init=rand_init
    new_grow_params%ref_node=ref_node
    new_grow_params%alpha=alpha
    new_grow_params%eta=eta
    new_grow_params%min_iter_move=min_iter_move
    new_grow_params%max_wait=max_wait
    new_grow_params%err_check_freq=err_check_freq
    new_grow_params%max_iter=max_iter
  end function new_grow_params


  !subroutine to add a new subset to a subset list
  !this subroutine attempts to avoid repeated allocations and deallocations by
  !instead increasing the allocations by factors of 2 when needed
  subroutine add_to_subset_list(subset_list,N,N_max,M,new_subset)

    implicit none
    integer, intent(inout), allocatable :: subset_list(:,:)
    integer, intent(inout) :: N, N_max
    integer, intent(in), allocatable :: new_subset(:)
    integer, intent(in) :: M
    integer, allocatable :: temp_subset_list(:,:)

    !if there are no elements, initialize the allocations
    if (N .eq. 0) then
       if (allocated(subset_list)) deallocate(subset_list)
       allocate(subset_list(1:M,1:1))
       subset_list(:,1)=new_subset
       N_max=1
    !if it is allocated, proceed
    elseif (allocated(subset_list)) then
       !if the new element is greater than the maximum index, double the allocation
       if (N+1 .gt. N_max) then
          N_max=2*N_max
          allocate(temp_subset_list(1:M,1:N))
          temp_subset_list=subset_list
          deallocate(subset_list)
          allocate(subset_list(1:M,1:N_max))
          subset_list(:,1:N)=temp_subset_list
          subset_list(:,N+1)=new_subset
          deallocate(temp_subset_list)
       !if the new element is within the current allocation, add the element
       else
          subset_list(:,N+1)=new_subset
       end if
    end if

    !increase the element count
    N=N+1
      
  end subroutine add_to_subset_list

  !subroutine to add a new direction pair to a direction pair list
  !this subroutine attempts to avoid repeated allocations and deallocations by
  !instead increasing the allocations by factors of 2 when needed
  subroutine add_to_direction_pairs_list(pair_list,dir_list,N,N_max,new_pair,new_dir)

    implicit none
    integer, intent(inout), allocatable :: pair_list(:,:), dir_list(:,:)
    integer, intent(inout) :: N, N_max
    integer, intent(in) :: new_pair(1:2), new_dir(1:2)
    integer, allocatable :: temp_pair_list(:,:), temp_dir_list(:,:)

    !if there are no elements, initialize the allocations
    if (N .eq. 0) then
       if (allocated(pair_list)) deallocate(pair_list)
       if (allocated(dir_list)) deallocate(dir_list)
       allocate(pair_list(1:2,1:1))
       allocate(dir_list(1:2,1:1))
       pair_list(:,1)=new_pair
       dir_list(:,1)=new_dir
       N_max=1
    !if it is allocated, proceed
    elseif (allocated(pair_list)) then
       !if the new element is greater than the maximum index, double the allocation
       if (N+1 .gt. N_max) then
          N_max=2*N_max
          allocate(temp_pair_list(1:2,1:N),temp_dir_list(1:2,1:N))
          temp_pair_list=pair_list
          temp_dir_list=dir_list
          deallocate(pair_list,dir_list)
          allocate(pair_list(1:2,1:N_max),dir_list(1:2,1:N_max))
          pair_list(:,1:N)=temp_pair_list
          pair_list(:,N+1)=new_pair
          dir_list(:,1:N)=temp_dir_list
          dir_list(:,N+1)=new_dir
          deallocate(temp_pair_list,temp_dir_list)
       !if the element is within the allocation, add the element
       else
          pair_list(:,N+1)=new_pair
          dir_list(:,N+1)=new_dir
       end if
    end if

    !increase the element count
    N=N+1
      
  end subroutine add_to_direction_pairs_list


  subroutine insert_node_pair(in_ring,x_pair,node_pair)

    implicit none

    type(ring), intent(inout) :: in_ring

    integer(ip), intent(in) :: x_pair(1:in_ring%dim,1:2)
    integer, intent(in) :: node_pair(1:2)

    integer :: i_fixed
    
    integer(ip) :: x_temp(1:in_ring%dim,1:in_ring%N_nodes)
    integer :: species_temp(1:in_ring%N_nodes)
    
    x_temp=in_ring%x
    species_temp=in_ring%species

    deallocate(in_ring%x,in_ring%species)

    in_ring%N_nodes=in_ring%N_nodes+2
    
    allocate(in_ring%x(1:in_ring%dim,1:in_ring%N_nodes))
    allocate(in_ring%species(1:in_ring%N_nodes))

    if (in_ring%fixed .eqv. .true.) then

       !determine the brackets of fixed regions
       
       ! allocate(fixed_node_brackets(1:2,1:1))

       ! i_fixed=1
       ! i_bracket=1
       
       ! fixed_node_brackets(1,i_bracket)=in_ring%fixed_nodes(i_fixed)
       ! fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)
       ! N_bracket=1
       
       ! do while (i_fixed .lt. in_ring%N_fixed)

       !    i_fixed=i_fixed+1

       !    if (fixed_node_brackets(1,i_bracket)+1 .eq. in_ring%fixed_nodes(i_fixed)) then

       !       fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)

       !    elseif (i_fixed .lt. in_ring%N_fixed) then

       !       allocate(fixed_node_brackets_temp(1:2,1:N_bracket))

       !       fixed_node_brackets_temp=fixed_node_brackets

       !       deallocate(fixed_node_brackets)
             
       !       N_bracket=N_bracket+1
       !       i_bracket=i_bracket+1

       !       allocate(fixed_node_brackets(1:2,1:N_bracket))

       !       fixed_node_brackets(:,1:N_bracket-1)=fixed_node_brackets_temp

       !       deallocate(fixed_node_brackets_temp)

       !       fixed_node_brackets(1,i_bracket)=in_ring%fixed_nodes(i_fixed)
       !       fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)
             
       !    end if
          
       ! end do

       ! !determine which brackets border the region in which the growth will occur

       ! do i_bracket=1,N_bracket

          
          
       ! end do
       
       !correct the fixed node indices

       do i_fixed=1,in_ring%N_fixed
          if (in_ring%fixed_nodes(i_fixed) .ge. node_pair(2)) then
             in_ring%fixed_nodes(i_fixed)=in_ring%fixed_nodes(i_fixed)+2
          end if
       end do
       
   end if

   !apply the growth
   
   if ((node_pair(1) .eq. 1) .and.&
        (node_pair(2) .eq. in_ring%N_nodes-2)) then

      in_ring%x(:,2:in_ring%N_nodes-1)=x_temp
      in_ring%species(2:in_ring%N_nodes-1)=species_temp

      in_ring%species(1)=species_temp(1)
      in_ring%species(in_ring%N_nodes)=species_temp(in_ring%N_nodes-2)

      in_ring%x(:,1)=x_pair(:,1)
      in_ring%x(:,in_ring%N_nodes)=x_pair(:,2)

   else
      
      if (node_pair(2) .eq. 2) then

         in_ring%x(:,4:in_ring%N_nodes)=x_temp(:,2:in_ring%N_nodes-2)
         in_ring%x(:,1)=x_temp(:,1)
         
         in_ring%species(4:in_ring%N_nodes)=species_temp(2:in_ring%N_nodes-2)
         in_ring%species(1)=species_temp(1)

      elseif (node_pair(1) .eq. in_ring%N_nodes-3) then

         in_ring%x(:,1:in_ring%N_nodes-3)=x_temp(:,1:in_ring%N_nodes-3)
         in_ring%x(:,in_ring%N_nodes)=x_temp(:,in_ring%N_nodes-2)

         in_ring%species(1:in_ring%N_nodes-3)=species_temp(1:in_ring%N_nodes-3)
         in_ring%species(in_ring%N_nodes)=species_temp(in_ring%N_nodes-2)

      else

         in_ring%x(:,1:node_pair(1))=x_temp(:,1:node_pair(1))
         in_ring%x(:,node_pair(1)+3:in_ring%N_nodes)=x_temp(:,node_pair(2):in_ring%N_nodes-2)
         
         in_ring%species(1:node_pair(1))=&
              species_temp(1:node_pair(1))
         in_ring%species(node_pair(1)+3:in_ring%N_nodes)=&
              species_temp(node_pair(2):in_ring%N_nodes-2)

      end if

      in_ring%species(node_pair(1)+1:node_pair(2)+1)=species_temp(node_pair(1):node_pair(2))
      !in_ring%x(:,node_pair(1)+1:node_pair(2)+1)=x_temp(:,node_pair(1):node_pair(2))
      in_ring%x(:,node_pair(1)+1:node_pair(2)+1)=x_pair
      
   end if

 end subroutine insert_node_pair

 !subroutine to insert multiple node pairs
 subroutine insert_node_pairs(in_ring,N,x_pairs,node_pair)

   implicit none

   type(ring), intent(inout) :: in_ring

   integer, intent(in) :: N
   integer(ip), intent(in) :: x_pairs(1:in_ring%dim,1:N)
   integer, intent(in) :: node_pair(1:2)

   integer :: i_fixed, i_node

   integer(ip) :: x_temp(1:in_ring%dim,1:in_ring%N_nodes)
   integer :: species_temp(1:in_ring%N_nodes)


   x_temp=in_ring%x
   species_temp=in_ring%species

   deallocate(in_ring%x,in_ring%species)

   in_ring%N_nodes=in_ring%N_nodes+N

   allocate(in_ring%x(1:in_ring%dim,1:in_ring%N_nodes))
   allocate(in_ring%species(1:in_ring%N_nodes))

   if (in_ring%fixed .eqv. .true.) then

      !determine the brackets of fixed regions

      ! allocate(fixed_node_brackets(1:2,1:1))

      ! i_fixed=1
      ! i_bracket=1

      ! fixed_node_brackets(1,i_bracket)=in_ring%fixed_nodes(i_fixed)
      ! fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)
      ! N_bracket=1

      ! do while (i_fixed .lt. in_ring%N_fixed)

      !    i_fixed=i_fixed+1

      !    if (fixed_node_brackets(1,i_bracket)+1 .eq. in_ring%fixed_nodes(i_fixed)) then

      !       fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)

      !    elseif (i_fixed .lt. in_ring%N_fixed) then

      !       allocate(fixed_node_brackets_temp(1:2,1:N_bracket))

      !       fixed_node_brackets_temp=fixed_node_brackets

      !       deallocate(fixed_node_brackets)

      !       N_bracket=N_bracket+1
      !       i_bracket=i_bracket+1

      !       allocate(fixed_node_brackets(1:2,1:N_bracket))

      !       fixed_node_brackets(:,1:N_bracket-1)=fixed_node_brackets_temp

      !       deallocate(fixed_node_brackets_temp)

      !       fixed_node_brackets(1,i_bracket)=in_ring%fixed_nodes(i_fixed)
      !       fixed_node_brackets(2,i_bracket)=in_ring%fixed_nodes(i_fixed)

      !    end if

      ! end do

      ! !determine which brackets border the region in which the growth will occur

      ! do i_bracket=1,N_bracket



      ! end do

      !correct the fixed node indices

      do i_fixed=1,in_ring%N_fixed
         if (in_ring%fixed_nodes(i_fixed) .ge. node_pair(2)) then
            in_ring%fixed_nodes(i_fixed)=in_ring%fixed_nodes(i_fixed)+N
         end if
      end do

   end if

   !apply the growth
   ! print *,'N_current=',in_ring%N_nodes-N
   ! print *,'node_pair'
   ! print *,node_pair
   ! print *,'x'
   ! print *,x_temp(:,node_pair(1))
   ! print *,x_temp(:,node_pair(2))
   ! print *,'x_pairs'
   ! do i_fixed=1,N
   !    print *,x_pairs(:,i_fixed)
   ! end do

   if ((node_pair(1) .eq. 1) .and.&
        (node_pair(2) .eq. in_ring%N_nodes-N)) then

      in_ring%x(:,N+1:in_ring%N_nodes)=x_temp
      do i_node=1,N
         in_ring%x(:,i_node)=x_pairs(:,N-i_node+1)
      end do
      in_ring%species(1:N/2)=species_temp(node_pair(2))
      in_ring%species(N/2+1:N)=species_temp(node_pair(1))
      ! in_ring%x(:,N/2+1:in_ring%N_nodes-N/2)=x_temp
      ! in_ring%species(N/2+1:in_ring%N_nodes-N/2)=species_temp

      ! in_ring%species(1:N/2)=species_temp(1)
      ! in_ring%species(in_ring%N_nodes-N/2+1:in_ring%N_nodes)=species_temp(in_ring%N_nodes-N)

      ! in_ring%x(:,1:N/2)=x_pairs(:,1:N/2)
      ! in_ring%x(:,in_ring%N_nodes-N/2+1:in_ring%N_nodes)=x_pairs(:,N/2+1:N)

   else

      in_ring%x(:,1:node_pair(1))=x_temp(:,1:node_pair(1))
      in_ring%x(:,node_pair(1)+1+N:in_ring%N_nodes)=x_temp(:,node_pair(2):in_ring%N_nodes-N)

      in_ring%species(1:node_pair(1))=&
           species_temp(1:node_pair(1))
      in_ring%species(node_pair(1)+1+N:in_ring%N_nodes)=&
           species_temp(node_pair(2):in_ring%N_nodes-N)

      in_ring%species(node_pair(1)+1:node_pair(1)+N/2)=species_temp(node_pair(1))
      in_ring%species(node_pair(1)+N/2+1:node_pair(1)+N)=species_temp(node_pair(2))
      !in_ring%x(:,node_pair(1)+1:node_pair(2)+1)=x_temp(:,node_pair(1):node_pair(2))
      in_ring%x(:,node_pair(1)+1:node_pair(1)+N)=x_pairs

   end if

 end subroutine insert_node_pairs
 
 !subroutine to convert obstacle types to an octree
 subroutine convert_obstacle_distribution(in_obstacle_distribution)

   implicit none
   
   type(obstacle_distribution), intent(inout) :: in_obstacle_distribution

   integer :: i_obst, i_placed, i_dim
   integer :: count, sites_per_type
   integer :: big_shift, shift
   integer :: small_dx(1:in_obstacle_distribution%dim), big_dx(1:in_obstacle_distribution%dim)
   integer :: i_dx, j_dx, k_dx

   write(unit_log,*)"[convert_obstacle_distribution] BEGIN"

   if (allocated(in_obstacle_distribution%obst_base%x)) deallocate(in_obstacle_distribution%obst_base%x)
   if (allocated(in_obstacle_distribution%obst_base%placed)) deallocate(in_obstacle_distribution%obst_base%placed)
   if (allocated(in_obstacle_distribution%obst_base%obst_types)) deallocate(in_obstacle_distribution%obst_base%obst_types)

   in_obstacle_distribution%obst_base%N=in_obstacle_distribution%N
   in_obstacle_distribution%obst_base%dim=in_obstacle_distribution%dim
   in_obstacle_distribution%obst_base%N_placed=in_obstacle_distribution%N_placed

   allocate(in_obstacle_distribution%obst_base%x(1:in_obstacle_distribution%obst_base%dim,&
        1:in_obstacle_distribution%obst_base%N))
   allocate(in_obstacle_distribution%obst_base%placed(1:in_obstacle_distribution%obst_base%N))
   allocate(in_obstacle_distribution%obst_base%obst_types(1:in_obstacle_distribution%obst_base%N))


   in_obstacle_distribution%obst_base%x=in_obstacle_distribution%x
   in_obstacle_distribution%obst_base%placed=0
   in_obstacle_distribution%obst_base%obst_types=in_obstacle_distribution%mapping
   do i_placed=1,in_obstacle_distribution%N_placed
      in_obstacle_distribution%obst_base%placed(in_obstacle_distribution%placed(i_placed))=1
   end do
   

   deallocate(in_obstacle_distribution%x)
   deallocate(in_obstacle_distribution%placed)
   deallocate(in_obstacle_distribution%mapping)


   !override all obstacles types to match what is given by obst_type
   if (in_obstacle_distribution%obst_type .ne. -1) then
      in_obstacle_distribution%obst_base%obst_types=in_obstacle_distribution%obst_type
   end if
      

   !determine the correct size of the converted obstacle arrays
   in_obstacle_distribution%N=0
   in_obstacle_distribution%N_placed=0
   do i_obst=1,in_obstacle_distribution%obst_base%N
      
      if ((in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 0) .or.&
           (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 1)) then
         sites_per_type=1
      elseif (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 2) then
         shift=1
         sites_per_type=1+in_obstacle_distribution%dim*2
      elseif (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 3) then
         shift=2
         sites_per_type=(shift**in_obstacle_distribution%dim)*(1+in_obstacle_distribution%dim*2)
      end if

      in_obstacle_distribution%N=in_obstacle_distribution%N+sites_per_type
      if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
         in_obstacle_distribution%N_placed=in_obstacle_distribution%N_placed+sites_per_type
      end if
      
   end do

   !allocate the obstacle arrays
   allocate(in_obstacle_distribution%x(1:in_obstacle_distribution%dim,1:in_obstacle_distribution%N))
   allocate(in_obstacle_distribution%placed(1:in_obstacle_distribution%N))
   allocate(in_obstacle_distribution%mapping(1:in_obstacle_distribution%N))


   count=1
   
   !loop over the set of base obstacles with the types given by obst_types and fill the obstacle arrays
   do i_obst=1,in_obstacle_distribution%obst_base%N

      if ((in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 0) .or.&
           (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 1)) then
         
         sites_per_type=1

         in_obstacle_distribution%x(:,count)=&
              in_obstacle_distribution%obst_base%x(:,i_obst)
         in_obstacle_distribution%placed(count)=&
              in_obstacle_distribution%obst_base%placed(i_obst)
         in_obstacle_distribution%mapping(count)=i_obst

         count=count+1
         
         
      elseif (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 2) then
         
         shift=1
         sites_per_type=1+in_obstacle_distribution%dim*2

         in_obstacle_distribution%x(:,count)=in_obstacle_distribution%obst_base%x(:,i_obst)
         in_obstacle_distribution%placed(count)=&
              in_obstacle_distribution%obst_base%placed(i_obst)
         in_obstacle_distribution%mapping(count)=i_obst

         count=count+1

         do i_dim=1,in_obstacle_distribution%dim

            big_dx=0
            
            do big_shift=-shift,shift,2*shift

               big_dx(i_dim)=big_shift

               in_obstacle_distribution%x(:,count)=in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx
               in_obstacle_distribution%placed(count)=&
                    in_obstacle_distribution%obst_base%placed(i_obst)
               in_obstacle_distribution%mapping(count)=i_obst
               
               count=count+1
               
            end do
            
         end do
         
      elseif (in_obstacle_distribution%obst_base%obst_types(i_obst) .eq. 3) then
         
         shift=2
         sites_per_type=(shift**in_obstacle_distribution%dim)*(1+in_obstacle_distribution%dim*2)

         !system is 2D
         if (in_obstacle_distribution%dim .eq. 2) then

            big_dx=0

            do i_dx=0,(shift-1)
               do j_dx=0,(shift-1)

                  small_dx=(/i_dx,j_dx/)

                  in_obstacle_distribution%x(:,count)=&
                       in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
                  in_obstacle_distribution%placed(count)=&
                       in_obstacle_distribution%obst_base%placed(i_obst)
                  in_obstacle_distribution%mapping(count)=i_obst

                  count=count+1

               end do !end loop over y shifts
            end do !end loop over y shifts

            do i_dim=1,in_obstacle_distribution%dim

               big_dx=0

               do big_shift=-shift,shift,2*shift

                  big_dx=0
                  big_dx(i_dim)=big_shift

                  do i_dx=0,(shift-1)
                     do j_dx=0,(shift-1)

                        small_dx=(/i_dx,j_dx/)

                        in_obstacle_distribution%x(:,count)=&
                             in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
                        in_obstacle_distribution%placed(count)=&
                             in_obstacle_distribution%obst_base%placed(i_obst)
                        in_obstacle_distribution%mapping(count)=i_obst

                        count=count+1

                     end do !end loop over y shifts
                  end do !end loop over x shifts

               end do !end loop over big shifts

            end do !end loop over dimensions

            !system is 3D
         elseif (in_obstacle_distribution%dim .eq. 3) then

            big_dx=0

            do i_dx=0,(shift-1)
               do j_dx=0,(shift-1)
                  do k_dx=0,(shift-1)

                     small_dx=(/i_dx,j_dx,k_dx/)

                     in_obstacle_distribution%x(:,count)=&
                          in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
                     in_obstacle_distribution%placed(count)=&
                          in_obstacle_distribution%obst_base%placed(i_obst)
                     in_obstacle_distribution%mapping(count)=i_obst                  

                     count=count+1

                  end do !end loop over z shifts
               end do !end loop over y shifts
            end do !end loop over y shifts

            do i_dim=1,in_obstacle_distribution%dim

               do big_shift=-shift,shift,2*shift

                  big_dx=0
                  big_dx(i_dim)=big_shift

                  do i_dx=0,(shift-1)
                     do j_dx=0,(shift-1)
                        do k_dx=0,(shift-1)

                           small_dx=(/i_dx,j_dx,k_dx/)

                           in_obstacle_distribution%x(:,count)=&
                                in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
                           in_obstacle_distribution%placed(count)=&
                                in_obstacle_distribution%obst_base%placed(i_obst)
                           in_obstacle_distribution%mapping(count)=i_obst
               
                           count=count+1

                        end do !end loop over z shifts
                     end do !end loop over y shifts
                  end do !end loop over x shifts

               end do !end loop over big shifts

            end do !end loop over dimensions

         end if !end condition for dimension of system
         
      end if !end conditional for type of obstacles
      
   end do !end loop over obstacles
   

   !obstacles are points
   ! if ((in_obstacle_distribution%obst_type .eq. 0) .or.&
   !      (in_obstacle_distribution%obst_type .eq. 1)) then

      
   !    ! do nothing
   !    in_obstacle_distribution%N=in_obstacle_distribution%obst_base%N
   !    in_obstacle_distribution%N_placed=in_obstacle_distribution%obst_base%N_placed

   !    allocate(in_obstacle_distribution%x(1:in_obstacle_distribution%dim,1:in_obstacle_distribution%N))
   !    allocate(in_obstacle_distribution%placed(1:in_obstacle_distribution%N))
   !    allocate(in_obstacle_distribution%mapping(1:in_obstacle_distribution%N))

   !    in_obstacle_distribution%x=in_obstacle_distribution%obst_base%x
   !    in_obstacle_distribution%placed=in_obstacle_distribution%obst_base%placed
   !    in_obstacle_distribution%mapping=(/(i, i=1, in_obstacle_distribution%N, 1)/)
      

   ! !obstacles are stars
   ! elseif (in_obstacle_distribution%obst_type .eq. 2) then

   !    shift=1
      
   !    sites_per_type=1+in_obstacle_distribution%dim*2

   !    in_obstacle_distribution%N=in_obstacle_distribution%obst_base%N*sites_per_type
   !    in_obstacle_distribution%N_placed=in_obstacle_distribution%obst_base%N_placed*sites_per_type

   !    allocate(in_obstacle_distribution%x(1:in_obstacle_distribution%dim,1:in_obstacle_distribution%N))
   !    allocate(in_obstacle_distribution%placed(1:in_obstacle_distribution%N))
   !    allocate(in_obstacle_distribution%mapping(1:in_obstacle_distribution%N))

   !    in_obstacle_distribution%placed=0

   !    count=1

   !    do i_obst=1,in_obstacle_distribution%obst_base%N

   !       in_obstacle_distribution%mapping(count)=i_obst
   !       in_obstacle_distribution%x(:,count)=in_obstacle_distribution%obst_base%x(:,i_obst)

   !       if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !          count_placed=count
   !          in_obstacle_distribution%placed(count_placed)=1
   !       end if

   !       count=count+1

   !       do i_dim=1,in_obstacle_distribution%dim

   !          big_dx=0
            
   !          do big_shift=-shift,shift,2*shift

   !             big_dx(i_dim)=big_shift

   !             in_obstacle_distribution%x(:,count)=in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx

   !             if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !                count_placed=count
   !                in_obstacle_distribution%placed(count_placed)=1
   !             end if
               
   !             count=count+1
               
   !          end do
            
   !       end do
         
   !    end do

   ! !obstacles are 2-site stars      
   ! elseif (in_obstacle_distribution%obst_type .eq. 3) then

   !    shift=2

   !    sites_per_type=(shift**in_obstacle_distribution%dim)*(1+in_obstacle_distribution%dim*2)

   !    in_obstacle_distribution%N=in_obstacle_distribution%obst_base%N*sites_per_type
   !    in_obstacle_distribution%N_placed=in_obstacle_distribution%obst_base%N_placed*sites_per_type

   !    allocate(in_obstacle_distribution%x(1:in_obstacle_distribution%dim,1:in_obstacle_distribution%N*sites_per_type))
   !    allocate(in_obstacle_distribution%placed(1:in_obstacle_distribution%N*sites_per_type))
   !    allocate(in_obstacle_distribution%mapping(1:in_obstacle_distribution%N*sites_per_type))

   !    in_obstacle_distribution%placed=0

   !    count=1
   !    count_placed=1

   !    !system is 2D
   !    if (in_obstacle_distribution%dim .eq. 2) then

   !       do i_obst=1,in_obstacle_distribution%obst_base%N

   !          big_dx=0

   !          do i_dx=0,(shift-1)
   !             do j_dx=0,(shift-1)
                  
   !                small_dx=(/i_dx,j_dx/)
                  
   !                in_obstacle_distribution%x(:,count)=&
   !                     in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
   !                in_obstacle_distribution%mapping(count)=i_obst

   !                if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !                   count_placed=count
   !                   in_obstacle_distribution%placed(count_placed)=1
   !                end if

   !                count=count+1
                  
   !             end do !end loop over y shifts
   !          end do !end loop over y shifts

   !          do i_dim=1,in_obstacle_distribution%dim
               
   !             big_dx=0

   !             do big_shift=-shift,shift,2*shift
                  
   !                big_dx=0
   !                big_dx(i_dim)=big_shift

   !                do i_dx=0,(shift-1)
   !                   do j_dx=0,(shift-1)
                        
   !                      small_dx=(/i_dx,j_dx/)

   !                      in_obstacle_distribution%x(:,count)=&
   !                           in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
   !                      in_obstacle_distribution%mapping(count)=i_obst

   !                      if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !                         count_placed=count
   !                         in_obstacle_distribution%placed(count_placed)=1
   !                      end if

   !                      count=count+1

   !                   end do !end loop over y shifts
   !                end do !end loop over x shifts
                  
   !             end do !end loop over big shifts

   !          end do !end loop over dimensions
   
   !       end do !end loop over nodes

   !    !system is 3D
   !    elseif (in_obstacle_distribution%dim .eq. 3) then

   !       do i_obst=1,in_obstacle_distribution%obst_base%N

   !          big_dx=0

   !          do i_dx=0,(shift-1)
   !             do j_dx=0,(shift-1)
   !                do k_dx=0,(shift-1)

   !                   small_dx=(/i_dx,j_dx,k_dx/)

   !                   in_obstacle_distribution%x(:,count)=&
   !                        in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
   !                   in_obstacle_distribution%mapping(count)=i_obst

   !                   if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !                      count_placed=count
   !                      in_obstacle_distribution%placed(count_placed)=1
   !                   end if

   !                   count=count+1

   !                end do !end loop over z shifts
   !             end do !end loop over y shifts
   !          end do !end loop over y shifts

   !          do i_dim=1,in_obstacle_distribution%dim
               

   !             do big_shift=-shift,shift,2*shift
                  
   !                big_dx=0
   !                big_dx(i_dim)=big_shift

   !                do i_dx=0,(shift-1)
   !                   do j_dx=0,(shift-1)
   !                      do k_dx=0,(shift-1)

   !                         small_dx=(/i_dx,j_dx,k_dx/)

   !                         in_obstacle_distribution%x(:,count)=&
   !                              in_obstacle_distribution%obst_base%x(:,i_obst)+big_dx+small_dx
   !                         in_obstacle_distribution%mapping(count)=i_obst

   !                         if (in_obstacle_distribution%obst_base%placed(i_obst) .eq. 1) then
   !                            count_placed=count
   !                            in_obstacle_distribution%placed(count_placed)=1
   !                         end if

   !                         count=count+1

   !                      end do !end loop over z shifts
   !                   end do !end loop over y shifts
   !                end do !end loop over x shifts
                  
   !             end do !end loop over big shifts

   !          end do !end loop over dimensions
   
   !       end do !end loop over nodes

   !    end if !end condition for dimension of system
      
   ! end if !end conditional for obstacle types

   ! do i_obst=1,in_obstacle_distribution%N
   !    print *,i_obst,in_obstacle_distribution%x(:,i_obst),in_obstacle_distribution%mapping(i_obst)
   ! end do

   write(unit_log,*)" [convert_obstacle_distribution] obst_base_N = ",in_obstacle_distribution%obst_base%N
   write(unit_log,*)" [convert_obstacle_distribution] obst_base_N_placed = ",in_obstacle_distribution%obst_base%N_placed
   write(unit_log,*)" [convert_obstacle_distribution] N = ",in_obstacle_distribution%N
   write(unit_log,*)" [convert_obstacle_distribution] N_placed = ",in_obstacle_distribution%N_placed

   write(unit_log,*)"[convert_obstacle_distribution] END"
   
 end subroutine convert_obstacle_distribution

 !subroutine to copy obstacle distributions
 subroutine copy_obst_dist(out_obstacle_distribution,in_obstacle_distribution)

   implicit none

   type(obstacle_distribution), intent(in) :: in_obstacle_distribution
   type(obstacle_distribution), intent(out) :: out_obstacle_distribution

   logical :: valid_tree

   out_obstacle_distribution=new_obstacle_distribution(in_obstacle_distribution%dim,&
        in_obstacle_distribution%N,in_obstacle_distribution%obst_type)
   out_obstacle_distribution%N_placed=in_obstacle_distribution%N_placed
   out_obstacle_distribution%x=in_obstacle_distribution%x
   out_obstacle_distribution%placed=in_obstacle_distribution%placed
   out_obstacle_distribution%mapping=in_obstacle_distribution%mapping

   if (allocated(out_obstacle_distribution%obst_base%x)) deallocate(out_obstacle_distribution%obst_base%x)
   if (allocated(out_obstacle_distribution%obst_base%placed)) deallocate(out_obstacle_distribution%obst_base%placed)
   if (allocated(out_obstacle_distribution%obst_base%obst_types)) deallocate(out_obstacle_distribution%obst_base%obst_types)
   
   out_obstacle_distribution%obst_base%dim=in_obstacle_distribution%obst_base%dim
   out_obstacle_distribution%obst_base%N=in_obstacle_distribution%obst_base%N
   out_obstacle_distribution%obst_base%N_placed=in_obstacle_distribution%obst_base%N_placed

   allocate(out_obstacle_distribution%obst_base%x(1:out_obstacle_distribution%obst_base%dim,&
        1:out_obstacle_distribution%obst_base%N))
   allocate(out_obstacle_distribution%obst_base%placed(1:out_obstacle_distribution%obst_base%N))
   allocate(out_obstacle_distribution%obst_base%obst_types(1:out_obstacle_distribution%obst_base%N))

   out_obstacle_distribution%obst_base%x=in_obstacle_distribution%obst_base%x
   out_obstacle_distribution%obst_base%placed=in_obstacle_distribution%obst_base%placed
   out_obstacle_distribution%obst_base%obst_types=in_obstacle_distribution%obst_base%obst_types

   call copy_tree(out_obstacle_distribution%obst_tree,in_obstacle_distribution%obst_tree)

   call validate_tree(valid_tree,out_obstacle_distribution)

   if (valid_tree .eqv. .false.) then
      stop
   end if
   
 end subroutine copy_obst_dist

 !subroutine to validate octree
 subroutine validate_tree(valid,in_obstacle_distribution)

   implicit none

   logical, intent(out) :: valid
   type(obstacle_distribution), intent(in) :: in_obstacle_distribution

   logical :: found_flag
   integer :: i_obst, point

   found_flag=.false.

   do i_obst=1,in_obstacle_distribution%N
      if (in_obstacle_distribution%placed(i_obst) .eq. 1) then
         call octree_search(found_flag,point,&
              in_obstacle_distribution%x(:,i_obst),&
              in_obstacle_distribution%N,&
              in_obstacle_distribution%x,&
              in_obstacle_distribution%obst_tree%root_node)
         ! if (i_obst .ne. point) then
         !    print *, in_obstacle_distribution%x(:,i_obst)
         !    print *, in_obstacle_distribution%x(:,point)
         ! end if
      end if
      if (found_flag .eqv. .false.) then
         valid=.false.
         return
      end if
   end do

   valid=.true.

 end subroutine validate_tree

end module ring_objects
