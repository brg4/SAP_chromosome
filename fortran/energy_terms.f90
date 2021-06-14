module energy_terms

  !this module contains the the derived type used for the energy terms and
  !subroutines and functions for evaluating the potentials

  use params
  use ring_objects
  use ring_procedures

  implicit none

  !derived type used to hold energy parameters
  type :: energy_params
     integer :: V_nn_type, N_nn_params
     integer :: V_bend_type, N_bend_params
     integer :: V_twist_type, N_twist_params
     integer :: V_coord_type, N_coord_params, N_V_coord
     integer, allocatable :: species_coord(:,:)
     integer :: V_obst_type, N_obst_params, N_V_obst
     integer, allocatable :: species_obst(:,:)
     integer :: V_spectwist_type, N_spectwist_params, N_spectwist
     integer, allocatable :: spectwist_active(:,:)
     integer :: V_pairwise_type, N_pair_params, N_species
     integer, allocatable :: pairs_active(:,:,:)
     real(rp) :: beta
     real(rp), allocatable :: nn_params(:), bend_params(:), twist_params(:)
     real(rp), allocatable :: coord_params(:,:), obst_params(:,:), spectwist_params(:,:)
     real(rp), allocatable :: pair_params(:,:,:)
  end type energy_params
  
contains

  !function to allocate new energy parameters
  function new_energy_params(beta,&
       V_nn_type,N_nn_params,&
       V_bend_type,N_bend_params,&
       V_twist_type,N_twist_params,&
       V_coord_type,N_coord_params,N_V_coord,&
       V_obst_type,N_obst_params,N_V_obst,&
       V_spectwist_type,N_spectwist_params,N_spectwist,&
       V_pairwise_type,N_pair_params,N_species)
    
    implicit none
    
    integer, intent(in) :: V_nn_type, V_bend_type, V_twist_type
    integer, intent(in) :: V_coord_type, V_spectwist_type, V_pairwise_type, V_obst_type
    integer, intent(in) :: N_nn_params, N_bend_params, N_twist_params
    integer, intent(in) :: N_coord_params, N_V_coord
    integer, intent(in) :: N_obst_params, N_V_obst
    integer, intent(in) :: N_spectwist_params, N_spectwist
    integer, intent(in) :: N_pair_params, N_species
    real(rp), intent(in) :: beta
    type(energy_params) :: new_energy_params


    new_energy_params%beta=beta

    new_energy_params%V_nn_type=V_nn_type
    new_energy_params%V_bend_type=V_bend_type
    new_energy_params%V_twist_type=V_twist_type
    new_energy_params%V_coord_type=V_coord_type
    new_energy_params%V_obst_type=V_obst_type
    new_energy_params%V_spectwist_type=V_spectwist_type
    new_energy_params%V_pairwise_type=V_pairwise_type

    new_energy_params%N_nn_params=N_nn_params
    new_energy_params%N_bend_params=N_bend_params
    new_energy_params%N_twist_params=N_twist_params
    
    new_energy_params%N_V_coord=N_V_coord
    new_energy_params%N_coord_params=N_coord_params
    
    new_energy_params%N_V_obst=N_V_obst
    new_energy_params%N_obst_params=N_obst_params

    new_energy_params%N_species=N_species
    new_energy_params%N_pair_params=N_pair_params

    !deallocate existing species interaction arrays and reallocate them
    if (allocated(new_energy_params%species_coord))&
         deallocate(new_energy_params%species_coord)
    if (allocated(new_energy_params%species_obst))&
         deallocate(new_energy_params%species_obst)
    if (allocated(new_energy_params%spectwist_active))&
         deallocate(new_energy_params%spectwist_active)
    if (allocated(new_energy_params%pairs_active))&
         deallocate(new_energy_params%pairs_active)

    allocate(new_energy_params%species_coord(1:N_species,1:N_V_coord))
    allocate(new_energy_params%species_obst(1:N_species,1:N_V_obst))
    allocate(new_energy_params%spectwist_active(1:2,1:N_species))
    allocate(new_energy_params%pairs_active(1:2,1:N_species,1:N_species))

    !deallocate existing parameter arrays and reallocate them
    if (allocated(new_energy_params%nn_params))&
         deallocate(new_energy_params%nn_params)
    if (allocated(new_energy_params%bend_params))&
         deallocate(new_energy_params%bend_params)
    if (allocated(new_energy_params%twist_params))&
         deallocate(new_energy_params%twist_params)
    if (allocated(new_energy_params%coord_params))&
         deallocate(new_energy_params%coord_params)
    if (allocated(new_energy_params%obst_params))&
         deallocate(new_energy_params%obst_params)
    if (allocated(new_energy_params%spectwist_params))&
         deallocate(new_energy_params%spectwist_params)
    if (allocated(new_energy_params%pair_params))&
         deallocate(new_energy_params%pair_params)

    allocate(new_energy_params%nn_params(1:N_nn_params))
    allocate(new_energy_params%bend_params(1:N_bend_params))
    allocate(new_energy_params%twist_params(1:N_twist_params))
    allocate(new_energy_params%coord_params(1:N_coord_params,1:N_V_coord))
    allocate(new_energy_params%obst_params(1:N_obst_params,1:N_V_obst))
    allocate(new_energy_params%spectwist_params(1:N_spectwist_params,1:N_species))
    allocate(new_energy_params%pair_params(1:N_pair_params,1:N_species,1:N_species))


    !initialize species interaction arrays to zero
    new_energy_params%species_coord=0
    new_energy_params%species_obst=0
    new_energy_params%spectwist_active=0
    new_energy_params%pairs_active=0

    !initial parameter arrays to zero
    new_energy_params%nn_params=0.0d0
    new_energy_params%bend_params=0.0d0
    new_energy_params%twist_params=0.0d0
    new_energy_params%coord_params=0.0d0
    new_energy_params%obst_params=0.0d0
    new_energy_params%spectwist_params=0.0d0
    new_energy_params%pair_params=0.0d0
    
    
  end function new_energy_params
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! POTENTIAL ROUTINES AND FUNCTIONS !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  pure function V_pairwise(dim,sep,i_spec,j_spec,x_i,x_j,in_energy_params)

    implicit none
    
    integer, intent(in) :: sep, i_spec, j_spec, dim
    integer(ip), intent(in) :: x_i(1:dim), x_j(1:dim)
    type(energy_params), intent(in) :: in_energy_params
    
    real(rp) :: V_pairwise, d_sqrd

    V_pairwise=0.0d0

    ! if ((in_energy_params%V_nn_type .eq. 1) .and.&
    !      (sep .eq. 1)) then
    !    V_pairwise=V_pairwise+in_energy_params%nn_params(1)
    ! end if
    
    if (in_energy_params%pairs_active(1,i_spec,j_spec) .eq. 1) then
       if (sep .eq. 1) then
          V_pairwise=in_energy_params%pair_params(1,&
               i_spec,j_spec)
       end if
    elseif (in_energy_params%pairs_active(1,i_spec,j_spec) .eq. 2) then
       d_sqrd=sum((x_i-x_j)**2.0d0)
       V_pairwise=in_energy_params%pair_params(1,i_spec,j_spec)*d_sqrd
    elseif (in_energy_params%pairs_active(1,i_spec,j_spec) .eq. 3) then
       d_sqrd=sum((x_i-x_j)**2.0d0)
       V_pairwise=in_energy_params%pair_params(1,i_spec,j_spec)*&
            (sqrt(d_sqrd)-in_energy_params%pair_params(2,i_spec,j_spec))**2.0d0
    end if

  end function V_pairwise

  pure function V_coord(dim,specie,x,in_ring,in_energy_params)

    implicit none
    
    integer, intent(in) :: specie, dim
    integer(ip), intent(in) :: x(1:dim)
    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_V
    
    real(rp) :: V_coord

    V_coord=0.0d0
    
    if (in_energy_params%V_pairwise_type .eq. 1) then
       do i_V=1,in_energy_params%N_V_coord
          if (in_energy_params%species_coord(specie,i_V) .eq. 1) then
             V_coord=V_coord+in_energy_params%coord_params(1,i_V)*&
                  sum((x-in_energy_params%coord_params(2:in_ring%dim+1,i_V))**2.0d0)
          end if
       end do
    end if

  end function V_coord

  pure function V_bend(dim,dx_pair,in_energy_params)

    implicit none

    integer, intent(in) :: dim
    integer(ip), intent(in) :: dx_pair(1:dim,1:2)
    
    type(energy_params), intent(in) :: in_energy_params
    
    real(rp) :: V_bend

    V_bend=-in_energy_params%bend_params(1)*dot_product(dx_pair(:,1),dx_pair(:,2))
    
  end function V_bend

  pure function V_pairwise_total(in_ring,in_energy_params)

    implicit none

    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, j_node, sep

    real(rp) :: V_pairwise_total

    V_pairwise_total=0.0d0
    

    if (in_energy_params%V_pairwise_type .ne. 0) then
       
       do i_node=1,in_ring%N_nodes-1
          do j_node=i_node+1,in_ring%N_nodes
             sep=sum(abs(in_ring%x(:,i_node)-in_ring%x(:,j_node)))
             V_pairwise_total=V_pairwise_total+&
                  V_pairwise(in_ring%dim,sep,&
                  in_ring%species(i_node),in_ring%species(j_node),&
                  in_ring%x(:,i_node),in_ring%x(:,i_node),&
                  in_energy_params)
          end do
       end do
       
    end if
    
  end function V_pairwise_total

  pure function V_nn_pairwise_total(in_ring,in_energy_params)

    implicit none

    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, j_node, sep, total_sep

    real(rp) :: V_nn_total, V_pairwise_total, V_nn_pairwise_total

    V_nn_total=0.0d0
    V_pairwise_total=0.0d0
    V_nn_pairwise_total=0.0d0

    total_sep=0
    
       
    do i_node=1,in_ring%N_nodes-1
       do j_node=i_node+1,in_ring%N_nodes
          sep=sum(abs(in_ring%x(:,i_node)-in_ring%x(:,j_node)))
          if (sep .eq. 1) total_sep=total_sep+sep
          if (in_energy_params%V_pairwise_type .ne. 0) then
             V_pairwise_total=V_pairwise_total+&
                  V_pairwise(in_ring%dim,sep,&
                  in_ring%species(i_node),in_ring%species(j_node),&
                  in_ring%x(:,i_node),in_ring%x(:,i_node),&
                  in_energy_params)
          end if
       end do
    end do

    if (in_energy_params%V_nn_type .ne. 0) then
       V_nn_total=total_sep*in_energy_params%nn_params(1)
    end if

    V_nn_pairwise_total=V_nn_total+V_pairwise_total

  end function V_nn_pairwise_total


  pure function V_bend_total(in_ring,in_energy_params)

    implicit none

    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    integer :: alignment(1:in_ring%N_nodes)

    real(rp) :: V_bend_total

    alignment=sum((cshift(in_ring%x,1,2)-in_ring%x)*&
         (in_ring%x-cshift(in_ring%x,-1,2)),1)

    V_bend_total=-in_energy_params%bend_params(1)*sum(alignment)

  end function V_bend_total


  pure function V_twist_total(in_ring,in_energy_params)

    implicit none

    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    integer :: i_node, i_quad
    integer :: temp_twist
    integer :: quad(1:4)
    integer :: p_tens(-3:3,-3:3)
    integer(ip) :: dx_triplet(in_ring%dim,3)

    real(rp) :: temp_param
    real(rp) :: V_twist_total

    V_twist_total=0.0d0

    p_tens=perm_tensor()
    
    do i_node=1,in_ring%N_nodes
       do i_quad=1,4
          quad(i_quad)=count_transform_alt(i_quad-2,i_node,in_ring%N_nodes)
       end do
       temp_param=calc_twist_param(in_ring%N_nodes,in_ring%species,&
            in_energy_params,quad)
       do i_quad=1,3
          dx_triplet(:,i_quad)=in_ring%x(:,quad(i_quad+1))-in_ring%x(:,quad(i_quad))
       end do
       call dihedral(dx_triplet,p_tens,temp_twist)
       V_twist_total=V_twist_total+temp_param*temp_twist
    end do
    
  end function V_twist_total

  subroutine minimize_ring_V_pairwise(in_ring,in_energy_params)

    implicit none

    type(ring), intent(inout) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    real(rp) :: V_min, V_eval
    integer :: i_shift, min_shift, shift_increment
    integer(ip) :: temp_coords(1:in_ring%dim,1:in_ring%N_nodes)

    V_min=0.0d0
    min_shift=0

    temp_coords=in_ring%x

    V_min=V_pairwise_total(in_ring,in_energy_params)

    shift_increment=ceiling(in_ring%N_nodes/1.0d2)
    
    !do i_shift=1,in_ring%N_nodes-1
    do i_shift=1,in_ring%N_nodes-1,shift_increment
       in_ring%x=cshift(temp_coords,i_shift,2)
       V_eval=V_pairwise_total(in_ring,in_energy_params)
       if (V_eval .lt. V_min) then
          V_min=V_eval
          min_shift=i_shift
       end if
    end do

    if (min_shift .ne. 0) then
       in_ring%x=cshift(temp_coords,min_shift,2)
    else
       in_ring%x=temp_coords
    end if
    
    
  end subroutine minimize_ring_V_pairwise

  pure function calc_twist_param(N,spec_list,in_energy_params,quad)

    implicit none

    integer, intent(in) :: N
    integer, intent(in) :: spec_list(1:N)
    type(energy_params), intent(in) :: in_energy_params
    integer, intent(in) :: quad(1:4)

    integer :: i, temp_spec
    
    real(rp) :: calc_twist_param

    calc_twist_param=0.0d0

    do i=1,4
       temp_spec=spec_list(quad(i))
       if (in_energy_params%spectwist_active(1,temp_spec) .eq. 1) then
          calc_twist_param=calc_twist_param+in_energy_params%spectwist_params(1,temp_spec)
       end if
    end do

    calc_twist_param=calc_twist_param/4.0d0

    
  end function calc_twist_param

  pure subroutine total_energy(energy,in_ring,in_energy_params)

    implicit none

    type(ring), intent(in) :: in_ring
    type(energy_params), intent(in) :: in_energy_params

    real(rp), intent(out) :: energy

    energy=0.0d0

    if ((in_energy_params%V_nn_type .ne. 0) .or.&
         (in_energy_params%V_pairwise_type .ne. 0)) then
       energy=energy+V_nn_pairwise_total(in_ring,in_energy_params)
    end if
    if (in_energy_params%V_bend_type .ne. 0) then
       energy=energy+V_bend_total(in_ring,in_energy_params)
    end if
    if (in_energy_params%V_twist_type .ne. 0) then
       energy=energy+V_twist_total(in_ring,in_energy_params)
    end if

  end subroutine

end module energy_terms
