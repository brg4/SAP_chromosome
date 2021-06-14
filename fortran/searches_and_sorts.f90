module searches_and_sorts

  use params
  
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SEARCHES AND SORTS !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !routine to perform a binary search for the indices in an array
  !that bracket the target, this is used to search for repeated samplings
  pure subroutine binary_search(lo,hi,A,N,tar)

    implicit none

    integer, intent(in) :: N, tar
    integer, intent(in) :: A(1:N)
    integer, intent(inout) :: lo, hi

    integer :: mid
    logical :: target_found, find

    if (N .gt. 1) then

       !select midpoint
       lo=1
       hi=N

       target_found=.false.

       do while((target_found .eqv. .false.) .and. (lo .lt. hi))

          mid=lo+(hi-lo)/2

          if (A(mid) .eq. tar) then
             target_found=.true.

             !target is above the midpoint
          elseif (A(mid) .lt. tar) then
             lo=mid+1
             !R=N

             !target is below the midpoint
          elseif (A(mid) .gt. tar) then
             !lo=1
             if (mid .ne. 1) then
                hi=mid-1
             else
                hi=1
             end if

          end if
       end do

       if ((lo .eq. hi) .and. (A(lo) .eq. tar)) then
          mid=lo
          target_found=.true.
       end if

       !target has been found, now search for full interval with target value
       if (target_found .eqv. .true.) then

          lo=mid
          hi=mid

          !search for the lower bound of the interval containing the target value
          find=.true.
          do while(find .eqv. .true.)
             if (lo .gt. 1) then
                if (A(lo-1) .eq. tar) then
                   lo=lo-1
                else
                   find=.false.
                end if
             else
                find=.false.
             end if
          end do

          !search for the upper bound of the interval containing the target value
          find=.true.
          do while(find .eqv. .true.)
             if (hi .lt. N) then
                if (A(hi+1) .eq. tar) then
                   hi=hi+1
                else
                   find=.false.
                end if
             else
                find=.false.
             end if
          end do

       else !if the target was not found, return -1's
          lo=-1
          hi=-1
       end if

    else

       !check if only element is the target
       if (A(1) .eq. tar) then
          lo=1
          hi=1
       else !if the target was not found, return -1's
          lo=-1
          hi=-1
       end if

    end if !end conditional checking if target was found


  end subroutine binary_search

  !routine to perform a shell sort of array A with N elements
  !this sort type was chosen because it seems to have the optimal performance
  !for a nearly sorted array
  !edit by BRG - this should potentially be replaced with a radix sort
  pure subroutine shell_sort(A,N)

    implicit none

    integer, intent(in) :: N
    integer, intent(inout) :: A(1:N)

    integer :: i, j, inc
    logical :: loop_continue
    integer :: temp

    !initial increment is half of the array size
    inc=N/2

    !while the increment is nonzero sort
    do while (inc .ge. 1)

       !iterate over the upper elements based on the increment
       do i=(inc+1),N
          j=i
          temp=A(i)
          !swap the elements
          loop_continue=.false.
          if (j .ge. inc+1) then
             if (A(j-inc) .gt. temp) then
                loop_continue=.true.
             end if
          end if
          do while (loop_continue .eqv. .true.)
             A(j)=A(j-inc)
             j=j-inc
             if (j .ge. inc+1) then
                if (A(j-inc) .le. temp) then
                   loop_continue=.false.
                end if
             else
                loop_continue=.false.
             end if
          end do
          A(j)=temp
       end do

       !reduce the increment
       if (inc .eq. 2) then
          inc=1
       else
          inc=inc*5/11
       end if

    end do !finish iterating over the increments

  end subroutine shell_sort


end module searches_and_sorts
