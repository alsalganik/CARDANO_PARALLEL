module task
  use mpi
  implicit none
contains
  subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
    implicit none
    real(8), intent(in), dimension(:,:) :: A
    integer(4), intent(out) :: x1, y1, x2, y2
    integer(4) :: n
    integer(4) :: L, R, Up, Down
    integer(4) :: m, tmp
    real(8), allocatable :: currcolumn(:), B(:,:)
    real(8) :: currsum, maxsum, global_maxsum
    logical :: istrans
    integer(4) :: mpiErr, mpiSize, mpiRank, tmp_mpiRank, global_mpiRank

    call mpi_comm_size(mpi_comm_world, mpiSize, mpiErr)

    call mpi_comm_rank(mpi_comm_world, mpiRank, mpiErr)

    m = size(A, dim=1)
    n = size(A, dim=2)
    istrans = .FALSE.

    if (m < n) then
       istrans = .TRUE.
       B = transpose(A)
       m = size(B, 1)
       n = size(B, 2)
    else
       B = A
    endif

    allocate(currcolumn(m))

    maxsum=B(1, 1)
    x1=1
    y1=1
    x2=1
    y2=1

    do L = mpiRank + 1, n, mpiSize
       currcolumn = B(:, L)
       do R=L, n

          if (R > L) then
             currcolumn = currcolumn + B(:, R)
          endif

          call FindMaxInArray(currcolumn, currsum, Up, Down)

          if (currsum > maxsum) then
             maxsum = currsum
             x1 = Up
             x2 = Down
             y1 = L
             y2 = R
          endif

       end do
    end do

    call mpi_reduce(maxsum, global_maxsum, 1, mpi_real8, mpi_max, 0, mpi_comm_world, mpiErr)

    call mpi_bcast(global_maxsum, 1, mpi_real8, 0, mpi_comm_world, mpiErr)

    tmp_mpiRank = -1
    if (maxsum .eq. global_maxsum) then
       tmp_mpiRank = mpiRank
    endif

    call mpi_allreduce(tmp_mpiRank, global_mpiRank, 1, mpi_integer4, mpi_max, mpi_comm_world, mpiErr)

    call mpi_bcast(x1, 1, mpi_real8, global_mpiRank, mpi_comm_world, mpiErr)
    call mpi_bcast(x2, 1, mpi_real8, global_mpiRank, mpi_comm_world, mpiErr)
    call mpi_bcast(y1, 1, mpi_real8, global_mpiRank, mpi_comm_world, mpiErr)
    call mpi_bcast(y2, 1, mpi_real8, global_mpiRank, mpi_comm_world, mpiErr)

    deallocate(currcolumn)

    if (istrans) then
       tmp = x1
       x1 = y1
       y1 = tmp

       tmp = y2
       y2 = x2
       x2 = tmp
    endif

  end subroutine GetMaxCoordinates


  subroutine FindMaxInArray(a, Sum, Up, Down)
    real(8), intent(in), dimension(:) :: A
    integer(4), intent(out) :: Up, Down
    real(8), intent(out) :: Sum
    real(8) :: cur_sum
    integer(4) :: minus_pos, i

    Sum = a(1)
    Up = 1
    Down = 1
    cur_sum = 0
    minus_pos = 0

    do i=1, size(A)
       cur_sum = cur_sum + A(i)
       if (cur_sum > Sum) then
          Sum = cur_sum
          Up = minus_pos + 1
          Down = i
       endif

       if (cur_sum < 0) then
          cur_sum = 0
          minus_pos = i
       endif

    enddo
  end subroutine FindMaxInArray

end module task
