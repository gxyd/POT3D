module mpi
    implicit none
    integer, parameter :: MPI_COMM_WORLD = 0
    integer, parameter :: MPI_COMM_TYPE_SHARED = 1
    integer, parameter :: MPI_INFO_NULL = 0
    integer, parameter :: MPI_REAL4 = 2
    integer, parameter :: MPI_REAL8 = 3
    integer, parameter :: MPI_THREAD_FUNNELED = 1
    integer, parameter :: MPI_INTEGER = 1
    integer, parameter :: MPI_IN_PLACE = 1
    integer, parameter :: MPI_SUM = 1
    integer, parameter :: MPI_PROC_NULL = -1
    integer, allocatable :: MPI_STATUSES_IGNORE(:)
    integer, parameter :: MPI_STATUS_IGNORE = 0

    interface
        function MPI_Wtime() result(time)
            real(kind=8) :: time
        end function

        subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
            type(*), dimension(..), intent(in) :: sendbuf
            type(*), dimension(..) :: recvbuf
            integer :: count, datatype, op, comm, ierror
        end subroutine

        subroutine MPI_Barrier(comm, ierror)
            integer :: comm, ierror
        end subroutine

        subroutine MPI_Bcast(buffer, count, datatype, root, comm, ierror)
            type(*), dimension(..) :: buffer
            integer, intent(in) :: count, root
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Finalize(ierror)
            integer, optional, intent(in) :: ierror
        end subroutine

        subroutine MPI_Init_thread(required, provided, ierror)
            integer, intent(in) :: required
            integer, intent(out) :: provided
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Comm_size(comm, size, ierror)
            integer, intent(in) :: comm
            integer, intent(out) :: size
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Comm_rank(comm, rank, ierror)
            integer, intent(in) :: comm
            integer, intent(out) :: rank
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Comm_split_type(comm, split_type, key, info, newcomm, ierror)
            integer :: comm
            integer, intent(in) :: split_type, key
            integer, intent(in) :: info
            integer, intent(out) :: newcomm
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
            integer, intent(in) :: comm_old
            integer, intent(in) :: ndims
            integer, intent(in) :: dims(ndims)
            logical, intent(in) :: periods(ndims), reorder
            integer, intent(out) :: comm_cart
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Cart_coords(comm, rank, maxdims, coords, ierror)
            integer, intent(in) :: comm
            integer, intent(in) :: rank, maxdims
            integer, intent(out) :: coords(maxdims)
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest, ierror)
            integer, intent(in) :: comm
            integer, intent(in) :: direction, disp
            integer, intent(out) :: rank_source, rank_dest
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Cart_sub(comm, remain_dims, newcomm, ierror)
            integer, intent(in) :: comm
            logical, intent(in) :: remain_dims(*)
            integer, intent(out) :: newcomm
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
            type(*), dimension(..), intent(in), asynchronous :: sendbuf
            type(*), dimension(..) :: recvbuf
            integer, intent(in) :: sendcount, recvcount
            integer, intent(in) :: sendtype, recvtype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Isend(buf, count, datatype, dest, tag, comm, request, ierror)
            type(*), dimension(..), intent(in), asynchronous :: buf
            integer, intent(in) :: count, dest, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, intent(out) :: request
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Recv(buf, count, datatype, source, tag, comm, status, ierror)
            type(*), dimension(..) :: buf
            integer, intent(in) :: count, source, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer :: status
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Irecv(buf, count, datatype, source, tag, comm, request, ierror)
            type(*), dimension(..), asynchronous :: buf
            integer, intent(in) :: count, source, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, intent(out) :: request
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Waitall(count, array_of_requests, array_of_statuses, ierror)
            integer, intent(in) :: count
            integer, intent(inout) :: array_of_requests(count)
            integer :: array_of_statuses(*)
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Ssend(buf, count, datatype, dest, tag, comm, ierror)
            type(*), dimension(..), intent(in) :: buf
            integer, intent(in) :: count, dest, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
        end subroutine
    end interface
end module
