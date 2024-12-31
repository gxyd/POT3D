module mpi
    use mpi_c_bindings
    implicit none
    integer, parameter :: MPI_COMM_WORLD = 0
    integer, parameter :: MPI_COMM_TYPE_SHARED = 1
    integer, parameter :: MPI_INFO_NULL = 0
    integer, parameter :: MPI_REAL4 = 2
    integer, parameter :: MPI_REAL8 = 3
    integer, parameter :: MPI_THREAD_FUNNELED = 1
    integer, parameter :: MPI_INTEGER = 1
    real(8), parameter :: MPI_IN_PLACE = 1
    integer, parameter :: MPI_SUM = 1
    integer, parameter :: MPI_PROC_NULL = -1
    integer, allocatable :: MPI_STATUSES_IGNORE(:)
    integer, parameter :: MPI_STATUS_IGNORE = 0

    interface MPI_Bcast
        module procedure MPI_Bcast_int
        module procedure MPI_Bcast_real
    end interface

    interface MPI_Allgather
        module procedure MPI_Allgather_int
        module procedure MPI_Allgather_real
    end interface

    interface MPI_Isend
        module procedure MPI_Isend_2d
        module procedure MPI_Isend_3d

    end interface

    interface MPI_Allreduce
        module procedure MPI_Allreduce_scalar
        module procedure MPI_Allreduce_1d
    end interface

    interface MPI_Wtime
        module procedure MPI_Wtime_proc
    end interface

    interface MPI_Barrier
        module procedure MPI_Barrier_proc
    end interface

    interface MPI_Finalize
        module procedure MPI_Finalize_proc
    end interface
    
    interface MPI_Init_thread
        module procedure MPI_Init_thread_proc
    end interface

    interface MPI_Comm_size
        module procedure MPI_Comm_size_proc
    end interface

    interface MPI_Comm_rank
        module procedure MPI_Comm_rank_proc
    end interface

    interface MPI_Comm_split_type
        module procedure MPI_Comm_split_type_proc
    end interface

    interface MPI_Cart_create
        module procedure MPI_Cart_create_proc
    end interface

    interface MPI_Cart_coords
        module procedure MPI_Cart_coords_proc
    end interface

    interface MPI_Cart_shift
        module procedure MPI_Cart_shift_proc
    end interface

    interface MPI_Cart_sub
        module procedure MPI_Cart_sub_proc
    end interface

    interface MPI_Recv
        module procedure MPI_Recv_proc
    end interface

    interface MPI_Irecv
        module procedure MPI_Irecv_proc
    end interface

    interface MPI_Waitall
        module procedure MPI_Waitall_proc
    end interface

    interface MPI_Ssend
        module procedure MPI_Ssend_proc
    end interface

    contains
        subroutine MPI_Bcast_int(buffer, count, datatype, root, comm, ierror)
            integer :: buffer
            integer, intent(in) :: count, root
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
            call c_mpi_bcast_int(buffer, count, datatype, root, comm, ierror)
        end subroutine

        subroutine MPI_Bcast_real(buffer, count, datatype, root, comm, ierror)
            real(8), dimension(:, :) :: buffer
            integer, intent(in) :: count, root
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
            call c_mpi_bcast_real(buffer, count, datatype, root, comm, ierror)
        end subroutine

        subroutine MPI_Allgather_int(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
            ! type(*), dimension(..), intent(in), asynchronous :: sendbuf
            integer, dimension(:), intent(in) :: sendbuf
            ! type(*), dimension(..) :: recvbuf
            integer, dimension(:, :) :: recvbuf
            integer, intent(in) :: sendcount, recvcount
            integer, intent(in) :: sendtype, recvtype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
            call c_mpi_allgather_int(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
        end subroutine

        subroutine MPI_Allgather_real(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
            ! type(*), dimension(..), intent(in), asynchronous :: sendbuf
            real(8), dimension(:), intent(in) :: sendbuf
            ! type(*), dimension(..) :: recvbuf
            real(8), dimension(:, :) :: recvbuf
            integer, intent(in) :: sendcount, recvcount
            integer, intent(in) :: sendtype, recvtype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
            call c_mpi_allgather_real(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
        end subroutine

        subroutine MPI_Isend_2d(buf, count, datatype, dest, tag, comm, request, ierror)
            real(8), dimension(:, :), intent(in) :: buf
            integer, intent(in) :: count, dest, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, intent(out) :: request
            integer, optional, intent(out) :: ierror
            call c_mpi_isend_2d(buf, count, datatype, dest, tag, comm, request, ierror)
        end subroutine

        subroutine MPI_Isend_3d(buf, count, datatype, dest, tag, comm, request, ierror)
            real(8), dimension(:, :, :), intent(in) :: buf
            integer, intent(in) :: count, dest, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, intent(out) :: request
            integer, optional, intent(out) :: ierror
            call c_mpi_isend_3d(buf, count, datatype, dest, tag, comm, request, ierror)
        end subroutine

        subroutine MPI_Allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror)
            real(8), intent(in) :: sendbuf
            real(8), intent(in) :: recvbuf
            integer :: count, datatype, op, comm, ierror
            call c_mpi_allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        end subroutine

        subroutine MPI_Allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror)
            real(8), intent(in) :: sendbuf
            real(8), dimension(:), intent(in) :: recvbuf
            integer :: count, datatype, op, comm, ierror
            call c_mpi_allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        end subroutine

        function MPI_Wtime_proc() result(time)
            real(8) :: time
            time = c_mpi_wtime()
        end function

        subroutine MPI_Barrier_proc(comm, ierror)
            integer :: comm, ierror
            call c_mpi_barrier(comm, ierror)
        end subroutine

        subroutine MPI_Finalize_proc(ierror)
            integer, optional, intent(in) :: ierror
            call c_mpi_finalize(ierror)
        end subroutine

        subroutine MPI_Init_thread_proc(required, provided, ierror)
            integer, intent(in) :: required
            integer, intent(out) :: provided
            integer, optional, intent(out) :: ierror
            call c_mpi_init_thread(required, provided, ierror)
        end subroutine

        subroutine MPI_Comm_size_proc(comm, size, ierror)
            integer, intent(in) :: comm
            integer, intent(out) :: size
            integer, optional, intent(out) :: ierror
            call c_mpi_comm_size(comm, size, ierror)
        end subroutine

        subroutine MPI_Comm_rank_proc(comm, rank, ierror)
            integer, intent(in) :: comm
            integer, intent(out) :: rank
            integer, optional, intent(out) :: ierror
            call c_mpi_comm_rank(comm, rank, ierror)
        end subroutine

        subroutine MPI_Comm_split_type_proc(comm, split_type, key, info, newcomm, ierror)
            integer :: comm
            integer, intent(in) :: split_type, key
            integer, intent(in) :: info
            integer, intent(out) :: newcomm
            integer, optional, intent(out) :: ierror
            call c_mpi_comm_split_type(comm, split_type, key, info, newcomm, ierror)
        end subroutine

        subroutine MPI_Cart_create_proc(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
            integer, intent(in) :: comm_old
            integer, intent(in) :: ndims
            integer, intent(in) :: dims(ndims)
            logical, intent(in) :: periods(ndims), reorder
            integer, intent(out) :: comm_cart
            integer, optional, intent(out) :: ierror
            ! call c_mpi_cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror)
        end subroutine

        ! from here
        subroutine MPI_Cart_coords_proc(comm, rank, maxdims, coords, ierror)
            integer, intent(in) :: comm
            integer, intent(in) :: rank, maxdims
            integer, intent(out) :: coords(maxdims)
            integer, optional, intent(out) :: ierror
            call c_mpi_cart_coords(comm, rank, maxdims, coords, ierror)
        end subroutine

        subroutine MPI_Cart_shift_proc(comm, direction, disp, rank_source, rank_dest, ierror)
            integer, intent(in) :: comm
            integer, intent(in) :: direction, disp
            integer, intent(out) :: rank_source, rank_dest
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Cart_sub_proc(comm, remain_dims, newcomm, ierror)
            integer, intent(in) :: comm
            logical, intent(in) :: remain_dims(*)
            integer, intent(out) :: newcomm
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Recv_proc(buf, count, datatype, source, tag, comm, status, ierror)
            real(8), dimension(..) :: buf
            integer, intent(in) :: count, source, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer :: status
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Irecv_proc(buf, count, datatype, source, tag, comm, request, ierror)
            real(8), dimension(..) :: buf
            integer, intent(in) :: count, source, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, intent(out) :: request
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Waitall_proc(count, array_of_requests, array_of_statuses, ierror)
            integer, intent(in) :: count
            integer, intent(inout) :: array_of_requests(count)
            integer :: array_of_statuses(*)
            integer, optional, intent(out) :: ierror
        end subroutine

        subroutine MPI_Ssend_proc(buf, count, datatype, dest, tag, comm, ierror)
            real(8), dimension(..), intent(in) :: buf
            integer, intent(in) :: count, dest, tag
            integer, intent(in) :: datatype
            integer, intent(in) :: comm
            integer, optional, intent(out) :: ierror
        end subroutine
end module

program main
    implicit none
end program main
