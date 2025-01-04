module mpi_c_bindings
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine c_mpi_bcast_int(buffer, count, datatype, root, comm, ierror) bind(C, name="MPI_Bcast")
            import c_int
            integer(c_int) :: buffer
            integer(c_int), intent(in) :: count, root
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_bcast_real(buffer, count, datatype, root, comm, ierror) bind(C, name="MPI_Bcast")
            import :: c_int, c_double
            real(c_double), dimension(:, :) :: buffer
            integer(c_int), intent(in) :: count, root
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allgather_int(sendbuf, sendcount, sendtype, recvbuf, &
                                       recvcount, recvtype, comm, ierror) bind(C, name="MPI_Allgather")
            import :: c_int, c_double
            integer(c_int), dimension(:), intent(in) :: sendbuf
            integer(c_int), dimension(*) :: recvbuf
            integer(c_int), intent(in) :: sendcount, recvcount
            integer(c_int), intent(in) :: sendtype, recvtype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allgather_real(sendbuf, sendcount, sendtype, recvbuf, &
                                        recvcount, recvtype, comm, ierror) bind(C, name="MPI_Allgather")
            import :: c_int, c_double
            real(c_double), dimension(:), intent(in) :: sendbuf
            real(c_double), dimension(*) :: recvbuf
            integer(c_int), intent(in) :: sendcount, recvcount
            integer(c_int), intent(in) :: sendtype, recvtype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_isend_2d(buf, count, datatype, dest, tag, comm, request, ierror) bind(C, name="MPI_Isend")
            import :: c_int, c_double
            real(c_double), dimension(:, :), intent(in) :: buf
            integer(c_int), intent(in) :: count, dest, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: request
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_isend_3d(buf, count, datatype, dest, tag, comm, request, ierror) bind(C, name="MPI_Isend")
            import :: c_int, c_double
            real(c_double), dimension(:, :, :), intent(in) :: buf
            integer(c_int), intent(in) :: count, dest, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: request
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror) bind(C, name="MPI_Allreduce")
            import :: c_int, c_double
            real(c_double), intent(in) :: sendbuf
            real(c_double), intent(in) :: recvbuf
            integer(c_int) :: count, datatype, op, comm, ierror
        end subroutine

        subroutine c_mpi_allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror) bind(C, name="MPI_Allreduce")
            import :: c_int, c_double
            real(c_double), intent(in) :: sendbuf
            real(c_double), dimension(:), intent(in) :: recvbuf
            integer(c_int) :: count, datatype, op, comm, ierror
        end subroutine

        function c_mpi_wtime() result(time) bind(C, name="MPI_Wtime")
            import :: c_double
            real(c_double) :: time
        end function

        subroutine c_mpi_barrier(comm, ierror) bind(C, name="MPI_Barrier")
            import :: c_int
            integer(c_int) :: comm, ierror
        end subroutine
        
        subroutine c_mpi_finalize(ierror) bind(C, name="MPI_Finalize")
            import :: c_int
            integer(c_int), optional, intent(in) :: ierror
        end subroutine

        subroutine c_mpi_init_thread(required, provided, ierror) bind(C, name="MPI_Init_thread")
            import :: c_int
            integer(c_int), intent(in) :: required
            integer(c_int), intent(out) :: provided
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_comm_size(comm, size, ierror) bind(C, name="MPI_Comm_size")
            import :: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: size
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_comm_rank(comm, rank, ierror) bind(C, name="MPI_Comm_rank_proc")
            import :: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: rank
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_comm_split_type(comm, split_type, key, info, newcomm, ierror) bind(C, name="MPI_Comm_split_type")
            import :: c_int
            integer(c_int) :: comm
            integer(c_int), intent(in) :: split_type, key
            integer(c_int), intent(in) :: info
            integer(c_int), intent(out) :: newcomm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierror) bind(C, name="MPI_Cart_create")
            import :: c_int
            integer(c_int), intent(in) :: comm_old
            integer(c_int), intent(in) :: ndims
            integer(c_int), intent(in) :: dims(ndims)
            integer(c_int), intent(in) :: periods(ndims), reorder  ! logical -> integer(c_int), I think that's correct
            integer(c_int), intent(out) :: comm_cart
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_cart_coords(comm, rank, maxdims, coords, ierror) bind(C, name="MPI_Cart_coords")
            import :: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(in) :: rank, maxdims
            integer(c_int), intent(out) :: coords(maxdims)
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_cart_shift(comm, direction, disp, rank_source, rank_dest, ierror) bind(C, name="MPI_Cart_shift")
            import :: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(in) :: direction, disp
            integer(c_int), intent(out) :: rank_source, rank_dest
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_cart_sub(comm, remain_dims, newcomm, ierror) bind(C, name="MPI_Cart_sub")
            import :: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(in) :: remain_dims(*)  ! logical -> integer(c_int)
            integer(c_int), intent(out) :: newcomm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_recv(buf, count, datatype, source, tag, comm, status, ierror) bind(C, name="MPI_Recv")
            import :: c_int, c_double
            real(c_double), dimension(..) :: buf
            integer(c_int), intent(in) :: count, source, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int) :: status
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror) bind(C, name="MPI_Irecv")
            import :: c_int, c_double
            real(c_double), dimension(..) :: buf
            integer(c_int), intent(in) :: count, source, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: request
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_waitall(count, array_of_requests, array_of_statuses, ierror) bind(C, name="MPI_Waitall")
            import :: c_int
            integer(c_int), intent(in) :: count
            integer(c_int), intent(inout) :: array_of_requests(count)
            integer(c_int) :: array_of_statuses(*)
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_ssend(buf, count, datatype, dest, tag, comm, ierror) bind(C, name="MPI_Ssend")
            import :: c_int, c_double
            real(c_double), dimension(..), intent(in) :: buf
            integer(c_int), intent(in) :: count, dest, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine
    end interface
end module
