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
    integer, parameter :: MPI_STATUSES_IGNORE = 0
    integer, parameter :: MPI_STATUS_IGNORE = 0

    interface
        function MPI_Wtime() result(time)
            real(kind=8) :: time
        end function
    end interface
end module
