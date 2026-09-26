program kmeans_scan_mpi
  !-------------------------------------------------------------------------
  ! MPI version of the k-means scan.
  !
  ! Parallelization strategy:
  !   For each (k, Nu) block, the Nu independent k-means restarts are split
  !   across MPI ranks (block distribution). Each rank runs its share
  !   locally and computes its own cost values; MPI_Gatherv collects them
  !   all onto rank 0, which writes the output files -- so the output
  !   format is IDENTICAL to the serial version.
  !
  !   This is "embarrassingly parallel": no communication happens during
  !   the actual k-means iterations, only the final gather per block.
  !
  ! Output (written by rank 0 only):
  !   mincost_vs_k.dat   -> columns: k   Nu   min_cost
  !   hist_nuXXXX.dat    -> columns: k   run   cost
  !
  ! Build:  mpif90 -O2 kmeans3_mpi.f90 functions.f90 -o kmeans3_mpi
  ! Run:    mpirun -np 4 ./kmeans3_mpi
  !-------------------------------------------------------------------------
  use functions_module   ! for fib_rnd(), rnd_seed
  use mpi
  implicit none

  ! ---- data ----
  integer :: n, idx, i, j, ios_read
  real, allocatable :: r(:,:)
  character(len=2) :: c2

  ! ---- k-means control ----
  integer, parameter :: kmax     = 30
  integer, parameter :: max_iter = 100
  real,    parameter :: eps      = 1.0e-5

  ! ---- experiment control ----
  integer, parameter :: n_nu_vals = 4
  integer :: nu_list(n_nu_vals) = (/ 1, 10, 100, 1000 /)
  integer :: nu, i_nu, k, irun
  real :: mincost, run_cost

  ! ---- MPI ----
  integer :: ierr, my_id, nproc
  integer, allocatable :: counts(:), displs(:)
  integer :: local_n
  real, allocatable :: local_costs(:), all_costs(:)

  integer :: unit_summary, unit_hist
  character(len=64) :: fname

  !-------------------------------------------------------------------------
  ! 0) MPI startup
  !-------------------------------------------------------------------------
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)

  ! give every rank an independent random stream (see note above)
  rnd_seed = rnd_seed + real(my_id, kind(rnd_seed)) * 104729.0d0

  !-------------------------------------------------------------------------
  ! 1) rank 0 reads the data, then broadcasts it to everyone
  !-------------------------------------------------------------------------
  if (my_id == 0) then
     print *, 'which data index?'
     read *, idx
     write(c2,'(i2.2)') idx

     open(unit=99, file='datatest1.dat', iostat=ios_read)
     if (ios_read /= 0) then
        print *, 'datatest1.dat could not be opened'
        call mpi_abort(mpi_comm_world, 1, ierr)
     end if

     n = 0
     do
        read(99,*,iostat=ios_read) i
        if (ios_read /= 0) exit
        n = n + 1
     end do
     rewind(99)

     allocate(r(2,n))
     do i = 1, n
        read(99,*) j, r(1,i), r(2,i)
     end do
     close(99)

     print *, 'read ', n, ' data points, running on ', nproc, ' ranks'
  end if

  call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
  if (my_id /= 0) allocate(r(2,n))
  call mpi_bcast(r, 2*n, mpi_real, 0, mpi_comm_world, ierr)

  !-------------------------------------------------------------------------
  ! 2) scan k = 1..30 for each Nu; parallelize the Nu restarts
  !-------------------------------------------------------------------------
  if (my_id == 0) then
     open(newunit=unit_summary, file='mincost_vs_k.dat', status='replace')
     write(unit_summary,'(a)') '# k   Nu   min_cost'
  end if

  do i_nu = 1, n_nu_vals
     nu = nu_list(i_nu)

     if (my_id == 0) then
        write(fname,'(a,i4.4,a)') 'hist_nu', nu, '.dat'
        open(newunit=unit_hist, file=trim(fname), status='replace')
        write(unit_hist,'(a)') '# k   run   cost'
     end if

     ! block-distribute the Nu restarts across ranks
     allocate(counts(0:nproc-1), displs(0:nproc-1))
     call compute_block_distribution(nu, nproc, counts, displs)
     local_n = counts(my_id)

     do k = 1, kmax

        if (allocated(local_costs)) deallocate(local_costs)
        allocate(local_costs(max(local_n,1)))

        do irun = 1, local_n
           call run_kmeans_once(r, n, k, max_iter, eps, run_cost)
           local_costs(irun) = run_cost
        end do

        if (my_id == 0) then
           if (allocated(all_costs)) deallocate(all_costs)
           allocate(all_costs(nu))
        end if

        call mpi_gatherv(local_costs, local_n, mpi_real,                 &
                          all_costs, counts, displs, mpi_real,           &
                          0, mpi_comm_world, ierr)

        if (my_id == 0) then
           mincost = minval(all_costs)
           write(unit_summary,'(i6,i8,es16.6)') k, nu, mincost
           do irun = 1, nu
              write(unit_hist,'(2i6,es16.6)') k, irun, all_costs(irun)
           end do
           print '(a,i3,a,i5,a,es12.4)', ' k=', k, '  Nu=', nu, '  min cost=', mincost
        end if

     end do

     deallocate(counts, displs)
     if (my_id == 0) close(unit_hist)

  end do

  if (my_id == 0) then
     close(unit_summary)
     print *, 'Done. See mincost_vs_k.dat and hist_nu*.dat'
  end if

  call mpi_finalize(ierr)

contains

  !-------------------------------------------------------------------------
  subroutine compute_block_distribution(total, nproc, counts, displs)
    integer, intent(in)  :: total, nproc
    integer, intent(out) :: counts(0:nproc-1), displs(0:nproc-1)
    integer :: base, remainder, p

    base      = total / nproc
    remainder = mod(total, nproc)

    do p = 0, nproc-1
       if (p < remainder) then
          counts(p) = base + 1
       else
          counts(p) = base
       end if
    end do

    displs(0) = 0
    do p = 1, nproc-1
       displs(p) = displs(p-1) + counts(p-1)
    end do
  end subroutine compute_block_distribution

  !-------------------------------------------------------------------------
  subroutine run_kmeans_once(r, n, k, max_iter, eps, cost_out)
    real,    intent(in)  :: r(:,:)
    integer, intent(in)  :: n, k, max_iter
    real,    intent(in)  :: eps
    real,    intent(out) :: cost_out

    real, allocatable :: centroid(:,:), new_centro(:,:)
    integer, allocatable :: indices(:), csize(:)
    real :: diff
    integer :: it

    call init_centroids(r, n, k, centroid)

    allocate(indices(n))
    allocate(new_centro(2,k))
    allocate(csize(k))

    do it = 1, max_iter
       call assign_clusters(r, n, k, centroid, indices)
       call update_centroids(r, n, k, indices, new_centro, csize)

       diff = maxval(abs(new_centro - centroid))
       centroid = new_centro

       if (diff < eps) exit
    end do

    call compute_cost(r, n, k, centroid, indices, cost_out)

    deallocate(centroid, new_centro, indices, csize)
  end subroutine run_kmeans_once

  !-------------------------------------------------------------------------
  subroutine init_centroids(r, n, k, centroid)
    real,    intent(in)  :: r(:,:)
    integer, intent(in)  :: n, k
    real, allocatable, intent(out) :: centroid(:,:)
    real :: xmin, xmax, ymin, ymax
    integer :: i

    allocate(centroid(2,k))
    xmin = minval(r(1,:)); xmax = maxval(r(1,:))
    ymin = minval(r(2,:)); ymax = maxval(r(2,:))

    do i = 1, k
       centroid(1,i) = xmin + (xmax - xmin) * fib_rnd()
       centroid(2,i) = ymin + (ymax - ymin) * fib_rnd()
    end do
  end subroutine init_centroids

  !-------------------------------------------------------------------------
  subroutine assign_clusters(r, n, k, centroid, indices)
    real,    intent(in)  :: r(:,:), centroid(:,:)
    integer, intent(in)  :: n, k
    integer, intent(out) :: indices(:)
    real :: d, dmin
    integer :: i, j

    do j = 1, n
       dmin = huge(1.0)
       indices(j) = 1
       do i = 1, k
          d = sum((r(:,j) - centroid(:,i))**2)
          if (d < dmin) then
             dmin = d
             indices(j) = i
          end if
       end do
    end do
  end subroutine assign_clusters

  !-------------------------------------------------------------------------
  subroutine update_centroids(r, n, k, indices, new_centro, csize)
    real,    intent(in)  :: r(:,:)
    integer, intent(in)  :: n, k, indices(:)
    real,    intent(out) :: new_centro(:,:)
    integer, intent(out) :: csize(:)
    integer :: i, j

    new_centro = 0.0
    csize = 0

    do j = 1, n
       i = indices(j)
       new_centro(:,i) = new_centro(:,i) + r(:,j)
       csize(i) = csize(i) + 1
    end do

    do i = 1, k
       if (csize(i) > 0) then
          new_centro(:,i) = new_centro(:,i) / real(csize(i))
       else
          ! empty cluster: re-seed on a random data point
          new_centro(:,i) = r(:, 1 + int(fib_rnd() * real(n-1)))
       end if
    end do
  end subroutine update_centroids

  !-------------------------------------------------------------------------
  subroutine compute_cost(r, n, k, centroid, indices, cost)
    real,    intent(in)  :: r(:,:), centroid(:,:)
    integer, intent(in)  :: n, k, indices(:)
    real,    intent(out) :: cost
    integer :: j

    cost = 0.0
    do j = 1, n
       cost = cost + sum((r(:,j) - centroid(:,indices(j)))**2)
    end do
    cost = cost / real(n)
  end subroutine compute_cost

end program kmeans_scan_mpi
