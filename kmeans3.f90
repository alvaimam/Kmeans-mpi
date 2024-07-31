program read_from_file
  use functions_module
  use mpi
  character(len=2) :: c2
  integer :: i, j, k, l
  real, dimension(:,:), allocatable :: r, centroid, new_centro, recv_buf
  real, dimension (:,:), allocatable :: sendcounts, displs
  real, dimension(:), allocatable :: cost
  integer,dimension(:),allocatable  :: indices,distancereg,cluster
  integer :: ios_read = 0
  integer :: n = 0 
  integer :: omega, n_buf, ierr, my_id, snd_id, rcv_id, n_proc

  print *, 'which data index?'
  read*, idx

  write(c2, '(i2.2)') idx   
  !open(unit=99, file='kmeans_data_distrib_'//c2//'.dat')
  ! open(unit=99, file='kmeans_data_distrib_'//c2//'_small.dat', iostat=ios_read)
!Here I use small dataset first for trial. 
   open(unit=99, file='datatest1.dat', iostat=ios_read)
  if (ios_read /= 0) then
    print *, "kmeans_data_distrib_"//c2//"_small.dat could not be opened"
    ! print
  end if

  !find the maximum lines
  do
    read(99, *, iostat=ios_read) i, x, y

    if (ios_read > 0) then
      print *, "something is wrong"
      stop
    else if (ios_read < 0) then
      print *, "end of file reached"
      exit
    else
      n = n+1
    end if
  end do

  rewind(99)

  !do i=1,n
  open(unit=98, file='rawdata.dat')
  allocate(r(2, n))
  do i = 1,n
    read(99, *, iostat=ios_read) j, x, y
    
    r(1, j) = x
    r(2, j) = y
    write(98, *) x, y
  end do

  close(99)   ! close kmeans
  close(98)   ! close rawdatai
  call mpi_init(ierr)
  
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  recv_buf = 0 
  allocate (recv_buf(2,nproc))
  if (my_id == 0) then 
        print*, 'this is root process'
        call mpi_scatterv(r,1,mpi_real,recv_buf,1,mpi_real,0,mpi_comm_world,ierr)
  end if
  print*, 'put k'
  read*, k
  call centroid_inits(r, n, k, centroid)
  do l = 1, 10
     ! change here the limit and control convergence by an eps-value
     ! eps > |c_Old - c_new|
     
     call min_distance(r, n, k, centroid, distance,indices,distancereg)
     call new_centroid(r,n,k,centroid,indices,new_centro,omega)
     call costfunction(r,n,k,distancereg,indices,new_centro,cluster,cost)
  enddo
  
  open(unit=99,file="kmeans3_test.dat")
  do i = 1, n
     write(99,"(2es14.5,i4)") r(:,i),indices(i)
  enddo
  close(99)
  
  
Contains
  subroutine centroid_inits(r,n,k,centroid)
      real,dimension (:,:),intent(in),allocatable :: r
      real,dimension (:,:),intent(out),allocatable:: centroid
      real,dimension(k),allocatable::xc(:) ,yc(:)
      integer,intent(in) :: n,k
      integer :: i
      real :: maks_x,maks_y,min_x,min_y
      allocate(centroid(2, k))
      allocate(xc(k))
      allocate(yc(k))

      maks_x = maxval(r(1,:))
      maks_y = maxval(r(2,:))
      min_x = minval(r(1,:))
      min_y = minval(r(2,:))
     ! print *, min_x, maks_x, min_y, maks_y

      do i = 1,k
        xc (i) = min_x + (maks_x - min_x) * fib_rnd()
        yc (i) = min_y + (maks_y - min_y) * fib_rnd()
        centroid (1,i) = xc(i)
        centroid (2,i) = yc(i)
      end do

      do i = 1,k 
        print *, centroid(:,i)
      end do
    end subroutine centroid_inits

  subroutine min_distance(r,n,k,centroid,distance,indices,distancereg)
        integer, intent(out):: n,k
        real,dimension(:,:),intent(in),allocatable::centroid
        real,dimension(:,:),intent(in),allocatable::r
        integer,dimension(:),intent(out),allocatable::indices,distancereg
        real ::d_min
        integer::y,i_min,j,i
        integer,parameter :: data_dim=2
        allocate (indices(n))
        allocate (distancereg(k))
        !cost=0.d0
        do j=1,n

           i_min = -1
           d_min=1.d6
           do i=1,k
              distance=0.d0
              distancereg(i)=0.d0
              do y=1,data_dim
                 distance = distance+abs(r(y,j)-centroid(y,i)) 
                 distancereg(i)=distancereg(i)+abs(r(y,j)-centroid(y,i))
              end do
              if (distance<d_min) then
                 d_min=distance
                 i_min=i
              end if
           end do
           if( i_min < 0 ) print*," found error by assigning k-index to particle ",j
           indices(j)=i_min
        end do

  end subroutine
 subroutine new_centroid(r,n,k,centroid,indices,new_centro,omega)
      integer, intent(in):: n
      real,dimension(:,:),intent(inout),allocatable ::centroid
      real,dimension(:,:),intent(in),allocatable ::r
      integer,dimension(:),intent(in),allocatable::indices
      real,dimension(:,:),intent(out),allocatable:: new_centro
      integer,intent(inout)::k
      integer :: t,y,j,k_ind
      integer,intent(out) :: omega
      real,dimension(:),allocatable :: summ
      allocate(summ(2))
      allocate (new_centro(2,k))
       t=2
      do k_ind=1,k
        
        omega = count(indices==k_ind)
        summ(1)=0
        summ(2)=0 
        do j=1,n 
            if (indices(j)==k_ind) then
               summ(1) =+ r(1,j) 
               summ(2) =+ r(2,j)
            end if 
        end do
      new_centro(1,k_ind) = summ(1)/omega
      new_centro(2,k_ind) = summ(2)/omega 
   end do

   !GS
   centroid = new_centro
   
      do k_ind=1,k
      print*, 'new centro',new_centro(:,k_ind)
      end do

 end subroutine           
 subroutine costfunction(r,n,k,distancereg,indices,new_centro,cluster,cost)
    integer, dimension (:), allocatable, intent(out)  :: distancereg, indices
    integer, dimension (:), allocatable, intent(out) :: cluster
    real, dimension (:,:), allocatable, intent(in) :: r
    real, dimension (:,:), intent(in), allocatable :: new_centro 
    real, dimension(:), intent(out), allocatable :: cost
    integer :: i,k
    allocate(cluster(k))
    allocate(cost(k))
    allocate(distancereg(k))       

    call min_distance(r,n,k,centroid,distance,indices,distancereg)
    cluster = 0
    do i=1,k
       cost(i)=0
       cluster(i)=count(indices==i)
       cost(i)=(1.0/cluster(i))*distancereg(i)
       print*,cost(i)
    end do

    print*," total sum of cluster members ",sum(cluster)," vs N ",n
    
 end subroutine

!subroutine mindistance_centro(r,n,k,centroid,distance,indices,distancereg)

call mpi_finalize(ierr)
end program read_from_file 













!#ifdef ALL
! subroutine get_cluster(r,n,k,omega,cluster)
!    integer, intent (in) :: n,k
!    real, dimension (:,:),intent(out),allocatable :: r
!    integer :: y,i,j,z
!    integer,intent(out):: omega
    
!    call(min_distance(r,n,k,new_centro,distance,indices))
!    do i=1,n
!      do j=1,k
!        if (r(1,n)-new_centro(1,k)
       
    
    
! end subroutine  
!#endif

 
! !   subroutine new_label (r,n,k,new_centro,distance
! !     integer, intent(in)::n,k
! !     real, dimension (:,:),intent(in),allocatable::centroid,r
! !     real, dimension (:,:),intent (out),allocatable:: distance
! !
! !     integer ::y,i,j,z
! !     allocate(distance(2,i)) 
! !     i = 0
! !     do j=1,n
! !        do y=1,k
! !           i = i+1
! !          distance(1,i) = abs(r(1,j)-centroid(1,y))
! !          distance(2,i) = abs(r(2,j)-centroid(2,y))
! !   
! !         enddo     
! !     end do    
  !
  !
! y=y+1      ! distance(2,i)=abs( r(2,i) - centroid (2,j) 
        !end do
      ! do z=1,i
      !  print*, 'this is distance', distance (1,z),distance (2,z)
      !end do
      !mindist_x=minval(distance(1,:))
      !mindist_y=minval(distance(2,:))  
!    end subroutine min_distance
!
!    subroutine get_cluster (r,n,k,centroid,distance,cluster) 
!      integer,intent(in)::n,k
!      real,dimension(:,:),intent(in),allocatable::centroid,r,distance
!
!      real,dimension(:),intent(out),allocatable :: cluster
!      
!      integer :: i,j
!      allocate (cluster(2,j))
!      !call min_distance(r,n,k,centroid,distance)
!      i=1
!      j=2
!      if (distance(:,i) < distance(:,j)) then
!         cluster(:,i) = r(:,i)
!         else 
!         cluster(:,j) = r(:,j)
!     end if
!      
!        
!
!    end subroutine get_cluster

    !subroutine find_cost
    !end subroutine find_cost

    !subroutine new_centroid
    !end subroutine new_centroid

    !subroutine new_cluster
    !end subroutine new_cluster

!end program read_from_file
