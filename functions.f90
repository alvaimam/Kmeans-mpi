      module functions_module

        use data_kind

        implicit none

        private
        
        public  ::  factorial
        public  ::  poisson_distrib, binomial_distrib, exponential_distrib
        public  ::  rnd, gasdev,gauss,cauchy_fct,cauchy_distrib
        public  ::  bernoulli_distrib, uniform_distrib, gauss_distrib
        public  ::  norm_hist, cum_norm_inv, cum_norm, gammln
        public  ::  levy, rnd_lcg, triang_distrib, triangle_fct
        public  ::  beta_distrib, empirical_distrib
        public  ::  log_norm_fct, sgasdev
        public  ::  triang_distrib_sig_mu,triang_distrib_xxxy
        public  ::  log2,p3d
        public  ::  fib_rnd, fib_srnd, srnd
        public  ::  bsp0, bsp1, bsp2, bsp3, bsp4
        public  ::  pot_bsp1,pot_bsp2,pot_bsp3,pot_bsp4
        public  ::  u_bsp, u_bsp0, u_bsp1, u_bsp2, u_bsp3, u_bsp4
        public  ::  e_bsp, e_bsp0, e_bsp1, e_bsp2, e_bsp3, e_bsp4
        public  ::  f_bsp, f_bsp0, f_bsp1, f_bsp2, f_bsp3, f_bsp4
        
        real (kind = real_8), parameter, private ::  pi = acos( -1.d0 )
        real (kind = real_8), parameter, private ::  d2p31m = 2147483647.0 ,        &
             d2p31  = 2147483648.0
        real( kind = real_8 ), save, public      ::  rnd_seed = 1.0

        character (len=*), parameter               ::  module_name = "functions_module"
        
        interface norm_hist
           module procedure norm_hist_1, norm_hist_2, norm_hist_3
        end interface

        interface triang_distrib
           module procedure triang_distrib_sig_mu,triang_distrib_xxxy
        end interface

        interface triangle_fct
           module procedure triangle_fct_x_sig_mu, triangle_fct_xxxy
        end interface


        interface rnd_lcg
           module procedure  rnd_lcg_r, rnd_lcg_i
        end interface
        interface factorial
           module procedure  factorial_i, factorial_r8!, factorial_r6
        end interface

      contains

        
        function f_bsp(r,rc,order)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of given order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          integer, intent(in)                   ::  order
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )
          character (len=*), parameter          ::  method_name = "f_bsp"

          select case( order )
!!$          case( -1 )
!!$             f = f_cos0(r,rc)

          case( 0 )
             f = f_bsp0(r,rc)

          case( 1 )
             f = f_bsp1(r,rc)

          case( 2 )
             f = f_bsp2(r,rc)

          case( 3 )
             f = f_bsp3(r,rc)

          case( 4 )
             f = f_bsp4(r,rc)

          case default
             write(6,"(//,a,i3,a,//,a,//)") module_name//"::"//method_name//": Order ",  &
                  order," spline not yet implemented as source term",    &
                  " >>>>>  program stops  <<<<<"

          end select

        end function f_bsp


        function f_bsp0(r,rc)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of 2nd order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r > 1.d0*rc )then
             f = 0.d0

          else if( r <= 1.d0*rc )then

             f = 3.d0 / (4.d0*pi*rc**3) 

          endif

        end function f_bsp0


        function f_bsp1(r,rc)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of 2nd order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )

          if( r >= rc )then
             f = 0.d0

          else if( r < rc )then
             f = 3.d0 / (pi*rc**3) * ( 1.d0 - r/rc)

          endif

        end function f_bsp1


        function f_bsp2(r,rc)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of 2nd order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )


          if( r >= rc )then
             f = 0.d0

          else if( r < rc/3.d0 )then
             f = 27.d0 / (4.d0*pi*rc**3) * ( 3.d0 - 9.d0*r**2/rc**2) / 4.d0

          else if( r >= rc/3.d0 .and. r < rc )then
             f = 27.d0 / (4.d0*pi*rc**3) * ( 9.d0 - 18.d0*r/rc + 9.d0*r**2/rc**2 ) / 8.d0

          endif

        end function f_bsp2



        function f_bsp3(r,rc)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of 2nd order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )


          if( r >= rc )then
             f = 0.d0

          else if( r < rc/2.d0 )then
             f = 12.d0 / (pi*rc**3) * (4.d0*r**3/rc**3 - 4.d0*r**2/rc**2 + 2.d0/3.d0)

          else if( r >= rc/2.d0 .and. r < rc )then
             f = 12.d0 / (pi*rc**3) * (-(4.d0/3.d0)*r**3/rc**3 + 4.d0*r**2/rc**2 - 4.d0*r/rc + 4.d0/3.d0)

          endif

        end function f_bsp3


        function f_bsp4(r,rc)       result(f)
          !----------------------------------------------------------------------------------
          !  B-Spline of 2nd order as shape function for charges
          !----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )


          if( r >= rc )then
             f = 0.d0

          else if( r < rc/5.d0 )then
             f = 75.d0 / (4.d0*pi*rc**3) * ((625.d0/64.d0)*r**4/rc**4 - (125.d0/32.d0)*r**2/rc**2 + 115.d0/192.d0)

          else if( r >= rc/5.d0 .and. r < 3.d0*rc/5.d0 )then
             f = 75.d0 / (4.d0*pi*rc**3) * (-(625.d0/96.d0)*r**4/rc**4 + (625.d0/48.d0)*r**3/rc**3   &
                  - (125.d0/16.d0)*r**2/rc**2 + (25.d0/48.d0)*r/rc + 55.d0/96.d0)

          else if( r >= 3.d0*rc/5.d0 .and. r < rc )then
             f = 75.d0 / (4.d0*pi*rc**3) * ((625.d0/384.d0)*r**4/rc**4 - (625.d0/96.d0)*r**3/rc**3   &
                  + (625.d0/64.d0)*r**2/rc**2 - (625.d0/96.d0)*r/rc + 625.d0/384.d0)

          endif

        end function f_bsp4
        
        
        function p3d(sigma)     result(rand)
          !----------------------------------------------------------------------
          
          real*8, intent(in)   ::  sigma
          real*8               ::  rand,xi,eta
          real*8, parameter    ::  sqrtpi2 = sqrt(0.5d0*acos(-1.d0))
          real*8, parameter    ::  pi8inv = sqrt(8.d0/acos(-1.d0))
          
          ! draw exponentially distributed random variate to restrict p
          do 
             eta = fib_rnd()*sigma !*pi8inv
             if( 1.d0-eta < 1e-5 ) cycle
!             xi = -sigma*sqrtpi2*log(1.d0-eta)
             xi = -sigma*sqrtpi2*log(eta)
!             print*,eta,1.d0-pi2inv*eta,log(1.d0-eta),xi
             exit
          enddo
          rand = xi
          
        end function p3d

        
        function bsp0(x,x0,R)     result(f)
          real (kind=real_8), intent(in)  ::  x,x0,R
          real (kind=real_8)              ::  f,z
          z = abs(x - x0)
          if( z < R )then
             f = 1.0/(2.0*R)
          else
             f = 0.0d0
          endif
        end function bsp0

        function bsp1(x,x0,R)     result(f)
          real (kind=real_8), intent(in)  ::  x,x0,R
          real (kind=real_8)              ::  f,z
          z = abs(x - x0)
          if( z < R )then
             f = (1.0 - z/R )/R
          else
             f = 0.0d0
          endif
        end function bsp1

        function bsp2(x,x0,R)     result(f)
          real (kind=real_8), intent(in)  ::  x,x0,R
          real (kind=real_8)              ::  f,z
          z = abs(x - x0)
          if( z < R/3.0d0 )then
             f = 9.0/4.0 - 27.0/4.0*z**2/R**2
          else if( z <= R )then
             f = 27.0d0/8.0d0 - 27.0/4.0*z/R + 27.0d0/8.0d0 * z**2/R**2 
          else
             f = 0.0d0
          endif
          f = f*1.0/(2.0*R)
        end function bsp2

        function bsp3(x,x0,R)     result(f)
          real (kind=real_8), intent(in)  ::  x,x0,R
          real (kind=real_8)              ::  f,z
          z = abs(x - x0)
          if( z < 0.5d0*R )then
             f = 4.0 - 24.0*z**2/R**2 + 24.0*z**3/R**3
          else if( z <= R )then
             f = 8.0d0 - 24.0*z/R + 24.0d0 * z**2/R**2 - 8.0d0 * z**3/R**3
          else
             f = 0.0d0
          endif
          f = f*1.0/(3.0*R)
        end function bsp3

        function bsp4(x,x0,R)     result(f)
          real (kind=real_8), intent(in)  ::  x,x0,R
          real (kind=real_8)              ::  f,z
          z = abs(x - x0)
          if( z < R/5.0d0 )then
             f = 5.0*(23.0*R**4 - 150.0*R**2*z**2 + 375.0*z**4) / (192.0*R**4)
          else if( z <= 3.0*R/5.0 )then
             f = 5.0*(11.0*R**4 + 10.0*R**3*z - 150.0*R**2*z**2 + 250.0*R*z**3 - 125.0*z**4) &
                  / (96.0*R**4)
          else if( z <= R )then
             f = 1.0/384.0 * (5.0-5.0*z/R)**4
          else
             f = 0.0d0
          endif
          f = f*5.0/(2.0*R)
        end function bsp4

        
        function log2(x)   result(f)
          real (kind = real_8), intent(in)  ::  x
          real (kind = real_8)              ::  f
          f = log(x)/log(2.0)
        end function log2

        function exponential_distrib(x,mu)     result(f)
          real (kind = real_8), intent(in)  ::  x,mu
          real (kind = real_8)              ::  f
          if( x < 0.0 ) return
          f = exp(-x/mu)/mu
        end function exponential_distrib

        function binomial_distrib(n,Nmax,a)   result(f)
          real (kind = real_8), intent(in)  ::  a
          integer, intent(in)               ::  n,Nmax
          real (kind = real_8)              ::  f
          integer                           ::  i
          f = (1.0-a)**Nmax
          do i = 1,n
             f = f * a*real(Nmax-i)/((1.0-a)*real(i+1))
          enddo
        end function binomial_distrib
          
        function poisson_distrib( n,x )    result(f)
          real (kind = real_8), intent(in)  ::  x
          integer, intent(in)               ::  n
          real (kind = real_8)              ::  f
          integer                           ::  i
          f = exp(-x)
          do i = 1,n
             f = x*f/real(i)
          enddo
        end function poisson_distrib

        function factorial_r8(x)  result(f)
          real (kind = real_8), intent(in)  ::  x
          real (kind = real_8)              ::  f
          integer                           ::  i
          f=1.0
          do i = 1,nint(x)
             f = f*real(i)
          enddo
        end function factorial_r8

!        function factorial_r6(x)  result(f)
!          real (kind = real_8), intent(in)  ::  x
!          real (kind = real_8)              ::  f
!          integer                           ::  i
!          f=1.0
!          do i = 1,nint(x)
!             f = f*real(i)
!          enddo
!        end function factorial_r6

        function factorial_i(n)    result(nf)
          integer, intent(in)       ::  n
          integer                   ::  i,nf
          nf = 1
          do i = 1,n
             nf = nf * i
          enddo
        end function factorial_i


        function rnd()             result( rand )
          real (kind = real_8)   ::  rand
          rnd_seed = modulo( 16807.0_real_8 * rnd_seed , d2p31m )
          rand     = real( rnd_seed / d2p31 , real_4 )
        end function rnd

        function srnd(seed)             result( rand )
          real (kind = real_8)   ::  rand
          real (kind = real_8), intent(inout)   ::  seed
          seed = modulo( 16807.0_real_8 * seed , d2p31m )
          rand     = real( seed / d2p31 , real_4 )
        end function srnd

        function fib_rnd()         result(rand)
          real (kind = real_8)   ::  rand
          integer, save          ::  i_ctrl = 0
          integer, save          ::  i=17, j=5
          integer                ::  k
          real (kind = real_8), dimension(17), save   ::  x
          if( i_ctrl == 0 )then
             do k = 1,17
                x(k) = rnd()
             enddo
             i_ctrl = 1
          endif
          rand = x(i) - x(j)
          if( rand < 0.0 ) rand = rand+1.0
          x(i) = rand
          i = i-1
          if( i == 0 ) i = 17
          j = j-1
          if( j == 0 ) j = 17
        end function fib_rnd

        function fib_srnd(seed)         result(rand)
          real (kind = real_8)   ::  rand
          real (kind = real_8), intent(inout)   ::  seed
          integer, save          ::  i_ctrl = 0
          integer, save          ::  i=17, j=5
          integer                ::  k
          real (kind = real_8), dimension(17), save   ::  x
          if( i_ctrl == 0 )then
             do k = 1,17
                x(k) = srnd(seed)
             enddo
             i_ctrl = 1
          endif
          rand = x(i) - x(j)
          if( rand < 0.0 ) rand = rand+1.0
          x(i) = rand
          i = i-1
          if( i == 0 ) i = 17
          j = j-1
          if( j == 0 ) j = 17
        end function fib_srnd

!!$        function fib_srnd_par(seed,tid,nth)         result(rand)
!!$          real (kind = real_8)   ::  rand
!!$          real (kind = real_8), intent(inout)   ::  seed
!!$          integer, save          ::  i_ctrl = 0
!!$          integer, save          ::  i=17, j=5
!!$          integer                ::  k
!!$          real (kind = real_8), dimension(17,nth), save   ::  x
!!$          if( i_ctrl == 0 )then
!!$             do k = 1,17
!!$                x(k,tid) = srnd(seed)
!!$             enddo
!!$             i_ctrl = 1
!!$          endif
!!$          rand = x(i,tid) - x(j,tid)
!!$          if( rand < 0.0 ) rand = rand+1.0
!!$          x(i) = rand
!!$          i = i-1
!!$          if( i == 0 ) i = 17
!!$          j = j-1
!!$          if( j == 0 ) j = 17
!!$        end function fib_srnd_par

        function fib_simple_rnd()         result(rand)
          real (kind = real_8)   ::  rand,z
          integer, parameter     ::  m = 16807,a=54321
          integer, save          ::  i_ctrl = 0
          integer, save          ::  i=2, j=1
          real (kind = real_8), dimension(2)   ::  x
          if( i_ctrl == 0 )then
             do i = 1,2
                x(i) = rnd()*m
             enddo
             i_ctrl = 1
          endif
          z = modulo(x(i) + x(j)+real(a),real(m,real_8))
          rand = z / real(m)
          i = i-1
          if( i == 0 ) i = 2
          j = j-1
          if( j == 0 ) j = 2          
          x(i) = z
        end function fib_simple_rnd

        function rnd_lcg_r(a,x,m,c)       result(rand)
          real (kind=real_8), intent(in)     ::  a,m,c
          real (kind=real_8), intent(inout)  ::  x
          real (kind = real_8)   ::  rand
          x = modulo(a*x+c,m)
          rand = x/m
        end function rnd_lcg_r

        function rnd_lcg_i(a,x,m,c)       result(rand)
          integer (kind=integer_8), intent(in)     ::  a,m,c
          integer (kind=integer_8), intent(inout)  ::  x
          real (kind = real_8)   ::  rand
          x = modulo(a*x+c,m)
          rand = real(x)/real(m)
        end function rnd_lcg_i

        FUNCTION gasdev()    result( rand )
          REAL (kind = real_8)  ::  rand
          INTEGER, save               ::  iset
          REAL (kind = real_8), save  ::  gset
          REAL (kind = real_8)        ::  fac,rsq,v1,v2
          DATA iset/0/
          if (iset.eq.0) then
             do 
                v1=2.0*rnd()-1.0
                v2=2.0*rnd()-1.0
                rsq=v1**2+v2**2
                if(rsq < 1.0 .and. rsq /= 0.0 ) exit
             enddo
             fac=sqrt(-2.*log(rsq)/rsq)
             gset=v1*fac
             rand=v2*fac
             iset=1
          else
             rand=gset
             iset=0
          endif
        END FUNCTION gasdev

        FUNCTION sgasdev(seed)    result( rand )
          REAL (kind = real_8)  ::  rand
          INTEGER, save               ::  iset
          REAL (kind = real_8), intent(inout)  ::  seed
          REAL (kind = real_8), save  ::  gset
          REAL (kind = real_8)        ::  fac,rsq,v1,v2
          DATA iset/0/
          if (iset.eq.0) then
             do 
                v1=2.0d0*fib_srnd(seed)-1.0d0
                v2=2.0d0*fib_srnd(seed)-1.0d0
                rsq=v1**2+v2**2
                if(rsq < 1.0d0 .and. rsq /= 0.0d0 ) exit
             enddo
             fac=sqrt(-2.d0*log(rsq)/rsq)
             gset=v1*fac
             rand=v2*fac
             iset=1
          else
             rand=gset
             iset=0
          endif
        END FUNCTION sgasdev

        function beta_distrib(x,a1,a2)     result( beta )
          real (kind = real_8), intent(in)    ::  x
          integer, intent(in)                 ::  a1,a2
          real (kind = real_8)                ::  beta,x1,x2
          x1 = real(a1,real_8)
          x2 = real(a2,real_8)
          beta = exp(gammln(x1+x2))/(exp(gammln(x1))*exp(gammln(x2)))*(x**(a1-1)*(1.0-x)**(a2-1))
        end function beta_distrib

        function cauchy_fct(x,x0,gamma)  result(f)
          real (kind = real_8), intent(in)    ::  x,x0,gamma
          real (kind = real_8)                ::  f
          f = 0.5*gamma/(pi*((x-x0)**2+(0.5*gamma)**2))
        end function cauchy_fct

        function cauchy_distrib()  result(x)
!!$        function cauchy_distrib( x0,gamma )  result(x)
!!$          real (kind = real_8), intent(in)  ::  x0,gamma
          real (kind = real_8)              ::  x0,gamma
          real (kind = real_8)              ::  x,u1,u2,c
          real (kind = real_8), parameter   ::  a = 100.0
          x0=0.0
          gamma=1.0
          c = 2.0/(pi*gamma)
          do 
             u1 = a*(2.0*rnd()-1.0)
             u2 = rnd()
             x = cauchy_fct( u1,x0,gamma )
             if( x >= c*u2 ) exit
          enddo
          x = u1
        end function cauchy_distrib

        function bernoulli_distrib(q)    result(f)
          ! gives discrete numbers p(-1)=1-q and p(1)=q
          real (kind=real_8), intent(in)  ::  q
          real (kind=real_8)              ::  z,f
          z = rnd()
          if( z == 1.0 ) z = 0.9999999
          f = 1.0 - 2.0*int(z+1.0-q)
        end function bernoulli_distrib

        function uniform_distrib()    result(f)
          real (kind=real_8)              ::  f
          f = 2.0*rnd()-1.0
        end function uniform_distrib

        function gauss_distrib()    result(f)
          real (kind=real_8)              ::  f
          f = gasdev()
        end function gauss_distrib
 

        function gauss(x,sigma,mu)      result(f)
          real (kind=real_8), intent(in)  ::  x,sigma,mu
          real (kind=real_8)              ::  f,arg
          arg = 0.5*(x-mu)**2/sigma**2
          if( arg < 30.0 )then
             f = exp(-0.5*(x-mu)**2/sigma**2)/(sqrt(2.0*pi)*sigma)          
          else
             f = 0.0
          endif
        end function gauss

        function log_norm_fct(x,sigma,mu)      result(f)
          real (kind=real_8), intent(in)  ::  x,sigma,mu
          real (kind=real_8)              ::  f,f1,f2
          f1 = exp(-0.5*(log(x)-mu)**2/sigma**2)
          f2 = (sqrt(2.0*pi)*sigma*x)          
          f = f1/f2
          if( f < 1.0e-15 ) f = 0.0
        end function log_norm_fct


        subroutine norm_hist_1( x,dx )
          real (kind = real_8), intent(inout),dimension(:)  ::  x
          real (kind = real_8), intent(in)     ::  dx
          real (kind = real_8)                 ::  sum
          integer                              ::  i
          sum = 0.0
          do i = 1,size(x)  !-size(x)/2,size(x)/2
             sum = sum + x(i)
          enddo
          x = x/(sum*dx)
        end subroutine norm_hist_1

        subroutine norm_hist_2( x,dx )
          real (kind = real_8), intent(inout),dimension(:,:)  ::  x
          real (kind = real_8), intent(in)     ::  dx
          real (kind = real_8)                 ::  sum
          integer                              ::  i,j
          do j = 1,size(x,1)
             sum = 0.0
             do i = 1,size(x,2)  !-size(x)/2,size(x)/2
                sum = sum + x(j,i)
             enddo
             x(j,:) = x(j,:)/(sum*dx)
          enddo
        end subroutine norm_hist_2

        subroutine norm_hist_3( x,dx )
          real (kind = real_8), intent(inout),dimension(:,:,:)  ::  x
          real (kind = real_8), intent(in)     ::  dx
          real (kind = real_8)                 ::  sum
          integer                              ::  i,j,k
          do k = 1,size(x,1)
             do j = 1,size(x,2)
                sum = 0.0
                do i = 1,size(x,3)  !-size(x)/2,size(x)/2
                   sum = sum + x(k,j,i)
                enddo
                x(k,j,:) = x(k,j,:)/(sum*dx)
             enddo
          enddo
        end subroutine norm_hist_3

        function cum_norm_inv(u)    result(x)
          real (kind = real_8), intent(in)     ::  u
          real (kind = real_8)                 ::  x,y,r
          real (kind = real_8), parameter  ::   &
               a0=2.50662823884, a1=-18.61500062529, a2=41.39119773534, a3=-25.44106049637, &
               b0=-8.47351093090, b1=23.08336743743, b2=-21.06224101826, b3=3.13082909833,  &
               c0=0.3374754822726147, c1=0.9761690190917186, c2=0.1607979714918209,    &
                      c3=0.0276438810333863, c4=0.0038405729373609, c5=0.0003951896511919,    &
                      c6=0.0000321767881768, c7=0.0000002888167364, c8=0.0000003960315187
          
          y = u-0.5
          if( abs(y) < 0.42 )then
             r = y**2
             x = y * (((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1.0)
          else
             r = u
             if( y > 0.0 ) r = 1.0-u
             r = log(-log(r))
             x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))))
             if( y < 0.0 ) x = -x
          endif
        end function cum_norm_inv

        function cum_norm(x)   result(y)
          real (kind = real_8), intent(in)     ::  x
          real (kind = real_8)                 ::  y,a,t,s
          real (kind = real_8), parameter  ::   &
               b1 = 0.31938153, b2 = -0.356563782, b3 = 1.781477937,  &
               b4 = -1.821255978, b5 = 1.330274429,  &
               p = 0.2316419, c=0.9189385333204672
          a = abs(x)
          t = 1.0/(1.0+a*p)
          s = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t
          y = s*exp(-0.5*x**2-c)
          if(x>0.0) y = 1.0-y
        end function cum_norm

        FUNCTION gammln(xx)          result(gamma)
          REAL (kind = real_8), intent(in) ::  xx
          REAL (kind = real_8)             ::  gamma
          INTEGER j
          REAL (kind = real_8)       ::  ser,tmp,x,y
          REAL (kind = real_8), save ::  stp
          REAL (kind = real_8), dimension(6), save ::  cof
          DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,       &
               24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,   &
               -.5395239384953d-5,2.5066282746310005d0/
          x=xx
          y=x
          tmp=x+5.5d0
          tmp=(x+0.5d0)*log(tmp)-tmp
          ser=1.000000000190015d0
          do j=1,6
             y=y+1.d0
             ser=ser+cof(j)/y
          enddo
          gamma=tmp+log(stp*ser/x)
        END FUNCTION gammln

        function levy(x,a,mu)      result(f)
          real (kind = real_8), intent(in)   ::  x,mu,a
          real (kind = real_8)               ::  f
          integer, parameter                 ::  n_max = 10
          integer                            ::  n
          f = 0.0
          do n = 1, n_max
             f = f + (-1.0)**(n+1)*a**n *gammln(1.0+n*mu)*sin(0.5*pi*mu*n)  /  &
                  (pi*factorial(n)*x**(1.0+n*mu))
          enddo
        end function levy

        function triangle_fct_x_sig_mu(x,sig,mu)       result(f)
          real (kind = real_8), intent(in)  ::  x,sig,mu
          real (kind = real_8)              ::  f

          if(abs(x) > 2.0*sig )then
             f = 0.0
          else if( x <= 0.0 )then
             f = x/sqrt(2.0*sig) + sqrt(2.0*sig)
             f = f/sqrt(8.0*sig**3)
          else
             f = -x/sqrt(2.0*sig) + sqrt(2.0*sig)
             f = f/sqrt(8.0*sig**3)
          endif
        end function triangle_fct_x_sig_mu

        function triangle_fct_xxxy(x,x1,x2,x3,y)       result(f)
          real (kind = real_8), intent(in)  ::  x,x1,x2,x3,y
          real (kind = real_8)              ::  f

          if( x < x1 .or. x > x3 )then
             f = 0.0
          else if( x >= x1 .and. x <= x2 )then
             f = y*(x-x1)/(x2-x1)
          else
             f = y * (x3-x)/(x3-x2)
          endif
        end function triangle_fct_xxxy

        function triang_distrib_sig_mu(sigma,mu)  result(x)
          real (kind = real_8), intent(in)  ::  sigma,mu
          real (kind = real_8)              ::  x,u1,u2,c
          c = 0.5/(sigma)
          do 
             u1 = 2.0*sigma*(2.0*rnd()-1.0)
             u2 = rnd()
             x = triangle_fct( u1,sigma,mu )
             if( x >= c*u2 ) exit
          enddo
          x = u1
        end function triang_distrib_sig_mu

        function triang_distrib_xxxy(x0,x1,x2,y)  result(x)
          real (kind = real_8), intent(in)  ::  x0,x1,x2,y
          real (kind = real_8)              ::  x,u1,u2,c
          do 
             u1 = x0 + rnd()*(x2-x0)
             u2 = rnd()
             x = triangle_fct( u1,x0,x1,x2,y )
             if( x >= y*u2 ) exit
          enddo
          x = u1
        end function triang_distrib_xxxy


        function empirical_distrib()      result(xi)

          integer, save  ::  i_ctrl = 0
        integer, parameter       ::  n_hist = 1000, n_rand = 10000
        real (kind = real_8)     ::  xi
        real (kind = real_8), dimension(-n_hist:n_hist)       ::  hist, hist_rej
        real (kind = real_8)     ::  x_emp_min, x_emp_max, p_emp_max, hist_max,var,x_avg,x_sq_avg,sig
        real (kind = real_8),save ::  x,y,z,xd,u1,u2,dx,sum1,sum2,p_max,sigma,sigma_inv,mu,c
                                     
        real (kind = real_8), parameter                       ::  de = 0.001
        real (kind = real_8), allocatable, dimension(:),save  ::  p_emp,x_emp
        real (kind = real_8), allocatable, dimension(:,:)     ::  tmp2
        integer                  ::  i,ii,j,k,n,n0,n_min,n_max,n_sav,n_size
        integer                  ::  ios_read, n_rej
          
          select case( i_ctrl )
          case default
             do
                u1 = triang_distrib(sigma,mu)
                u2 = rnd()
                i = nint(u1/dx)
                x = p_emp(i)
                if( x > c*u2 )then
                   exit
                else
                   cycle
                endif
             enddo
             xi = (u1-x_avg) * sigma_inv

             case( 0 )
                ! initialization step
                open( unit=99,file="rand_empirical.his", status="unknown", action="read" )
                n = 0
                x_emp_min = huge(x_emp_min)
                x_emp_max = -huge(x_emp_max)
                do 
                   read(unit=99,fmt=*,iostat=ios_read) xd,xd,x,z
                   if(z > 0.0 )then
                      if( x < x_emp_min )then
                         x_emp_min = x                 
                         n_min = n
                      endif
                      if( x > x_emp_max )then
                         x_emp_max = x
                         n_max = n
                      endif
                      if( x == 0.0 ) n0 = n
                   endif
                   if( ios_read /= 0 ) exit
                   n = n+1
                enddo
                rewind(99)
                n_sav = n
                print*," << Empirical distribution setup >>"
                print*," found ",n_sav," input lines"
                print*," minimum x: ",n_min,x_emp_min
                print*," maximum x: ",n_max,x_emp_max

                n_size = max(n_max-n0,n0-n_min)
                print*," size of hist: ",n
                allocate( tmp2(2,1:n_sav))
                allocate( p_emp(-n:n))
                allocate( x_emp(-n:n))
                do i = 1,n_sav
                   read(unit=99,fmt=*,iostat=ios_read) xd,xd,tmp2(1:2,i)
                enddo
                close(99)
                n = 0
                do i = 1,n_sav
                   if(tmp2(1,i) >= 0.0 )then
                      x_emp(n) = tmp2(1,i)
                      p_emp(n) = tmp2(2,i)
                      n = n+1
                      if( n > n_size ) exit
                   endif
                enddo
                n = 0
                do i = n_sav,1,-1
                   if(tmp2(1,i) < 0.0 )then
                      n = n-1
                      x_emp(n) = tmp2(1,i)
                      p_emp(n) = tmp2(2,i)
                      if( n < -n_size ) exit
                   endif
                enddo

                p_max = 0.0
                do i = -n_size,n_size
                   if( p_emp(i) > p_max ) p_max = p_emp(i)
                enddo
                dx = x_emp(1)-x_emp(0)
                p_emp = p_emp/(sum(p_emp)*dx)

                ! calculation of moments
                x_avg = 0.0
                x_sq_avg = 0.0
                var = 0.0
                do i=-n_size,n_size
                   x_avg = x_avg + x_emp(i)*p_emp(i)*dx
                   x_sq_avg = x_sq_avg + x_emp(i)**2*p_emp(i)*dx
                enddo
                var = x_sq_avg - x_avg**2
                sigma_inv = 1.0/sqrt(var)

                print*," avg. value: ",x_avg
                print*," variance:   ",var
                print*," standard deviation: ",sqrt(var)

                ! empirical distribution with acceptance-rejection
                open(unit=99,file="rand_empirical.dat",status="unknown",action="write")
                sigma = 0.5*x_emp(n_size)
                mu = 0.0
                hist = 0.0
                hist_rej = 0.0
                c = p_max 
                print*," triang: c = ",c
                n = 0
                n_rej = 0
                do 
                   u1 = triang_distrib(sigma,mu)
                   u2 = rnd()
                   i = nint(u1/dx)
                   x = p_emp(i)
                   if( x <= c*u2 )then
                      n_rej = n_rej+1
                      j = nint(u1/dx)
                      if( iabs(j) > n_hist ) cycle
                      hist_rej(j) = hist_rej(j)+1.0
                      cycle
                   endif
                   j = nint(u1/dx)
                   if( iabs(j) > n_hist ) cycle
                   hist(j) = hist(j)+1.0
                   n = n+1
                   if( n >= n_rand ) exit
                enddo
                sum1 = 0.0
                sum2 = 0.0
                do i = -n_hist,n_hist           
                   sum1 = sum1 + hist(i)*dx
                   sum2 = sum2 + hist_rej(i)*dx
                enddo
                hist(-n_hist) = 2.0*hist(-n_hist)
                hist(n_hist) = 2.0*hist(n_hist)
                hist_rej(-n_hist) = 2.0*hist_rej(-n_hist)
                hist_rej(n_hist) = 2.0*hist_rej(n_hist)
                print*," Triang: ratio acc/rej = ",real(n_rand)/real(n_rej),sum1/sum2
!!$        call norm_hist( hist,dx )
!!$        call norm_hist( hist_rej,dx )
                x = sum(hist)*dx
                hist = hist/x
                hist_rej = hist_rej/x
                sum1 = 0.0

                p_emp_max = 0.0
                hist_max  = 0.0
                do i = -n_hist,n_hist
                   if( hist(i) > hist_max ) hist_max = hist(i)
                   j = nint(x/dx)
                   if( j >= -n_sav .and. j <= n_sav )then
                      if( p_emp(j) > p_emp_max ) p_emp_max = p_emp(j)
                   endif
                enddo
                hist_max = max(hist_max,p_emp_max)

                sig = sqrt(var)
                mu =0.0

                do i = -n_hist,n_hist
                   x = real(i)*dx
                   j = nint(x/dx)
                   if( j >= -n_sav .and. j <= n_sav )then
                      y = p_emp(j)
                   else
                      y = 0.0
                   endif
                   write(91,"(10e14.6)") x,hist(i),gauss(x,sig,mu),y,hist_rej(i),   &
                                         2.0*sigma*hist_max*triangle_fct(x,sigma,mu),sum1
                   sum1 = sum1 + hist(i)*dx
                enddo
                close(91)

             end select

             i_ctrl = 9999

           end function empirical_distrib


        function u_bsp(r,rc,order)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of given order as shape function for charges
!----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          integer, intent(in)                   ::  order
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )
          character (len=*), parameter          ::  method_name = "u_bsp"

          select case( order )
!!$          case( -1 )
!!$             f = u_cos0(r,rc)

          case( 0 )
             f = u_bsp0(r,rc)

          case( 1 )
             f = u_bsp1(r,rc)

          case( 2 )
             f = u_bsp2(r,rc)

          case( 3 )
             f = u_bsp3(r,rc)

          case( 4 )
             f = u_bsp4(r,rc)

          case default
             write(6,"(//,a,i3,a,//,a,//)") module_name//"::"//method_name//": Order ",  &
                  order," spline not yet implemented as source term",    &
                  " >>>>>  program stops  <<<<<"

          end select

        end function u_bsp


        function u_bsp0(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc)then
             f = (1.5d0*rc**2 - 0.5d0*r**2) / rc**3

          else if( r > rc )then
             f = 1.d0 / r

          endif

        end function u_bsp0

        function u_bsp1(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc)then
             f = 12.d0*r**2 * (1.d0/3.d0-(1.d0/4.d0)*r/rc)/rc**3   &
                  + 12.d0*(1.d0/2.d0*(1.d0-r**2/rc**2)              &
                  - 1.d0/3.d0*(1.d0-r**3/rc**3))/rc
          else if( r > rc )then
             f = 1.d0 / r
          endif

        end function u_bsp1

        function u_bsp2(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc/3.d0)then
             f = (243.d0/80.d0)*r**4/rc**5 - (27.d0/8.d0)*r**2/rc**3   &
                  + 39.d0/(16.d0*rc)
          else if( r > rc/3.d0 .and. r <= rc )then
             f = -1.d0/(80.d0*r)-(243.d0/160.d0)*r**4/rc**5                  &
                  + (81.d0/16.d0)*r**3/rc**4-(81.d0/16.d0)*r**2/rc**3      &
                  + 81.d0/(32.d0*rc)
          else if( r > rc )then
             f = 1.d0 / r
          endif

        end function u_bsp2


        function u_bsp3(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc/2.d0)then
             f = -(32.d0/5.d0)*r**5/rc**6 + (48.d0/5.d0)*r**4/rc**5   &
                  - (16.d0/3.d0)*r**2/rc**3 + 14.d0/(5.d0*rc)
          else if( r > rc/2.d0 .and. r <= rc )then
             f = -1.d0/(15.d0*r) + (32.d0/15.d0)*r**5/rc**6           &
                  - (48.d0/5.d0)*r**4/rc**5 + 16.d0*r**3/rc**4        &
                  - (32.d0/3.d0)*r**2/rc**3 + 16.d0/(5.d0*rc)
          else if( r > rc )then
             f = 1.d0 / r
          endif

        end function u_bsp3


        function u_bsp4(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 4th order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc/5.d0)then
             f = -(15625.d0/896.d0)*r**6/rc**7 + (1875.d0/128.d0)*r**4/rc**5  &
                  - (2875.d0/384.d0)*r**2/rc**3 + 1199.d0/(384.d0*rc)
          else if( r > rc/5.d0 .and. r <= 3.d0*rc/5.d0 )then
             f = 1.d0/(6720.d0*r) + (15625.d0/1344.d0)*r**6/rc**7             &
                  - (3125.d0/96.d0)*r**5/rc**6 + (1875.d0/64.d0)*r**4/rc**5   &
                  - (625.d0/192.d0)*r**3/rc**4 - (1375.d0/192.d0)*r**2/rc**3  &
                  + 599.d0/(192.d0*rc)
          else if( r > 3.d0*rc/5.d0 .and. r <= rc )then
             f = -437.d0/(2688.d0*r) - (15625.d0/5376.d0)*r**6/rc**7 +        &
                  (3125.d0/192.d0)*r**5/rc**6 - (9375.d0/256.d0)*r**4/rc**5   &
                  + (15625.d0/384.d0)*r**3/rc**4                              &
                  - (15625.d0/768.d0)*r**2/rc**3 + 3125.d0/(768.d0*rc)
          else if( r > rc )then
             f = 1.d0 / r
          endif

        end function u_bsp4




        function e_bsp(r,rc,order)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of given order as shape function for charges
!----------------------------------------------------------------------------------

          real (kind = real_8), intent(in)    ::  r, rc
          integer, intent(in)                   ::  order
          real (kind = real_8)                ::  f
          real (kind = real_8), parameter     ::  pi = acos( -1.d0 )
          character (len=*), parameter          ::  method_name = "u_bsp"

          select case( order )
!!$          case( -1 )
!!$             f = u_cos0(r,rc)

          case( 0 )
             f = e_bsp0(r,rc)

          case( 1 )
             f = e_bsp1(r,rc)

          case( 2 )
             f = e_bsp2(r,rc)

          case( 3 )
             f = e_bsp3(r,rc)

          case( 4 )
             f = e_bsp4(r,rc)

          case default
             write(6,"(//,a,i3,a,//,a,//)") module_name//"::"//method_name//": Order ",  &
                  order," spline not yet implemented as source term",    &
                  " >>>>>  program stops  <<<<<"

          end select

        end function e_bsp


        function e_bsp0(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r < rc )then
             f = r / R**3
                          
          else if( r > rc )then
             f = 1.d0 / r**2
             
          endif

        end function e_bsp0

        function e_bsp1(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r < rc )then
             f = r * (4.d0*rc - 3.d0*r) / rc**4
                          
          else if( r > rc )then
             f = 1.d0 / r**2
             
          endif

        end function e_bsp1
        

        function e_bsp2(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r <= rc/3.d0)then
             f = 27.d0*r*( 5.d0*rc**2 - 9.d0*r**2 ) / (20.d0*rc**5)
             
          else if( r > rc/3.d0 .and. r <= rc )then
             f = - ( rc**5 - 810.d0*rc**2 * r**3 + 1215.d0 *r**4 * rc - 486.d0*r**5 ) &
                  / ( 80.d0*r**2 * rc**5)
             
          else if( r > rc )then
             f = 1.d0 / r**2
             
          endif

        end function e_bsp2

        
        function e_bsp3(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 3rd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r < rc/2.d0)then
             f = ( 160.d0*rc**3*r**3 - 576.d0*rc*r**5 + 480.d0*r**6 ) / ( 15.d0*r**2*rc**6 )
             
          else if( r >= rc/2.d0 .and. r < rc )then
             f = ( -rc**6 + 320.d0*rc**3*r**3 - 720.d0*rc**2*r**4 + 576.d0*rc*r**5 - 160.d0*r**6)  &
                  / ( 15.d0*r**2*rc**6 )
             
          else if( r > rc )then
             f = 1.d0 / r**2
             
          endif

        end function e_bsp3
        

        function e_bsp4(r,rc)       result(f)
!----------------------------------------------------------------------------------
!  B-Spline of 2nd order as shape function for charges
!----------------------------------------------------------------------------------
          
          real (kind = real_8), intent(in)    ::  r, rc
          real (kind = real_8)                ::  f

          if( r >= 0.d0 .and. r < rc/5.d0)then
             f = -( -40250.d0 * rc**4 * r + 157500.d0 * rc**2 * r**3 - 281250.d0 * r**5 ) &
                  / (2688.d0 * rc**7)
             
          else if( r >= rc/5.d0 .and. r < 3.d0*rc/5.d0 )then
             f = - (20965.d0*rc**6 - 144375.d0*rc**4*r**2 - 87500.d0*rc**3*r**3   &
                  + 984375.d0*rc**2*r**4 - 1312500.d0*rc*r**5 + 546875.d0*r**6)   &
                  / (6720.d0*rc**7*r) &
                  + (rc**7 + 20965.d0*rc**6*r - 48125.d0*rc**4*r**3 - 21875.d0*rc**3*r**4  &
                  + 196875.d0*rc**2*r**5 - 218750.d0*rc*r**6 + 78125.d0*r**7)  &
                  / (6720.d0*rc**7*r**2)
             
          else if( r >= 3.d0*rc/5.d0 .and. r < rc )then
             f = -(21875.d0*rc**6 - 328125.d0*rc**4*r**2 + 875000.d0*rc**3*r**3   &
                  - 984375.d0*rc**2*r**4 + 525000.d0*rc*r**5 - 109375.d0*r**6)    &
                  / (5376.d0*rc**7*r)      &
                  + (-874.d0*rc**7 + 21875.d0*rc**6*r - 109375.d0*rc**4*r**3      &
                  + 218750.d0*rc**3*r**4 - 196875.d0*rc**2*r**5 + 87500.d0*rc*r**6  &
                  - 15625.d0*r**7)    &
                  / (5376.d0*rc**7*r**2)
             
          else if( r >= rc )then
             f = 1.d0 / r**2
             
          endif

        end function e_bsp4








!!$
!!$        function bsp0(x,x0,R)     result(f)
!!$          real (kind=real_8), intent(in)  ::  x,x0,R
!!$          real (kind=real_8)              ::  f,z
!!$          z = abs(x - x0)
!!$          if( z < R )then
!!$             f = 1.0/(2.0*R)
!!$          else
!!$             f = 0.0d0
!!$          endif
!!$        end function bsp0
!!$
!!$        function bsp1(x,x0,R)     result(f)
!!$          real (kind=real_8), intent(in)  ::  x,x0,R
!!$          real (kind=real_8)              ::  f,z
!!$          z = abs(x - x0)
!!$          if( z < R )then
!!$             f = (1.0 - z/R )/R
!!$          else
!!$             f = 0.0d0
!!$          endif
!!$        end function bsp1
!!$
!!$        function bsp2(x,x0,R)     result(f)
!!$          real (kind=real_8), intent(in)  ::  x,x0,R
!!$          real (kind=real_8)              ::  f,z
!!$          z = abs(x - x0)
!!$          if( z < R/3.0d0 )then
!!$             f = 9.0/4.0 - 27.0/4.0*z**2/R**2
!!$          else if( z <= R )then
!!$             f = 27.0d0/8.0d0 - 27.0/4.0*z/R + 27.0d0/8.0d0 * z**2/R**2 
!!$          else

!!$             f = 0.0d0
!!$          endif
!!$          f = f*1.0/(2.0*R)
!!$        end function bsp2
!!$
!!$        function bsp3(x,x0,R)     result(f)
!!$          real (kind=real_8), intent(in)  ::  x,x0,R
!!$          real (kind=real_8)              ::  f,z
!!$          z = abs(x - x0)
!!$          if( z < 0.5d0*R )then
!!$             f = 4.0 - 24.0*z**2/R**2 + 24.0*z**3/R**3
!!$          else if( z <= R )then
!!$             f = 8.0d0 - 24.0*z/R + 24.0d0 * z**2/R**2 - 8.0d0 * z**3/R**3
!!$          else
!!$             f = 0.0d0
!!$          endif
!!$          f = f*1.0/(3.0*R)
!!$        end function bsp3
!!$
!!$        function bsp4(x,x0,R)     result(f)
!!$          real (kind=real_8), intent(in)  ::  x,x0,R
!!$          real (kind=real_8)              ::  f,z
!!$          z = abs(x - x0)
!!$          if( z < R/5.0d0 )then
!!$             f = 5.0*(23.0*R**4 - 150.0*R**2*z**2 + 375.0*z**4) / (192.0*R**4)
!!$          else if( z <= 3.0*R/5.0 )then
!!$             f = 5.0*(11.0*R**4 + 10.0*R**3*z - 150.0*R**2*z**2 + 250.0*R*z**3 - 125.0*z**4) &
!!$                  / (96.0*R**4)
!!$          else if( z <= R )then
!!$             f = 1.0/384.0 * (5.0-5.0*z/R)**4
!!$          else
!!$             f = 0.0d0
!!$          endif
!!$          f = f*5.0/(2.0*R)
!!$        end function bsp4


        function pot_bsp1(x,r)   result(f)
!-----------------------------------------------------------------------
! potential of 1st order B-Spline
!-----------------------------------------------------------------------
          real (kind=real_8), intent(in)  ::  x,r
          real (kind=real_8)              ::  f,z

          z = abs(x)
          ! potential of 1st order B-Spline
          if (z > r) then
             f = 1.0d0 / z
          else 
             f = (2.0*r**3 - 2.0*r*z**2 + z**3) / r**4
          end if
        end function pot_bsp1

        function pot_bsp2(x,r)   result(f)
!-----------------------------------------------------------------------
! potential of 2nd order B-Spline
!-----------------------------------------------------------------------
          real (kind=real_8), intent(in)  ::  x,r
          real (kind=real_8)              ::  f,z

          z = abs(x)
          ! potential of 2nd order B-Spline
          if (z > r) then
             f = 1.0d0 / z
          elseif ( z <= r .and. z > r/3.0d0) then
             f = ( -243.0d0 * z**5                 &
                  + 810.0d0 * z**4 * r             &
                  - 810.0d0 * z**3 * r**2          &
                  + 405.0d0 * z    * r**4          &
                  - 2.0d0*r**5) / (160.0d0 * z * r**5)
          else 
             f = 3.0d0 * (81.0d0*z**4 - 90d0*z**2*r**2         &
                  + 65.0d0 * r**4) / (80.0d0 * r**5)
          end if
        end function pot_bsp2

        function pot_bsp3(x,r)   result(f)
!-----------------------------------------------------------------------
! potential of 3rd order B-Spline
!-----------------------------------------------------------------------

          real (kind=real_8), intent(in)  ::  x,r
          real (kind=real_8)              ::  f,z

          z = abs(x)
          ! potential of 3rd order B-Spline
          if ( z > r ) then
             f = 1.0d0 / z
          elseif ( z > 0.5d0*r ) then
             f = (19.0d0/(30.0d0*z) + (-1344.0d0*r**6 +    &
                  6144.0d0*r**5*z -                    &
                  20480.0d0*r**3*z**3 +                &
                  30720.0d0*r**2*z**4 -                &
                  18432.0d0*r*z**5 +                   &
                  4096.0d0*z**6)/(1920.0d0*r**6*z))
          else
             f = (-5760.0d0*z**5 + 8640.0d0*z**4*r       &
                  -4800.0d0*z**2*r**3 + 2520.0d0*r**5    &
                  ) / (900.0d0*r**6)
          end if
        end function pot_bsp3

        function pot_bsp4(x,r)   result(f)
!-----------------------------------------------------------------------
! potential of 4th order B-Spline
!-----------------------------------------------------------------------

          real (kind=real_8), intent(in)  ::  x,r
          real (kind=real_8)              ::  f,z

          z = abs(x)
          ! potential of 3rd order B-Spline
          if( z < r/5.0d0 )then
             f = (8393.0d0*r**6 - 20125.0d0*r**4*z**2 +            &
                  39375.0d0*r**2*z**4 - 46875.0d0*z**6)            &
                  / (2688.0d0*r**7)
          else if( z < 3.0d0*r/5.0d0 )then
             f = (r**7 + 20965.0d0*r**6*z - 48125.0d0*r**4*z**3 - 21875.0d0*r**3*z**4 &
                  + 196875.0d0*r**2*z**5 - 218750.0d0*r*z**6 + 78125.0d0*z**7) &
                  / (6720.0d0*r**7*z)
          else if( z < r )then
             f = -(874.0d0*r**7 - 21875.0d0*r**6*z + 109375.0d0*r**4*z**3 - 218750.0d0*r**3*z**4 &
                  + 196875.0d0*r**2*z**5 - 87500.0d0*r*z**6 + 15625.0d0*z**7  ) &
                  / ( 5376.0d0*r**7*z)
          else
             f = 1.0d0 / z
          end if
        end function pot_bsp4
           
      end module functions_module





