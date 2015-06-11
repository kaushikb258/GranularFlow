            ! KAUSHIK CREATED THESE SUBROUTINES
!----------------------------------------------------------

           subroutine compute_g0(alps,g0)

          implicit none

            real(kind=8) :: alps, g0, t, alps1
            real(kind=8) :: alps_max, erest, rhos, dias 


            call solid_constants(alps_max,erest,rhos,dias)

 
            ! make sure g0 doesn't blow up!
            alps1 = max(min(alps,0.65d0),0.0d0)
 

            t = 1.0d0 - (alps1/alps_max)**0.333d0   

            g0 = 1.0d0/t


           return
           end subroutine

    
!----------------------------------------------------------

           subroutine compute_cs(alps,thetas,cs)

           implicit none

           real(kind=8) :: alps, thetas, cs, g0, alps1 
           real(kind=8) :: t1, t2, t3, t4
           real(kind=8) :: alps_max, erest, rhos, dias 

            call solid_constants(alps_max,erest,rhos,dias)
           
           call compute_g0(alps,g0) 


           ! make sure g0 doesn't blow up!
            alps1 = max(min(alps,0.65d0),0.0d0)
 


           t1 = 1.0d0 + 4.0d0*(1.0d0 + erest)*g0*alps1
           t1 = t1*thetas

           t2 = 2.0d0*(1.0d0+erest)*thetas*g0*g0/3.0d0* & 
                    (alps1**1.33333d0)/(alps_max**0.3333d0) 

           t3 = 1.0d0 + 2.0d0*(1.0d0+erest)*g0*alps1
           t4 = 2.0d0/3.0d0*thetas*t3*t3 
               
           cs = t1 + t2 + t4

           if(cs.lt.0.0d0) then
             write(*,*) 'error in cs ', cs, alps, thetas, alps1 
             !call ps_stop('solid routines 2 ') 
           endif

           cs = sqrt(cs)

           return
           end subroutine 

!----------------------------------------------------------

      subroutine cons_2_prim(N0,N5,M0,M5,K0,K5,NV,gamma,rhos, &
          Q,rhog,ug,vg,wg,eg,pg,alps,alpg,us,vs,ws,thetas,es)

      implicit none   

      real(kind=8), parameter :: alps_cut = 1.0d-5 

      integer :: N0, N5, M0, M5, K0, K5, NV, i, j, k
      real(kind=8) :: gamma, rhos, ke
          
      real(kind=8), dimension(N0:N5,M0:M5,K0:K5,1:NV) :: Q  
      real(kind=8), dimension(N0:N5,M0:M5,K0:K5) :: rhog, ug, vg, wg, eg, pg
      real(kind=8), dimension(N0:N5,M0:M5,K0:K5) :: alps, alpg
      real(kind=8), dimension(N0:N5,M0:M5,K0:K5) :: us, vs, ws, thetas, es


      do i = N0, N5  
      do j = M0, M5  
      do k = K0, K5  
        
       alps(i,j,k) = Q(i,j,k,7)/rhos


               if(alps(i,j,k).lt.0.0d0.or.alps(i,j,k).gt.1.0d0) then
                !write(*,*) 'alps bug ', alps(i,j,k), rhos, i, j, k
                !call ps_stop()
                alps(i,j,k) = 0.0d0
               endif


       alps(i,j,k) = min(max(alps(i,j,k),0.0d0),1.0d0)
       !alps(i,j,k) = min(alps(i,j,k),0.7d0)
       Q(i,j,k,7) = alps(i,j,k)*rhos

       alpg(i,j,k) = 1.0d0 - alps(i,j,k)  

       rhog(i,j,k) = Q(i,j,k,1)/alpg(i,j,k) 
       ug(i,j,k) = Q(i,j,k,2)/Q(i,j,k,1) 
       vg(i,j,k) = Q(i,j,k,3)/Q(i,j,k,1) 
       wg(i,j,k) = Q(i,j,k,4)/Q(i,j,k,1) 

       ke = 0.5d0*(ug(i,j,k)**2.0d0 + vg(i,j,k)**2.0d0 + wg(i,j,k)**2.0d0)
       eg(i,j,k) = Q(i,j,k,5)/Q(i,j,k,1) - ke 

       if(eg(i,j,k).le.0.0d0) then
!        write(*,*) 'eg < 0 in cons_2_prim ', i,j,k,eg(i,j,k),  & 
!                   rhog(i,j,k), alpg(i,j,k), ke, Q(i,j,k,:)

          ! eg reference!
          eg(i,j,k) = 7.5d5 
             !call ps_stop('solid routines 3 ') 
       endif  

       pg(i,j,k) = (gamma-1.0d0)*rhog(i,j,k)*eg(i,j,k)  
         

       if(alps(i,j,k).gt.alps_cut) then
        us(i,j,k) = Q(i,j,k,8)/Q(i,j,k,7)
        vs(i,j,k) = Q(i,j,k,9)/Q(i,j,k,7)
        ws(i,j,k) = Q(i,j,k,10)/Q(i,j,k,7)
        thetas(i,j,k) = Q(i,j,k,11)/Q(i,j,k,7)*2.0d0/3.0d0
        es(i,j,k) = Q(i,j,k,12)/Q(i,j,k,7) 
        
        thetas(i,j,k) = max(thetas(i,j,k),0.0d0)
        es(i,j,k) = max(es(i,j,k),0.0d0)

       else
        us(i,j,k) = 0.0d0
        vs(i,j,k) = 0.0d0
        ws(i,j,k) = 0.0d0
        thetas(i,j,k) = 0.0d0
        es(i,j,k) = 0.0d0 
       endif


      enddo
      enddo
      enddo


       return
       end subroutine

!-------------------------------------------------------------    
       subroutine solid_constants(alps_max,erest,rhos,dias)


           real(kind=8) :: alps_max, erest, rhos, dias

           alps_max = 0.7405d0
           erest = 0.8d0

           ! Moon
           !rhos = 3100.0d0
           !dias = 70.0d-6 

           ! Mars 
           rhos = 1650.0d0 
           dias = 400.0d-6  


       return
       end subroutine

!-------------------------------------------------------------   
