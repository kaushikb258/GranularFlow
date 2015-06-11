!------------------------------------------------------------------
!                  KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!------------------------------------------------------------------


       subroutine nozzling_terms(dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,alps, & 
                     P,U,V,W,dflux)

        implicit none

        integer :: N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: dflux
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: alps, U, V, W, P
        real(kind=8) :: dx, dy, dz
        integer :: ii, jj, kk, l, m  
        real(kind=8), dimension(1:3) :: grad_alphas, grad_alphag, gradP
        real(kind=8) :: N_x, N_y, N_z, t1  

         dflux = 0.0d0



         do ii = N2-2, N3+2 
          do jj = M2-2, M3+2
           do kk = K2-2, K3+2


             grad_alphas(1) = (alps(ii-2,jj,kk) - 8.0d0*alps(ii-1,jj,kk) + & 
                        8.0d0*alps(ii+1,jj,kk) - alps(ii+2,jj,kk))/(12.0d0*dx)
              
             grad_alphas(2) = (alps(ii,jj-2,kk) - 8.0d0*alps(ii,jj-1,kk) + & 
                        8.0d0*alps(ii,jj+1,kk) - alps(ii,jj+2,kk))/(12.0d0*dy)
              
             grad_alphas(3) = (alps(ii,jj,kk-2) - 8.0d0*alps(ii,jj,kk-1) + & 
                        8.0d0*alps(ii,jj,kk+1) - alps(ii,jj,kk+2))/(12.0d0*dz)


             grad_alphag = -grad_alphas



             ! x-mom gas 
             t1 = P(ii,jj,kk) 
             dflux(ii,jj,kk,2) = dflux(ii,jj,kk,2) + t1*grad_alphag(1) 

             ! y-mom gas 
             t1 = P(ii,jj,kk) 
             dflux(ii,jj,kk,3) = dflux(ii,jj,kk,3) + t1*grad_alphag(2) 

             ! z-mom gas 
             t1 = P(ii,jj,kk)
             dflux(ii,jj,kk,4) = dflux(ii,jj,kk,4) + t1*grad_alphag(3)

             ! energy gas
             N_x = P(ii,jj,kk)*U(ii,jj,kk) 
             N_y = P(ii,jj,kk)*V(ii,jj,kk) 
             N_z = P(ii,jj,kk)*W(ii,jj,kk) 

             dflux(ii,jj,kk,5) = dflux(ii,jj,kk,5) + N_x*grad_alphag(1) & 
                               + N_y*grad_alphag(2) + N_z*grad_alphag(3)



!------------------------------------------------------------------------


             gradP(1) = (P(ii-2,jj,kk) - 8.0d0*P(ii-1,jj,kk) + & 
                        8.0d0*P(ii+1,jj,kk) - P(ii+2,jj,kk))/(12.0d0*dx)
             gradP(2) = (P(ii,jj-2,kk) - 8.0d0*P(ii,jj-1,kk) + & 
                        8.0d0*P(ii,jj+1,kk) - P(ii,jj+2,kk))/(12.0d0*dy)
             gradP(3) = (P(ii,jj,kk-2) - 8.0d0*P(ii,jj,kk-1) + & 
                        8.0d0*P(ii,jj,kk+1) - P(ii,jj,kk+2))/(12.0d0*dz)


             dflux(ii,jj,kk,8) = dflux(ii,jj,kk,8) - alps(ii,jj,kk)*gradP(1)
             dflux(ii,jj,kk,9) = dflux(ii,jj,kk,9) - alps(ii,jj,kk)*gradP(2)
             dflux(ii,jj,kk,10) = dflux(ii,jj,kk,10) - alps(ii,jj,kk)*gradP(3)
 

           enddo 
          enddo    
         enddo
         

!-----------------------------------------
!  df/dx(i) = (f(i-2)-8f(i-1)+8f(i+1)-f(i+2))/(12*dx)
!-----------------------------------------


         return

        end subroutine 

!-------------------------------------------------------------
