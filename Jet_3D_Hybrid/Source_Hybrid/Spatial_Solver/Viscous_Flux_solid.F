!------------------------------------------------------------------
!                  KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!------------------------------------------------------------------


       subroutine viscous_flux_solid(dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,alps,U,V,W,thetas, & 
                   es,dflux)


        implicit none

        integer :: N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: dflux
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: alps, U, V, W, thetas, es
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: alpg
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: Pk, Pc
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: PP
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: qk_plus_qc

        real(kind=8) :: dx, dy, dz
        integer :: ii, jj, kk, l, m  
        real(kind=8) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(kind=8) :: g0, musc, xis, kappa_dilute
        real(kind=8) :: alps_max, erest, t1, t2, kappa_eff, rhos, dias
        real(kind=8), dimension(1:3) :: grad_theta
        REALTYPE, dimension(1:3,1:3)    :: du_dx_symm 
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3,1:3)    :: du_dx 
        REALTYPE, dimension(1:3,1:3)    :: Pc1, Pc2, Pc3 
        real(kind=8) :: taux, tauy, tauz, qx, qy, qz, visc_diss



         call solid_constants(alps_max,erest,rhos,dias) 


         dflux = 0.0d0

         Pk = 0.0d0
         Pc = 0.0d0
         PP = 0.0d0
         du_dx = 0.0d0


         alpg = 1.0d0 - alps    


         do ii = N2-2, N3+2 
          do jj = M2-2, M3+2
           do kk = K2-2, K3+2


            thetas(ii,jj,kk) = max(thetas(ii,jj,kk),0.0d0)            


            dudx = (U(ii-2,jj,kk) - 8.0d0*U(ii-1,jj,kk) + & 
              8.0d0*U(ii+1,jj,kk) - U(ii+2,jj,kk))/(12.0d0*dx)   
                         
            dudy = (U(ii,jj-2,kk) - 8.0d0*U(ii,jj-1,kk) + & 
              8.0d0*U(ii,jj+1,kk) - U(ii,jj+2,kk))/(12.0d0*dy)   
 
            dudz = (U(ii,jj,kk-2) - 8.0d0*U(ii,jj,kk-1) + & 
              8.0d0*U(ii,jj,kk+1) - U(ii,jj,kk+2))/(12.0d0*dz)   

            dvdx = (V(ii-2,jj,kk) - 8.0d0*V(ii-1,jj,kk) + & 
              8.0d0*V(ii+1,jj,kk) - V(ii+2,jj,kk))/(12.0d0*dx)   
                         
            dvdy = (V(ii,jj-2,kk) - 8.0d0*V(ii,jj-1,kk) + & 
              8.0d0*V(ii,jj+1,kk) - V(ii,jj+2,kk))/(12.0d0*dy)   
 
            dvdz = (V(ii,jj,kk-2) - 8.0d0*V(ii,jj,kk-1) + & 
              8.0d0*V(ii,jj,kk+1) - V(ii,jj,kk+2))/(12.0d0*dz)   

            dwdx = (W(ii-2,jj,kk) - 8.0d0*W(ii-1,jj,kk) + & 
              8.0d0*W(ii+1,jj,kk) - W(ii+2,jj,kk))/(12.0d0*dx)   
                         
            dwdy = (W(ii,jj-2,kk) - 8.0d0*W(ii,jj-1,kk) + & 
              8.0d0*W(ii,jj+1,kk) - W(ii,jj+2,kk))/(12.0d0*dy)   
 
            dwdz = (W(ii,jj,kk-2) - 8.0d0*W(ii,jj,kk-1) + & 
              8.0d0*W(ii,jj,kk+1) - W(ii,jj,kk+2))/(12.0d0*dz)   



            Pk(ii,jj,kk,1,1) = alps(ii,jj,kk)*rhos*thetas(ii,jj,kk)
            Pk(ii,jj,kk,2,2) = Pk(ii,jj,kk,1,1)
            Pk(ii,jj,kk,3,3) = Pk(ii,jj,kk,1,1)



            call compute_g0(alps(ii,jj,kk),g0)

            musc = 4.0d0/5.0d0*alps(ii,jj,kk)*rhos*dias*g0 & 
                *(1.0d0+erest)*sqrt(thetas(ii,jj,kk)/3.14159265d0)
            xis = 4.0d0/3.0d0*alps(ii,jj,kk)*rhos*dias*g0 & 
                *(1.0d0+erest)*sqrt(thetas(ii,jj,kk)/3.14159265d0) 

            Pc1 = 0.0d0
            Pc2 = 0.0d0
            Pc3 = 0.0d0
 
            do l = 1, 3
             Pc1(l,l) = 2.0d0*(1.0d0+erest)*g0*rhos* & 
                 (alps(ii,jj,kk)**2.0d0)*thetas(ii,jj,kk)
            enddo

            du_dx(ii,jj,kk,1,1) = dudx
            du_dx(ii,jj,kk,1,2) = dudy
            du_dx(ii,jj,kk,1,3) = dudz
            du_dx(ii,jj,kk,2,1) = dvdx
            du_dx(ii,jj,kk,2,2) = dvdy
            du_dx(ii,jj,kk,2,3) = dvdz
            du_dx(ii,jj,kk,3,1) = dwdx
            du_dx(ii,jj,kk,3,2) = dwdy
            du_dx(ii,jj,kk,3,3) = dwdz

            do l = 1, 3 
            do m = 1, 3            
             du_dx_symm(l,m) = 0.5d0*(du_dx(ii,jj,kk,l,m)+du_dx(ii,jj,kk,m,l))
            enddo
            enddo

            Pc2 = -2.0d0*alps(ii,jj,kk)*musc*du_dx_symm

            Pc3(1,1) = -alps(ii,jj,kk)*xis*dudx  
            Pc3(2,2) = -alps(ii,jj,kk)*xis*dvdy  
            Pc3(3,3) = -alps(ii,jj,kk)*xis*dwdz  

            Pc(ii,jj,kk,:,:) = Pc1 + Pc2 + Pc3
 

            grad_theta(1) = (thetas(ii-2,jj,kk) - 8.0d0*thetas(ii-1,jj,kk) + & 
              8.0d0*thetas(ii+1,jj,kk) - thetas(ii+2,jj,kk))/(12.0d0*dx)   
            grad_theta(2) = (thetas(ii,jj-2,kk) - 8.0d0*thetas(ii,jj-1,kk) + & 
              8.0d0*thetas(ii,jj+1,kk) - thetas(ii,jj+2,kk))/(12.0d0*dy)   
            grad_theta(3) = (thetas(ii,jj,kk-2) - 8.0d0*thetas(ii,jj,kk-1) + & 
              8.0d0*thetas(ii,jj,kk+1) - thetas(ii,jj,kk+2))/(12.0d0*dz)   

            kappa_dilute = 75.0d0/384.0d0*sqrt(3.14159265d0)*rhos*dias*sqrt(thetas(ii,jj,kk))
 
            t1 = 2.0d0/g0/(1.0d0+erest)*(1.0d0+1.2d0*alps(ii,jj,kk) & 
                        *g0*(1.0d0+erest))**2.0d0
            t2 = 2.0d0*rhos*alps(ii,jj,kk)*alps(ii,jj,kk)*dias*(1.0d0+erest)*g0

            kappa_eff = t1*kappa_dilute + t2*sqrt(thetas(ii,jj,kk)/3.14159265d0)

            qk_plus_qc(ii,jj,kk,:) = -kappa_eff*grad_theta  

           enddo 
          enddo    
         enddo
         

!---------------------------------------

              PP = Pk + Pc

!-----------------------------------------
!  df/dx(i) = (f(i-2)-8f(i-1)+8f(i+1)-f(i+2))/(12*dx)
!-----------------------------------------


         do ii = N2, N3
          do jj = M2, M3
           do kk = K2, K3

            ! viscous flux, x-direction
            l = 1

            m = 1
            taux = -(PP(ii-2,jj,kk,l,m) - 8.0d0*PP(ii-1,jj,kk,l,m) + & 
               8.0d0*PP(ii+1,jj,kk,l,m) - PP(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = -(PP(ii,jj-2,kk,l,m) - 8.0d0*PP(ii,jj-1,kk,l,m) + & 
               8.0d0*PP(ii,jj+1,kk,l,m) - PP(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = -(PP(ii,jj,kk-2,l,m) - 8.0d0*PP(ii,jj,kk-1,l,m) + & 
               8.0d0*PP(ii,jj,kk+1,l,m) - PP(ii,jj,kk+2,l,m))/(12.0d0*dz)   

  
            dflux(ii,jj,kk,8) = dflux(ii,jj,kk,8) + taux + tauy + tauz

 

            ! viscous flux, y-direction
            l = 2

            m = 1
            taux = -(PP(ii-2,jj,kk,l,m) - 8.0d0*PP(ii-1,jj,kk,l,m) + & 
               8.0d0*PP(ii+1,jj,kk,l,m) - PP(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = -(PP(ii,jj-2,kk,l,m) - 8.0d0*PP(ii,jj-1,kk,l,m) + & 
               8.0d0*PP(ii,jj+1,kk,l,m) - PP(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = -(PP(ii,jj,kk-2,l,m) - 8.0d0*PP(ii,jj,kk-1,l,m) + & 
               8.0d0*PP(ii,jj,kk+1,l,m) - PP(ii,jj,kk+2,l,m))/(12.0d0*dz)   
 
 
            dflux(ii,jj,kk,9) = dflux(ii,jj,kk,9) + taux + tauy + tauz 


 
            ! viscous flux, z-direction
            l = 3
  
            m = 1
            taux = -(PP(ii-2,jj,kk,l,m) - 8.0d0*PP(ii-1,jj,kk,l,m) + & 
               8.0d0*PP(ii+1,jj,kk,l,m) - PP(ii+2,jj,kk,l,m))/(12.0d0*dx)   
            m = 2  
            tauy = -(PP(ii,jj-2,kk,l,m) - 8.0d0*PP(ii,jj-1,kk,l,m) + & 
               8.0d0*PP(ii,jj+1,kk,l,m) - PP(ii,jj+2,kk,l,m))/(12.0d0*dy)   
            m = 3
            tauz = -(PP(ii,jj,kk-2,l,m) - 8.0d0*PP(ii,jj,kk-1,l,m) + & 
               8.0d0*PP(ii,jj,kk+1,l,m) - PP(ii,jj,kk+2,l,m))/(12.0d0*dz)   
 

            dflux(ii,jj,kk,10) = dflux(ii,jj,kk,10) + taux + tauy + tauz 



            ! thermal conduction (solid)           

            m = 1
            qx = -(qk_plus_qc(ii-2,jj,kk,m) - 8.0d0*qk_plus_qc(ii-1,jj,kk,m) + & 
               8.0d0*qk_plus_qc(ii+1,jj,kk,m) - qk_plus_qc(ii+2,jj,kk,m))/(12.0d0*dx)   

            m = 2
            qy = -(qk_plus_qc(ii,jj-2,kk,m) - 8.0d0*qk_plus_qc(ii,jj-1,kk,m) + & 
               8.0d0*qk_plus_qc(ii,jj+1,kk,m) - qk_plus_qc(ii,jj+2,kk,m))/(12.0d0*dy)   

            m = 3
            qz = -(qk_plus_qc(ii,jj,kk-2,m) - 8.0d0*qk_plus_qc(ii,jj,kk-1,m) + & 
               8.0d0*qk_plus_qc(ii,jj,kk+1,m) - qk_plus_qc(ii,jj,kk+2,m))/(12.0d0*dz)   

  
            dflux(ii,jj,kk,11) = dflux(ii,jj,kk,11) + qx + qy + qz 




            ! viscous dissipation

            visc_diss = 0.0d0
            do l = 1, 3
            do m = 1, 3
              visc_diss = visc_diss + PP(ii,jj,kk,l,m)*du_dx(ii,jj,kk,l,m)  
            enddo     
            enddo     

            dflux(ii,jj,kk,11) = dflux(ii,jj,kk,11) + visc_diss



           enddo
          enddo
         enddo 


                   


         return

        end subroutine 

!-------------------------------------------------------------

