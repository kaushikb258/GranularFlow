!------------------------------------------------------------------
!                  KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!------------------------------------------------------------------


      subroutine source_terms(dt, gamma,mu,kappa,dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV, & 
          rhog,ug,vg,wg,pg,Tg,alps,us,vs,ws,thetas,es,source)

        implicit none

        real(kind=8) :: gamma, mu, kappa, dt
        integer :: N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: source
        REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: rhog,ug,vg,wg,pg,Tg & 
                                        ,alps,alpg, us,vs,ws,thetas,es
        real(kind=8) :: dx, dy, dz
        integer :: ii, jj, kk, l, m  
        real(kind=8) :: dudx, dvdy, dwdz, div_vel, vel, Ts
        real(kind=8) :: g0, fdrag(1:3)
        real(kind=8) :: alps_max, erest, rhos, dias
        real(kind=8) :: t1, t2, t3, Reyn, cd, qconv, Ksg, nuss, hsg 
        real(kind=8) :: C1, C2  
        real(kind=8) :: Q5, source5, ke, eg, time_scale, time_step, cg
        real(kind=8) :: cvsolid, Tlow, Thigh, Tmid, ee
        integer :: iter
        
        
        ! = 1 for moon; = 2 for mars
        integer, parameter :: moon_mars = 2 !1
        


        real(kind=8), parameter :: grav = 3.71d0 !1.62d0
        real(kind=8), parameter :: prandtl = 0.71d0

 



         call solid_constants(alps_max,erest,rhos,dias)


         source = 0.0d0

         alpg = 1.0d0 - alps


         do ii = N2-2, N3+2 
          do jj = M2-2, M3+2
           do kk = K2-2, K3+2

            source(ii,jj,kk,4) = alpg(ii,jj,kk)*rhog(ii,jj,kk)*grav
            source(ii,jj,kk,5) = alpg(ii,jj,kk)*rhog(ii,jj,kk)*wg(ii,jj,kk)*grav   

                       if(alps(ii,jj,kk).gt.1.0d-4) then

!-------------------------------------------------

            dudx = (us(ii-2,jj,kk) - 8.0d0*us(ii-1,jj,kk) + & 
              8.0d0*us(ii+1,jj,kk) - us(ii+2,jj,kk))/(12.0d0*dx)   
                         
            dvdy = (vs(ii,jj-2,kk) - 8.0d0*vs(ii,jj-1,kk) + & 
              8.0d0*vs(ii,jj+1,kk) - vs(ii,jj+2,kk))/(12.0d0*dy)   
 
            dwdz = (ws(ii,jj,kk-2) - 8.0d0*ws(ii,jj,kk-1) + & 
              8.0d0*ws(ii,jj,kk+1) - ws(ii,jj,kk+2))/(12.0d0*dz)   

            div_vel = dudx + dvdy + dwdz
            

!-------------------------------------------------

            t1 = ug(ii,jj,kk) - us(ii,jj,kk)
            t2 = vg(ii,jj,kk) - vs(ii,jj,kk)
            t3 = wg(ii,jj,kk) - ws(ii,jj,kk)

            vel = sqrt(t1*t1 + t2*t2 + t3*t3)


            Reyn = rhog(ii,jj,kk)*vel*dias/mu


!-------------------------------------------------

            ! Drag coefficient 


            if(Reyn.gt.5.0d0) then

              C1 = 24.0d0/Reyn + 4.4d0/sqrt(Reyn) + 0.42d0
              C2 = 4.0d0/3.0d0/alpg(ii,jj,kk)*  & 
               (1.75d0 + 150.0d0*alps(ii,jj,kk)/alpg(ii,jj,kk)/Reyn)
              
              if(alps(ii,jj,kk).le.0.08d0) then
               cd = C1
              else if(alps(ii,jj,kk).ge.0.45d0) then
               cd = C2
              else
               t1 = (alps(ii,jj,kk)-0.08d0)*C2 
               t2 = (0.45d0-alps(ii,jj,kk))*C1
               cd = (t1 + t2)/0.37d0
              endif  

            else
              cd = 0.0d0
            endif  

             ! To avoid unusually large values that can cause the code
             ! to crash


              cd = min(cd, 1.0d0) 



             Ksg = 0.75d0*cd*alps(ii,jj,kk)*rhog(ii,jj,kk)*vel/dias



            fdrag(1) = Ksg*(ug(ii,jj,kk)-us(ii,jj,kk))
            fdrag(2) = Ksg*(vg(ii,jj,kk)-vs(ii,jj,kk))
            fdrag(3) = Ksg*(wg(ii,jj,kk)-ws(ii,jj,kk))  

!-------------------------------------------------

            ! Nusselt number

            if(Reyn.le.200.0d0) then
             nuss = 2.0d0 + 0.106d0*Reyn*(prandtl**0.333d0) 
            else
             nuss = 2.274d0 + 0.6d0*(Reyn**0.67d0)*(prandtl**0.333d0)
            endif  


!            if(nuss.gt.500.0d0) then 
!             write(*,*) 'nusselt error ', nuss, Reyn, dias, alps(ii,jj,kk), vel
!            endif 


            nuss = min(nuss, 100.0d0)    


            hsg = 6.0d0*alps(ii,jj,kk)*kappa*nuss/dias/dias
       

!-------------------------------------------------

       if(moon_mars.eq.1) then

          ! Moon

            cvsolid = 766.0d0


            if(es(ii,jj,kk).le.1.2d0) then
              Ts = 10.0d0
              goto 832
            else if(es(ii,jj,kk).ge.0.1793d6) then
              cvsolid = 875.0d0
              Ts = 380.0d0 + (es(ii,jj,kk)-0.1793d6)/cvsolid 
              goto 832
            else
              ! do bisection

              iter = 0

              Tlow = 10.0d0
              Thigh = 380.0d0
 
831           continue
 
              Tmid = 0.5d0*(Tlow + Thigh)               

              ee = -1.649d-5*(Tmid**4.0d0)/4.0d0 & 
                   + 6.694d-3*(Tmid**3.0d0)/3.0d0 & 
                   + 2.233d0*Tmid*Tmid/2.0d0 - 48.41d0*Tmid

              if(ee.ge.es(ii,jj,kk)) then
                Thigh = Tmid
              else
                Tlow = Tmid 
              endif 

              if(abs(1.0d0 - Tlow/Thigh).le.1.0d-4) then
                Ts = Tmid
                goto 832
              endif 

              iter = iter + 1
 
              if(iter.le.25) then
               goto 831
              else
               ! not converging
                cvsolid = 766.0d0
                Ts = es(ii,jj,kk)/cvsolid 
                goto 832
              endif  


            endif



       else if(moon_mars.eq.2) then
           
         ! Mars   

         ee = (79.42d0)**2.0d0 + 4.0d0*1.15d0*es(ii,jj,kk)
         if(ee.ge.0.0d0) then
          Ts = (-79.42d0 + sqrt(ee))/(2.0d0*1.15d0)  
         else
          ! compute at 300 K
          cvsolid = 79.42d0 + 2.3d0*(300.0d0)
          Ts = es(ii,jj,kk)/cvsolid 
         endif
  
         goto 832

       else
         write(*,*) 'wrong moon_mars ', moon_mars
       endif





832         continue
!-------------------------------------------------

           ! one more check
           if(Ts.le.25.0d0.or.Ts.ge.10000.0d0) then 
            cvsolid = 766.0d0
            Ts = es(ii,jj,kk)/cvsolid
           endif
           Ts = max(Ts,25.0d0)
           Ts = min(Ts,5000.0d0)




            qconv = hsg*(Tg(ii,jj,kk) - Ts) 


!-------------------------------------------------

                ! Test for energy time scale

             ke = 0.5d0*(ug(ii,jj,kk)**2.0d0 + vg(ii,jj,kk)**2.0d0 & 
                       + wg(ii,jj,kk)**2.0d0)
 
             eg = pg(ii,jj,kk)/rhog(ii,jj,kk)/(gamma-1.0d0)

             Q5 = (eg + ke)*alpg(ii,jj,kk)*rhog(ii,jj,kk)

            t1 = fdrag(1)*us(ii,jj,kk) + fdrag(2)*vs(ii,jj,kk) & 
               + fdrag(3)*ws(ii,jj,kk)  

            source5 = abs(qconv + t1)


            time_scale = Q5/(source5 + 1.0d-14)
 

            cg = sqrt(gamma*pg(ii,jj,kk)/rhog(ii,jj,kk))
 
            t2 = sqrt(ug(ii,jj,kk)**2.0d0 + vg(ii,jj,kk)**2.0d0 & 
                    + wg(ii,jj,kk)**2.0d0)

            time_step = dt !0.25d0*dx/(t2 + cg)  


            if(time_step/time_scale.gt.0.25d0) then
              !write(*,*) 'stiff source terms ', time_step, time_scale, time_step/time_scale
              fdrag = 0.0d0
              Ksg = 0.0d0
              qconv = 0.0d0
            endif

    
!-------------------------------------------------

            ! momentum equations 
            source(ii,jj,kk,2) = -fdrag(1)
            source(ii,jj,kk,3) = -fdrag(2)
            source(ii,jj,kk,4) = source(ii,jj,kk,4) -fdrag(3) 
 
            source(ii,jj,kk,8) = fdrag(1)
            source(ii,jj,kk,9) = fdrag(2)
            source(ii,jj,kk,10) = fdrag(3) + alps(ii,jj,kk)*rhos*grav


!-------------------------------------------------

            t1 = fdrag(1)*us(ii,jj,kk) + fdrag(2)*vs(ii,jj,kk) & 
               + fdrag(3)*ws(ii,jj,kk)  
            source(ii,jj,kk,5) = source(ii,jj,kk,5) -qconv - t1    




            call compute_g0(alps(ii,jj,kk),g0)
            
            t1 = 3.0d0*(erest**2.0d0 - 1.0d0)*(alps(ii,jj,kk)**2.0d0)*rhos
            t2 = g0*thetas(ii,jj,kk)*(4.0d0/dias*sqrt(thetas(ii,jj,kk)/3.14159265d0) - div_vel)  
            
 

            source(ii,jj,kk,11) = t1*t2 - 3.0d0*Ksg*thetas(ii,jj,kk)
            source(ii,jj,kk,12) = qconv  

!-------------------------------------------------

                       endif

           enddo 
          enddo    
         enddo
         

         return

        end subroutine 

!-------------------------------------------------------------

