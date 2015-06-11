!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

               subroutine WENO_cons_flux(gamma,rhos,index,uu,flux) 
 
         implicit none

            integer, parameter :: nprim = 13

            real(kind=8) :: gamma, rhos
            real(kind=8), dimension(1:2) :: rhog,ug,vg,wg,eg,pg    
            real(kind=8), dimension(1:2) :: alpg, alps
            real(kind=8), dimension(1:2) :: us,vs,ws,thetas,es    
            real(kind=8) :: ykl, ykr      
            real(kind=8) :: uu(1:6,1:nprim), uiphm(1:nprim), uiphp(1:nprim)
            real(kind=8) :: uuu(1:2,1:nprim)
            real(kind=8) :: flux(1:12), cl, cr, kel, ker
            integer :: index
            real(kind=8) :: nx, ny, nz, ql, qr 


             ! Rusanov Riemann solver (very dissipative)
 
             
              ! uu is primitive variable


 
                   if(gamma.lt.1.32d0.or.gamma.gt.1.34d0) then
                     write(*,*) 'gamma error ', gamma
                   endif


!---------------------------------------------------------------

            call weno5_gas(gamma,nprim,uu,uiphm,uiphp)
            


            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif  


           rhog(1) = uiphm(1); rhog(2) = uiphp(1)
           ug(1) = uiphm(2); ug(2) = uiphp(2)
           vg(1) = uiphm(3); vg(2) = uiphp(3)
           wg(1) = uiphm(4); wg(2) = uiphp(4)
           pg(1) = uiphm(6); pg(2) = uiphp(6)
           ykl = uiphm(7); ykr = uiphp(7)


           ! interpolation may have violated EOS; ensure energy is conserved 
           eg(1) = pg(1)/rhog(1)/(gamma-1.0d0)
           eg(2) = pg(2)/rhog(2)/(gamma-1.0d0)
           uiphm(5) = eg(1)
           uiphp(5) = eg(2)


           alps(1) = uiphm(8); alps(2) = uiphp(8)
           alps(1) = min(max(alps(1),0.0d0),1.0d0)
           alps(2) = min(max(alps(2),0.0d0),1.0d0)
           alpg(1) = 1.0d0 - alps(1)
           alpg(2) = 1.0d0 - alps(2)

           us(1) = uiphm(9); us(2) = uiphp(9)
           vs(1) = uiphm(10); vs(2) = uiphp(10)
           ws(1) = uiphm(11); ws(2) = uiphp(11)
           thetas(1) = uiphm(12); thetas(2) = uiphp(12)
           es(1) = uiphm(13); es(2) = uiphp(13)
           

!-------------------------------------------------

              uuu(1,1:nprim) = uiphm
              uuu(2,1:nprim) = uiphp

!-------------------------------------------------


              if(rhog(1).le.0.0d0.or.rhog(2).le.0.0d0) then
               write(*,*) 'rho bug ', rhog
               write(*,*) 'uiphm ', uiphm
               write(*,*) 'uiphp ', uiphp
               write(*,*) 'uu1 ', uu(:,1)
               !call ps_stop('weno error 1 ')
              endif   
              if(pg(1).le.0.0d0.or.pg(2).le.0.0d0) then
               write(*,*) 'p bug ', pg, uu(:,6)
              ! use first order + Rusanov here
              uuu(1,1:nprim) = uu(3,:)
              uuu(2,1:nprim) = uu(4,:)
              endif   




!-------------------------------------------------

                   ! RUSANOV FLUX 

           call rusanov_flux(gamma,rhos,index,uuu,flux)


!-------------------------------------------------


           return

          end subroutine

!---------------------------------------------------------------------------

            subroutine weno5_gas(gamma,nprim,uu,uiphm,uiphp)

           implicit none

           integer :: nprim  
           real(kind=8) :: uu(1:6,1:nprim), uiphm(1:nprim), uiphp(1:nprim)
           real(kind=8) :: d0, d1, d2, dt2, dt1, dt0
           real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, & 
                               omega0, omega1, omega2
           real(kind=8) :: alphat0, alphat1, alphat2
           real(kind=8) :: sum1, p0, p1, p2, eps
           integer :: i, l

           real(kind=8) :: phi, phi1, phi2, ke, gamma
           integer :: imonotonic


              ! uu is primitive variable


!---------------------------------------------------              
                  ! first order

                    !uiphm = uu(3,:)
                    !uiphp = uu(4,:)
                     !return
                     !write(*,*) 'should not be here in weno'
                     !stop
!---------------------------------------------------              



                  d0 = 3.0d0/10.0d0
                  d1 = 6.0d0/10.0d0
                  d2 = 1.0d0/10.0d0
                  dt2 = 3.0d0/10.0d0
                  dt1 = 6.0d0/10.0d0
                  dt0 = 1.0d0/10.0d0

                  eps = 1.0d-12


                  uiphm = 0.0d0
                  uiphp = 0.0d0

!-----------------------

              ! compute uu(i+1/2)-

      do l = 1, nprim

                  i = 3
             if(uu(i,l).ne.uu(i-1,l)) then                 

           beta0 = 13.0d0/12.0d0*(uu(i,l) - 2.0d0*uu(i+1,l) + &
            uu(i+2,l))**2.0d0 + 1.0d0/4.0d0*(3.0d0*uu(i,l) - 4.0d0*uu(i+1,l) + uu(i+2,l))**2.0d0

           beta1 = 13.0d0/12.0d0*(uu(i-1,l) - 2.0d0*uu(i,l) + &
                uu(i+1,l))**2.0d0 + 1.0d0/4.0d0*(uu(i-1,l) - uu(i+1,l))**2.0d0

           beta2 = 13.0d0/12.0d0*(uu(i-2,l) - 2.0d0*uu(i-1,l) + uu(i,l))**2.0d0 + &
                 1.0d0/4.0d0*(uu(i-2,l) - 4.0d0*uu(i-1,l) + 3.0d0*uu(i,l))**2.0d0

           alpha0 = d0/(eps + beta0)**2.0d0
           alpha1 = d1/(eps + beta1)**2.0d0
           alpha2 = d2/(eps + beta2)**2.0d0

           sum1 = alpha0 + alpha1 + alpha2

           omega0 = alpha0/sum1
           omega1 = alpha1/sum1
           omega2 = alpha2/sum1


           p0 = 1.0d0/3.0d0*uu(i,l) + 5.0d0/6.0d0*uu(i+1,l) - 1.0d0/6.0d0*uu(i+2,l)
           p1 = -1.0d0/6.0d0*uu(i-1,l) + 5.0d0/6.0d0*uu(i,l) + 1.0d0/3.0d0*uu(i+1,l)
           p2 = 1.0d0/3.0d0*uu(i-2,l) - 7.0d0/6.0d0*uu(i-1,l) + 11.0d0/6.0d0*uu(i,l)


           uiphm(l) = omega0*p0 + omega1*p1 + omega2*p2


                      ! slope limiter

               phi1 = (uu(i+1,l)-uu(i,l))/(uu(i,l)-uu(i-1,l))
               phi2 = (uiphm(l)-uu(i,l))/(uu(i,l)-uu(i-1,l))       
               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                                       
               uiphm(l) = uu(i,l) + 0.5d0*(uu(i,l) - uu(i-1,l))*phi
             else
               uiphm(l) = uu(i,l)
             endif

      enddo

!-----------------------

                  ! compute uu(i+1/2)+


      do l = 1, nprim

                  i = 4
             if(uu(i,l).ne.uu(i+1,l)) then               

           beta0 = 13.0d0/12.0d0*(uu(i,l) - 2.0d0*uu(i+1,l) + &
            uu(i+2,l))**2.0d0 + 1.0d0/4.0d0*(3.0d0*uu(i,l) - 4.0d0*uu(i+1,l) + uu(i+2,l))**2.0d0

           beta1 = 13.0d0/12.0d0*(uu(i-1,l) - 2.0d0*uu(i,l) + &
                uu(i+1,l))**2.0d0 + 1.0d0/4.0d0*(uu(i-1,l) - uu(i+1,l))**2.0d0

           beta2 = 13.0d0/12.0d0*(uu(i-2,l) - 2.0d0*uu(i-1,l) + uu(i,l))**2.0d0 + &
                 1.0d0/4.0d0*(uu(i-2,l) - 4.0d0*uu(i-1,l) + 3.0d0*uu(i,l))**2.0d0

           alphat0 = dt0/(eps + beta0)**2.0d0
           alphat1 = dt1/(eps + beta1)**2.0d0
           alphat2 = dt2/(eps + beta2)**2.0d0

           sum1 = alphat0 + alphat1 + alphat2

           omega0 = alphat0/sum1
           omega1 = alphat1/sum1
           omega2 = alphat2/sum1

           p0 = 11.0d0/6.0d0*uu(i,l) - 7.0d0/6.0d0*uu(i+1,l) + 1.0d0/3.0d0*uu(i+2,l)
           p1 = 1.0d0/3.0d0*uu(i-1,l) + 5.0d0/6.0d0*uu(i,l) - 1.0d0/6.0d0*uu(i+1,l)
           p2 = -1.0d0/6.0d0*uu(i-2,l) + 5.0d0/6.0d0*uu(i-1,l) + 1.0d0/3.0d0*uu(i,l)


           uiphp(l) = omega0*p0 + omega1*p1 + omega2*p2




                      ! slope limiter
 
 
               phi1 = (uu(i-1,l)-uu(i,l))/(uu(i,l)-uu(i+1,l))
               phi2 = (uiphp(l)-uu(i,l))/(uu(i,l)-uu(i+1,l))        
               phi = max(0.0d0,min(2.0d0,2.0d0*phi1,2.0d0*phi2))
                        
               uiphp(l) = uu(i,l) - 0.5d0*(uu(i+1,l) - uu(i,l))*phi
             else
               uiphp(l) = uu(i,l)
             endif 
                      

      enddo

!-----------------------

              ! uu is primitive variable



           uiphm(5) = uiphm(6)/uiphm(1)/(gamma-1.0d0)
           uiphp(5) = uiphp(6)/uiphp(1)/(gamma-1.0d0)


              ! 1 = gas density 
              if(uiphm(1).le.0.0d0.or.uiphp(1).le.0.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif


              ! 5 = gas internal energy 
              if(uiphm(5).le.0.0d0.or.uiphp(5).le.0.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif

              ! 6 = gas pressure
              if(uiphm(6).le.0.0d0.or.uiphp(6).le.0.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif


              ! 8 = solid volume fraction
              if(uiphm(8).ge.1.0d0.or.uiphp(8).ge.1.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif
              if(uiphm(8).lt.0.0d0.or.uiphp(8).lt.0.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif



              ! 12 = solid granular temperature
              if(uiphm(12).lt.0.0d0.or.uiphp(12).lt.0.0d0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
              endif


!--------------------------------------------------

                ! check for monotonicity

             imonotonic = 1

             do l = 1, nprim
              if(uiphm(l).lt.minval(uu(:,l)).or.uiphm(l).gt.maxval(uu(:,l))) then
               imonotonic = 0
              endif
              if(uiphp(l).lt.minval(uu(:,l)).or.uiphp(l).gt.maxval(uu(:,l))) then
               imonotonic = 0
              endif
             enddo

             if(imonotonic.eq.0) then
                 ! something has gone awry; switch to first order here
                  ! first order
                    uiphm = uu(3,:)
                    uiphp = uu(4,:)
             endif


!--------------------------------------------------

            
            
           return
           end subroutine

!---------------------------------------------------------------------------
