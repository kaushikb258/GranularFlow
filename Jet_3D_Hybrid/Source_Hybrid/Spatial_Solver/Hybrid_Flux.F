#if defined (VARIABLE_GAMMA)
#define _EULER_GAMMA	Gamma_Gas, Gamma 
#else
#define _EULER_GAMMA	Gamma
#endif

  Subroutine Hybrid_Flux (dt, Q, D_Flux, Flag, Flag_WENO)

  USE PseudoPack 
  ! KAUSHIK ADDED THIS USE STATEMENT


  implicit none

  REALTYPE            :: D_F_Epsilon
  REALTYPE :: dt

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV)           :: Q, D_Flux
  integer , dimension(N0:N5,M0:M5,K0:K5,3) , OPTIONAL :: Flag_WENO

  integer , dimension(MR_N0:MR_N5,MR_M0:MR_M5,MR_K0:MR_K5,3) , OPTIONAL :: Flag

  integer  :: Index

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: F, D_F
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: rhog, ug, vg, wg, eg, pg, Tg 
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: alpg, alps
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: us, vs, ws, thetas, es 

  REALTYPE, dimension(M0:M5,K0:K5,NV)       :: Lambda_x
  REALTYPE, dimension(N0:N5,K0:K5,NV)       :: Lambda_y
  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Lambda_z

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:3)    :: q_cond

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,1:NV)    :: xflux, yflux, zflux
  integer   :: ii, jj, kk
  
  REALTYPE, dimension(1:3) :: sw_x, sw_y, sw_z  
  REALTYPE, dimension(N0:N5) :: switchx, switchx1  
  REALTYPE, dimension(M0:M5) :: switchy, switchy1  
  REALTYPE, dimension(K0:K5) :: switchz, switchz1  

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Yk
  REALTYPE :: thetai, thetaip1, term1, term2

!-----------------------------------------------------------
  real(kind=8) :: alps_max, erest, rhos, dias

  real(kind=8), parameter :: Rgas = 400.0d0
 
!-----------------------------------------------------------

  integer, parameter :: sponge_on = 0
    
  REALTYPE, PARAMETER :: eps_sw =  1.0d-2 
  REALTYPE, dimension(1:6,1:13) :: uu 
  integer, parameter :: n_cells = 5 

  REALTYPE, dimension(1:NV)           :: dQdx, dQdy, dQdz
  

#if defined (VARIABLE_GAMMA)
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: Gamma

  call PS_Specific_Heat_Gamma (Gamma_Formulation, Q, _EULER_GAMMA)
#endif

#include "Adjust_D_Flux_Epsilon.i"


!-------------------------------------------------------------------------


                call solid_constants(alps_max,erest,rhos,dias)

 
!-------------------------------------------------------------------------

!  call PS_Euler_Primitive_Variable (       N0, N5, M0, M5, K0, K5, NV,  &
!                                    Gamma, Q, U, V, W, P, H, C)

  call cons_2_prim(N0,N5,M0,M5,K0,K5,NV,gamma,rhos, &
          Q,rhog,ug,vg,wg,eg,pg,alps,alpg,us,vs,ws,thetas,es)  

 
       Tg(:,:,:) = pg(:,:,:)/rhog(:,:,:)/Rgas

!--------------------------------------------------------------
     ! EXAMPLE: 0-24 interior cells, 6 ghost cells either side

             ! N0 = -6
             ! N1 = -1
             ! N2 =  0
             ! N3 =  24
             ! N4 =  25
             ! N5 =  30
!--------------------------------------------------------------
!--------------------------------------------------------------

           D_Flux = 0.0d0


           do kk = K0, K5
            do jj = M0, M5
             do ii = N0, N5

               if(eg(ii,jj,kk).le.0.0d0) then
                write(*,*) 'e < 0 in Hybrid_Flux ', eg(ii,jj,kk), ii, jj, kk
                !call ps_stop('Hybrid flux error1 ')
               endif


              Yk(ii,jj,kk) = Q(ii,jj,kk,6)/Q(ii,jj,kk,1)

                if(Yk(ii,jj,kk).lt.0.0d0.or.Yk(ii,jj,kk).gt.1.0d0) then
                  write(*,*) 'Yk out of bounds', Yk(ii,jj,kk), ii, jj, kk
                !call ps_stop('Hybrid flux error2 ')
                endif

             enddo
            enddo
           enddo
!--------------------------------------------------



           ! I-flux
           xflux = 0.0d0


           do kk = K2, K3
            do jj = M2, M3

             switchx = 0.0d0

             do ii = N2-1, N3+1
                
              term1 = alpg(ii+1,jj,kk) - 2.0d0*alpg(ii,jj,kk) + alpg(ii-1,jj,kk)
              term2 = alpg(ii+1,jj,kk) + 2.0d0*alpg(ii,jj,kk) + alpg(ii-1,jj,kk)
              sw_x(1) = abs(term1/term2) 

              term1 = pg(ii+1,jj,kk) - 2.0d0*pg(ii,jj,kk) + pg(ii-1,jj,kk)
              term2 = pg(ii+1,jj,kk) + 2.0d0*pg(ii,jj,kk) + pg(ii-1,jj,kk)
              sw_x(2) = abs(term1/term2) 

              term1 = rhog(ii+1,jj,kk) - 2.0d0*rhog(ii,jj,kk) + rhog(ii-1,jj,kk)
              term2 = rhog(ii+1,jj,kk) + 2.0d0*rhog(ii,jj,kk) + rhog(ii-1,jj,kk)
              sw_x(3) = abs(term1/term2) 


              if(maxval(sw_x(1:3)).ge.eps_sw) then
               switchx(ii+1) = 1.0d0 
              else
               switchx(ii+1) = 0.0d0
              endif

             enddo



             switchx1 = 0.0d0
             do ii = N2-1, N3+1
              if(switchx(ii+1).eq.1.0d0) then
               switchx1(ii+1-n_cells:ii+1+n_cells) = 1.0d0
              endif
             enddo



             do ii = N2-1, N3+1

              uu(1:6,1) = rhog(ii-2:ii+3,jj,kk)
              uu(1:6,2) = ug(ii-2:ii+3,jj,kk)
              uu(1:6,3) = vg(ii-2:ii+3,jj,kk)
              uu(1:6,4) = wg(ii-2:ii+3,jj,kk)
              uu(1:6,5) = eg(ii-2:ii+3,jj,kk)
              uu(1:6,6) = pg(ii-2:ii+3,jj,kk)
              uu(1:6,7) = Yk(ii-2:ii+3,jj,kk)
              uu(1:6,8) = alps(ii-2:ii+3,jj,kk)
              uu(1:6,9) = us(ii-2:ii+3,jj,kk)
              uu(1:6,10) = vs(ii-2:ii+3,jj,kk)
              uu(1:6,11) = ws(ii-2:ii+3,jj,kk)
              uu(1:6,12) = thetas(ii-2:ii+3,jj,kk)
              uu(1:6,13) = es(ii-2:ii+3,jj,kk)

           if(switchx1(ii+1).eq.1.0d0) then
               call WENO_cons_flux(gamma,rhos,1,uu,xflux(ii+1,jj,kk,1:NV))
           else
               ! DANGER
               call WENO_cons_flux(gamma,rhos,1,uu,xflux(ii+1,jj,kk,1:NV))
               !call central_flux(gamma,rhos,1,uu,xflux(ii+1,jj,kk,1:NV))
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = (xflux(ii,jj,kk,:) - xflux(ii+1,jj,kk,:))/dx
             enddo
            enddo    
           enddo


!--------------------------------------------------


           ! J-flux
           yflux = 0.0d0


           do kk = K2, K3
            do ii = N2, N3

             switchy = 0.0d0

             do jj = M2-1, M3+1

              term1 = alpg(ii,jj+1,kk) - 2.0d0*alpg(ii,jj,kk) + alpg(ii,jj-1,kk)
              term2 = alpg(ii,jj+1,kk) + 2.0d0*alpg(ii,jj,kk) + alpg(ii,jj-1,kk)
              sw_y(1) = abs(term1/term2) 

              term1 = pg(ii,jj+1,kk) - 2.0d0*pg(ii,jj,kk) + pg(ii,jj-1,kk)
              term2 = pg(ii,jj+1,kk) + 2.0d0*pg(ii,jj,kk) + pg(ii,jj-1,kk)
              sw_y(2) = abs(term1/term2) 

              term1 = rhog(ii,jj+1,kk) - 2.0d0*rhog(ii,jj,kk) + rhog(ii,jj-1,kk)
              term2 = rhog(ii,jj+1,kk) + 2.0d0*rhog(ii,jj,kk) + rhog(ii,jj-1,kk)
              sw_y(3) = abs(term1/term2) 



              if(maxval(sw_y(1:3)).ge.eps_sw) then
               switchy(jj+1) = 1.0d0 
              else
               switchy(jj+1) = 0.0d0
              endif

             enddo



             switchy1 = 0.0d0
             do jj = M2-1, M3+1
              if(switchy(jj+1).eq.1.0d0) then
               switchy1(jj+1-n_cells:jj+1+n_cells) = 1.0d0
              endif
             enddo



             do jj = M2-1, M3+1

              uu(1:6,1) = rhog(ii,jj-2:jj+3,kk)
              uu(1:6,2) = ug(ii,jj-2:jj+3,kk)
              uu(1:6,3) = vg(ii,jj-2:jj+3,kk)
              uu(1:6,4) = wg(ii,jj-2:jj+3,kk)
              uu(1:6,5) = eg(ii,jj-2:jj+3,kk)
              uu(1:6,6) = pg(ii,jj-2:jj+3,kk)
              uu(1:6,7) = Yk(ii,jj-2:jj+3,kk)
              uu(1:6,8) = alps(ii,jj-2:jj+3,kk)
              uu(1:6,9) = us(ii,jj-2:jj+3,kk)
              uu(1:6,10) = vs(ii,jj-2:jj+3,kk)
              uu(1:6,11) = ws(ii,jj-2:jj+3,kk)
              uu(1:6,12) = thetas(ii,jj-2:jj+3,kk)
              uu(1:6,13) = es(ii,jj-2:jj+3,kk)

           if(switchy1(jj+1).eq.1.0d0) then
               call WENO_cons_flux(gamma,rhos,2,uu,yflux(ii,jj+1,kk,1:NV)) 
           else
               ! DANGER
               call WENO_cons_flux(gamma,rhos,2,uu,yflux(ii,jj+1,kk,1:NV)) 
               !call central_flux(gamma,rhos,2,uu,yflux(ii,jj+1,kk,1:NV))
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = D_Flux(ii,jj,kk,:) & 
          + (yflux(ii,jj,kk,:) - yflux(ii,jj+1,kk,:))/dy
             enddo
            enddo    
           enddo

 
!--------------------------------------------------


           ! K-flux
           zflux = 0.0d0


           do ii = N2, N3
            do jj = M2, M3

             switchz = 0.0d0

             do kk = K2-1, K3+1

              term1 = alpg(ii,jj,kk+1) - 2.0d0*alpg(ii,jj,kk) + alpg(ii,jj,kk-1)
              term2 = alpg(ii,jj,kk+1) + 2.0d0*alpg(ii,jj,kk) + alpg(ii,jj,kk-1)
              sw_z(1) = abs(term1/term2) 

              term1 = pg(ii,jj,kk+1) - 2.0d0*pg(ii,jj,kk) + pg(ii,jj,kk-1)
              term2 = pg(ii,jj,kk+1) + 2.0d0*pg(ii,jj,kk) + pg(ii,jj,kk-1)
              sw_z(2) = abs(term1/term2) 

              term1 = rhog(ii,jj,kk+1) - 2.0d0*rhog(ii,jj,kk) + rhog(ii,jj,kk-1)
              term2 = rhog(ii,jj,kk+1) + 2.0d0*rhog(ii,jj,kk) + rhog(ii,jj,kk-1)
              sw_z(3) = abs(term1/term2) 


              if(maxval(sw_z(1:3)).ge.eps_sw) then
               switchz(kk+1) = 1.0d0 
              else
               switchz(kk+1) = 0.0d0
              endif

             enddo



             switchz1 = 0.0d0
             do kk = K2-1, K3+1
              if(switchz(kk+1).eq.1.0d0) then
               switchz1(kk+1-n_cells:kk+1+n_cells) = 1.0d0
              endif
             enddo



             do kk = K2-1, K3+1

              uu(1:6,1) = rhog(ii,jj,kk-2:kk+3)
              uu(1:6,2) = ug(ii,jj,kk-2:kk+3)
              uu(1:6,3) = vg(ii,jj,kk-2:kk+3)
              uu(1:6,4) = wg(ii,jj,kk-2:kk+3)
              uu(1:6,5) = eg(ii,jj,kk-2:kk+3)
              uu(1:6,6) = pg(ii,jj,kk-2:kk+3)
              uu(1:6,7) = Yk(ii,jj,kk-2:kk+3)
              uu(1:6,8) = alps(ii,jj,kk-2:kk+3)
              uu(1:6,9) = us(ii,jj,kk-2:kk+3)
              uu(1:6,10) = vs(ii,jj,kk-2:kk+3)
              uu(1:6,11) = ws(ii,jj,kk-2:kk+3)
              uu(1:6,12) = thetas(ii,jj,kk-2:kk+3)
              uu(1:6,13) = es(ii,jj,kk-2:kk+3)

           if(switchz1(kk+1).eq.1.0d0) then
               call WENO_cons_flux(gamma,rhos,3,uu,zflux(ii,jj,kk+1,1:NV)) 
           else
               ! DANGER
               call WENO_cons_flux(gamma,rhos,3,uu,zflux(ii,jj,kk+1,1:NV)) 
               !call central_flux(gamma,rhos,3,uu,zflux(ii,jj,kk+1,1:NV))
           endif

             enddo                     



            enddo    
           enddo


           do kk = K2, K3
            do jj = M2, M3
             do ii = N2, N3
              D_Flux(ii,jj,kk,:) = D_Flux(ii,jj,kk,:)  & 
          + (zflux(ii,jj,kk,:) - zflux(ii,jj,kk+1,:))/dz
             enddo
            enddo    
           enddo

!--------------------------------------------------
  
                 ! VISCOUS FLUXES (GAS)

      D_F = 0.0d0

       call viscous_flux_gas(gamma,mu_0,kappa_0,dx,dy,dz, & 
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,alpg,ug,vg,wg,pg,Tg,D_F)


    D_Flux = D_Flux + D_F

!--------------------------------------------------
  
                 ! VISCOUS FLUXES (SOLID)

      D_F = 0.0d0

       call viscous_flux_solid(dx,dy,dz, & 
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,alps,us,vs,ws,thetas,es,D_F)


    D_Flux = D_Flux + D_F


!----------------------------------------------------------

                 ! NOZZLING TERMS 

      D_F = 0.0d0

       call nozzling_terms(dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,alps,pg,us,vs,ws,D_F)


    D_Flux = D_Flux + D_F 

!----------------------------------------------------------

                  ! LES TERMS

      D_F = 0.0d0
 

        ! LES FOR GRANULAR
      call les_smagorinsky(gamma,dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV,rhog,ug,vg,wg,pg,Tg,alpg,D_F)


    D_Flux = D_Flux + D_F


!-------------------------------------------------------------------------

                    ! SOURCE TERMS

      D_F = 0.0d0

      call source_terms(dt, gamma,mu_0,kappa_0,dx,dy,dz, &
          N0,N5,N2,N3,M0,M5,M2,M3,K0,K5,K2,K3,NV, & 
          rhog,ug,vg,wg,pg,Tg,alps,us,vs,ws,thetas,es,D_F)


    D_Flux = D_Flux + D_F


!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

             if(sponge_on.eq.1) then



!      KAUSHIK CREATED THIS FOR SPONGE BC IN Z

         D_F = 0.0d0



      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(z(ii,jj,kk).ge.0.25d0) then

          ! use backward differences for gradients of Q 
          if(kk.eq.K0) then
           dQdz = 0.0d0
          else 
           dQdz(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj,kk-1,1:NV))/dz
          endif

          call sponge_zmax(gamma,z(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdz,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo



          D_Flux = D_Flux + D_F



!-------------------------------------------------------------------------
!      KAUSHIK CREATED THIS FOR SPONGE BC IN X

         D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(x(ii,jj,kk).le.-0.04d0) then

          if(ii.eq.N5) then
           dQdx(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii-1,jj,kk,1:NV))/dx
          else
           dQdx(1:NV) = (Q(ii+1,jj,kk,1:NV)-Q(ii,jj,kk,1:NV))/dx
          endif

          call sponge_xmin(gamma,x(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdx,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F

           
          D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(x(ii,jj,kk).ge.0.04d0) then

          if(ii.eq.N0) then
           dQdx(1:NV) = (Q(ii+1,jj,kk,1:NV)-Q(ii,jj,kk,1:NV))/dx
          else
           dQdx(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii-1,jj,kk,1:NV))/dx
          endif

          call sponge_xmax(gamma,x(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdx,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F 



!------------------------------------------------------------------------- 

!      KAUSHIK CREATED THIS FOR SPONGE BC IN Y

         D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(y(ii,jj,kk).le.-0.04d0) then

          if(jj.eq.M5) then
           dQdy(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj-1,kk,1:NV))/dy
          else
           dQdy(1:NV) = (Q(ii,jj+1,kk,1:NV)-Q(ii,jj,kk,1:NV))/dy
          endif

          call sponge_ymin(gamma,y(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdy,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F

           
          D_F = 0.0d0

      do ii = N0, N5
       do jj = M0, M5
        do kk = K0, K5

         if(y(ii,jj,kk).ge.0.04d0) then

          if(jj.eq.M0) then
           dQdy(1:NV) = (Q(ii,jj+1,kk,1:NV)-Q(ii,jj,kk,1:NV))/dy
          else
           dQdy(1:NV) = (Q(ii,jj,kk,1:NV)-Q(ii,jj-1,kk,1:NV))/dy
          endif

          call sponge_ymax(gamma,y(ii,jj,kk),Q(ii,jj,kk,1:NV),dQdy,D_F(ii,jj,kk,1:NV))
         endif

        enddo
       enddo
      enddo

          D_Flux = D_Flux + D_F 


             endif
       ! endif for sponge_on

!------------------------------------------------------------------------- 


  
            return

!--------------------------------------------------------------


  END Subroutine Hybrid_Flux

#if defined (_EULER_GAMMA)
#undef _EULER_GAMMA
#endif


!-----------------------------------------------------------------------------------------------------
