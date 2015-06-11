  Subroutine Uniform_Flow(time,dt3)

  REALTYPE :: U_Gas_1, C_Gas_1, f
  real(kind=8) :: time, dt3 
  integer :: ii, jj, kk
  real(kind=8) :: alps_max, erest, rhos, dias 

  real(kind=8), parameter :: alps_init = 4.0d-1
  real(kind=8) :: alpg_init, egas
  real(kind=8), parameter :: cvsolid = 766.0d0
  real(kind=8), parameter :: thetas_init = 10.0d0

  real(kind=8) :: r, t1, t2, t3, t4, t5, t6, t7, t8, zint, phi
  real(kind=8), parameter :: LL = 4.0d0

 

    alpg_init = 1.0d0 - alps_init
    call solid_constants(alps_max,erest,rhos,dias)


  C_Gas_1 = SQRT(Gamma_Gas_1*P_Gas_1/Rho_Gas_1)
  U_Gas_1 = Mach*C_Gas_1

#if defined (ECHO_FLOW_CONDITION)
#include "Echo_Uniform_Flow_Condition.i"
#endif

 

  ! GAS-ONLY REGIONS
  Q(:,:,:,1) = Rho_Gas_1
  Q(:,:,:,2) = 0.0d0
  Q(:,:,:,3) = 0.0d0
  Q(:,:,:,4) = 0.0d0
  Q(:,:,:,5) = P_Gas_1/(Gamma_Gas_1-ONE) 
  Q(:,:,:,6) = Rho_Gas_1


  Q(:,:,:,7:12) = 0.0d0

 
  ! GRANULAR REGIONS
  do ii = N0, N5
  do jj = M0, M5
  do kk = K0, K5



            call random_number(r)
            phi = 0.8d0 + 0.4d0*r 
            t1 = 0.06d0*phi*sin(2.0d0*x(ii,jj,kk)/LL)*sin(2.0d0*y(ii,jj,kk)/LL)

            call random_number(r)
            phi = 0.8d0 + 0.4d0*r 
            t2 = 0.03d0*phi*sin(4.0d0*x(ii,jj,kk)/LL)*cos(4.0d0*y(ii,jj,kk)/LL)
 
            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t3 = 0.072d0*phi*cos(4.0d0*x(ii,jj,kk)/LL)*sin(8.0d0*y(ii,jj,kk)/LL)

            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t4 = 0.1d0*phi*sin(8.0d0*x(ii,jj,kk)/LL)*sin(4.0d0*y(ii,jj,kk)/LL)

            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t5 = 0.074d0*phi*cos(8.0d0*x(ii,jj,kk)/LL)*cos(8.0d0*y(ii,jj,kk)/LL)
       
            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t6 = 0.091d0*phi*cos(16.0d0*x(ii,jj,kk)/LL)*cos(8.0d0*y(ii,jj,kk)/LL)

            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t7 = 0.097d0*phi*sin(8.0d0*x(ii,jj,kk)/LL)*cos(16.0d0*y(ii,jj,kk)/LL)

            call random_number(r)
            phi = 0.8d0 + 0.4d0*r
            t8 = 0.075d0*phi*cos(16.0d0*x(ii,jj,kk)/LL)*cos(16.0d0*y(ii,jj,kk)/LL)

            zint = 10.0d0 + t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 


    if(z(ii,jj,kk).ge.zint) then
     Q(ii,jj,kk,1) = alpg_init*Rho_Gas_1
     Q(ii,jj,kk,2:4) = 0.0d0
     Q(ii,jj,kk,5) = alpg_init*P_Gas_1/(Gamma_Gas_1-ONE) 
     Q(ii,jj,kk,6) = Q(ii,jj,kk,1)
     Q(ii,jj,kk,7) = alps_init*rhos
     Q(ii,jj,kk,8:10) = 0.0d0
     Q(ii,jj,kk,11) = 3.0d0/2.0d0*alps_init*rhos*thetas_init
     Q(ii,jj,kk,12) = alps_init*rhos*cvsolid*300.0d0
    endif   

  enddo 
  enddo 
  enddo 


  call BC (Q,time,dt3)

  Specific_Heat = Gamma

  Q0_N0 = Q(N0:N2,:,:,:) ; Q0_N5 = Q(N3:N5,:,:,:)
  Q0_M0 = Q(:,M0:M2,:,:) ; Q0_M5 = Q(:,M3:M5,:,:)
  Q0_K0 = Q(:,:,K0:K2,:) ; Q0_K5 = Q(:,:,K3:K5,:)

  Q_0(:,0) = Q(N0,M0,K0,:) 
  Q_0(:,1) = Q(N5,M0,K0,:) 

  END Subroutine Uniform_Flow
