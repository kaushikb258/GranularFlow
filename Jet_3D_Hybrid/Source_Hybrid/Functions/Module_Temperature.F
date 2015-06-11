#include "pseudopack.h"


MODULE Temperature

  USE Gas_Properties
  USE Gas_Mixture_Formulation

implicit NONE

INTERFACE PS_Temperature
  MODULE PROCEDURE Temperature_0D
  MODULE PROCEDURE Temperature_1D
  MODULE PROCEDURE Temperature_2D
  MODULE PROCEDURE Temperature_3D
END INTERFACE

PRIVATE
PUBLIC :: PS_Temperature

CONTAINS
!
!===================================================================
!
  Subroutine Temperature_0D (NV, Q, P, T, EOS)

  integer                       :: NV 

  REALTYPE, dimension(NV)                   :: Q
  REALTYPE                                  :: P, T
  REALTYPE                                  :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state
  
              write(*,*) 'should not be in Temperature_0D '

             ! KAUSHIK CHANGED THIS
          T = P/Q(1)/287.0d0   
          return  
         

  END Subroutine Temperature_0D
!
!===================================================================
!
  Subroutine Temperature_1D (N0,N5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, NV 

  REALTYPE, dimension(N0:N5,NV)             :: Q
  REALTYPE, dimension(N0:N5)                :: P, T
  REALTYPE, dimension(N0:N5)                :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state


              write(*,*) 'should not be in Temperature_1D '

             ! KAUSHIK CHANGED THIS
          T = P/Q(:,1)/287.0d0
          return



  END Subroutine Temperature_1D
!
!===================================================================
!
  Subroutine Temperature_2D (N0,N5, M0,M5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, M0,M5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,NV)       :: Q
  REALTYPE, dimension(N0:N5,M0:M5)          :: P, T
  REALTYPE, dimension(N0:N5,M0:M5)          :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state


              write(*,*) 'should not be in Temperature_2D '
 
             ! KAUSHIK CHANGED THIS
          T = P/Q(:,:,1)/287.0d0
          return



  END Subroutine Temperature_2D
!
!===================================================================
!
  Subroutine Temperature_3D (N0,N5, M0,M5, K0,K5, NV, Q, P, T, EOS)

  integer                       :: N0,N5, M0,M5, K0,K5, NV 

  REALTYPE, dimension(N0:N5,M0:M5,K0:K5,NV) :: Q
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: P, T, alps, alpg, rhog
  REALTYPE, dimension(N0:N5,M0:M5,K0:K5)    :: f_1, f_2, MW

  integer                        , OPTIONAL :: EOS
  integer                                   :: Equation_Of_state
  integer :: ii, jj, kk
  real(kind=8) :: rhos


           rhos = 3100.0d0 

       do ii = N0, N5
       do jj = M0, M5
       do kk = K0, K5
         alps(ii,jj,kk) = Q(ii,jj,kk,7)/rhos
         alpg(ii,jj,kk) = 1.0d0 - alps(ii,jj,kk) 
         rhog(ii,jj,kk) = Q(ii,jj,kk,1)/alpg(ii,jj,kk)
       enddo 
       enddo 
       enddo 

              
             ! KAUSHIK CHANGED THIS
          T(:,:,:) = P(:,:,:)/rhog(:,:,:)/400.0d0
          return 


  END Subroutine Temperature_3D

END MODULE Temperature
