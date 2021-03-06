#include "pseudopack.h"


Module WENO_Euler_Short

  USE WENO_Coefficients
  USE WENO_Option

implicit NONE

PRIVATE
PUBLIC  :: PS_WENO_Euler_Short

CONTAINS

  Subroutine PS_WENO_Euler_Short (Order, N0, N5, i0, i1, F, F_Half,  &
                                         WENO_Epsilon)

  integer                    :: Order, N0, N5, i0, i1

  REALTYPE, dimension(N0:N5,-(Order-1)/2:(Order-1)/2,2) :: F
  REALTYPE, dimension(N0:N5                           ) :: F_Half

  REALTYPE, OPTIONAL         :: WENO_Epsilon
  REALTYPE                   :: Epsilon_WENO

  integer                    :: i, k, l, r

  REALTYPE, dimension(0:(Order-1)/2                            ) :: V, a, b
  REALTYPE, dimension(0:(Order-1)/2                            ) :: d
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2              ) :: C
  REALTYPE, dimension(0:(Order-1)/2,0:(Order-1)/2,0:(Order-1)/2) :: S

                             Epsilon_WENO = PS_WENO_Option%Epsilon
  if (PRESENT(WENO_Epsilon)) Epsilon_WENO = WENO_Epsilon

  call PS_WENO_Coefficient (Order, d, C, S)

  k = (Order-1)/2

  F_Half = ZERO

  do l = 1,2

  SELECT CASE (Order)
    CASE (1)
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-0,l)*SUM(S(0,0:k,r)*F(i,r-0:r,l)) 
        enddo

        a = d/(Epsilon_WENO+ b)**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/Sum(a)
      enddo
!$OMP END PARALLEL DO

    CASE (3)
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-1,l)*SUM(S(0,0:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(1,1:k,r)*F(i,r-0:r,l)) 
        enddo

        a = d/(Epsilon_WENO + b)**2 

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (5)
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-2,l)*SUM(S(0,0:k,r)*F(i,r-2:r,l)) + &
                 F(i,r-1,l)*SUM(S(1,1:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(2,2:k,r)*F(i,r-0:r,l)) 
        enddo

        a = d/(Epsilon_WENO + b)**2 

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (7) 
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-3,l)*SUM(S(0,0:k,r)*F(i,r-3:r,l)) + &
                 F(i,r-2,l)*SUM(S(1,1:k,r)*F(i,r-2:r,l)) + &
                 F(i,r-1,l)*SUM(S(2,2:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(3,3:k,r)*F(i,r-0:r,l))
        enddo

        a = d/(Epsilon_WENO + b)**2 

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (9) 
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-4,l)*SUM(S(0,0:k,r)*F(i,r-4:r,l)) + &
                 F(i,r-3,l)*SUM(S(1,1:k,r)*F(i,r-3:r,l)) + &
                 F(i,r-2,l)*SUM(S(2,2:k,r)*F(i,r-2:r,l)) + &
                 F(i,r-1,l)*SUM(S(3,3:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(4,4:k,r)*F(i,r-0:r,l))
        enddo

        a = d/(Epsilon_WENO + b)**2 

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

    CASE (11)
!$OMP PARALLEL DO
      do i = i0,i1  
        do r = 0,k
          V(r) = SUM(C(r,:)*F(i,r-k:r,l))

          b(r) = F(i,r-5,l)*SUM(S(0,0:k,r)*F(i,r-5:r,l)) + &
                 F(i,r-4,l)*SUM(S(1,1:k,r)*F(i,r-4:r,l)) + &
                 F(i,r-3,l)*SUM(S(2,2:k,r)*F(i,r-3:r,l)) + &
                 F(i,r-2,l)*SUM(S(3,3:k,r)*F(i,r-2:r,l)) + &
                 F(i,r-1,l)*SUM(S(4,4:k,r)*F(i,r-1:r,l)) + &
                 F(i,r-0,l)*SUM(S(5,5:k,r)*F(i,r-0:r,l))
        enddo

        a = d/(Epsilon_WENO + b)**2

#if defined (MAPPED_WENO)
        if (PS_WENO_Option%Weight_Type == 1) a = W(d,a)
#endif
        
        F_Half(i) = F_Half(i) + SUM(a*V)/SUM(a)
      enddo
!$OMP END PARALLEL DO

  END SELECT

  enddo

  END Subroutine PS_WENO_Euler_Short
!
! ========================================================
!
  FUNCTION W(d, W_k)
 
  REALTYPE                     :: S_W_k
  REALTYPE, dimension(:)       :: d, W_k
  REALTYPE, dimension(SIZE(d)) :: W, Omega

  S_W_k = SUM(W_k)
  Omega = W_k/S_W_k
  W     = Omega*(d*(d+ONE)+Omega*(Omega-THREE*d))/                 &
                (d*d      +Omega*(ONE  -  TWO*d))

  END FUNCTION 

END Module WENO_Euler_Short
