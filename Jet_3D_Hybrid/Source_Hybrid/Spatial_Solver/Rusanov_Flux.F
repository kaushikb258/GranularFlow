!-----------------------------------------------------------------------------------------------------
!       KAUSHIK CREATED THE FOLLOWING SUBROUTINES
!-----------------------------------------------------------------------------------------------------

       subroutine rusanov_flux(gamma,rhos,index,uu,flux) 
 
         implicit none

            real(kind=8) :: gamma, rhos 
            real(kind=8), dimension(1:2) :: rhog, ug, vg, wg, eg, pg 
            real(kind=8), dimension(1:2) :: alpg, alps 
            real(kind=8), dimension(1:2) :: us, vs, ws, thetas, es 
            real(kind=8) :: ykl, ykr      
            real(kind=8), dimension(1:2,1:13) :: uu
            real(kind=8) :: flux(1:12), cl, cr, lamtil, kel, ker
            real(kind=8) :: fl, fr, consl, consr   
            integer :: index  
            real(kind=8) :: nx, ny, nz, ql, qr 


            
            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif  


!--------------------------------------------

             ! 1: left
             ! 2: right


           rhog(1) = uu(1,1); rhog(2) = uu(2,1)

           ug(1) = uu(1,2); ug(2) = uu(2,2)
           vg(1) = uu(1,3); vg(2) = uu(2,3)
           wg(1) = uu(1,4); wg(2) = uu(2,4)

           eg(1) = uu(1,5); eg(2) = uu(2,5)
           pg(1) = uu(1,6); pg(2) = uu(2,6)
           ykl = uu(1,7); ykr = uu(2,7)


           alps(1) = uu(1,8); alps(2) = uu(2,8)
           alpg(1) = 1.0d0 - alps(1)
           alpg(2) = 1.0d0 - alps(2)

           us(1) = uu(1,9); us(2) = uu(2,9)
           vs(1) = uu(1,10); vs(2) = uu(2,10)
           ws(1) = uu(1,11); ws(2) = uu(2,11)  

           thetas(1) = uu(1,12); thetas(2) = uu(2,12)           
           es(1) = uu(1,13); es(2) = uu(2,13)           


!--------------------------------------------------


                     ! GAS PHASE


               ql = ug(1)*nx + vg(1)*ny + wg(1)*nz
               qr = ug(2)*nx + vg(2)*ny + wg(2)*nz


             cl = sqrt(gamma*pg(1)/rhog(1))
             cr = sqrt(gamma*pg(2)/rhog(2))



           lamtil = max(abs(ql)+cl,abs(qr)+cr)


           ! continuity
           fl = alpg(1)*rhog(1)*ql
           fr = alpg(2)*rhog(2)*qr
           consl = alpg(1)*rhog(1)
           consr = alpg(2)*rhog(2)
           flux(1) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! x-momentum
           fl = alpg(1)*rhog(1)*ql*ug(1) + alpg(1)*pg(1)*nx
           fr = alpg(2)*rhog(2)*qr*ug(2) + alpg(2)*pg(2)*nx
           consl = alpg(1)*rhog(1)*ug(1)
           consr = alpg(2)*rhog(2)*ug(2)
           flux(2) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! y-momentum
           fl = alpg(1)*rhog(1)*ql*vg(1) + alpg(1)*pg(1)*ny
           fr = alpg(2)*rhog(2)*qr*vg(2) + alpg(2)*pg(2)*ny
           consl = alpg(1)*rhog(1)*vg(1)
           consr = alpg(2)*rhog(2)*vg(2)
           flux(3) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)      

           ! z-momentum
           fl = alpg(1)*rhog(1)*ql*wg(1) + alpg(1)*pg(1)*nz
           fr = alpg(2)*rhog(2)*qr*wg(2) + alpg(2)*pg(2)*nz
           consl = alpg(1)*rhog(1)*wg(1)
           consr = alpg(2)*rhog(2)*wg(2)
           flux(4) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! energy 
           kel = 0.5d0*(ug(1)**2.0d0 + vg(1)**2.0d0 + wg(1)**2.0d0)
           ker = 0.5d0*(ug(2)**2.0d0 + vg(2)**2.0d0 + wg(2)**2.0d0)
           fl = alpg(1)*rhog(1)*ql*(eg(1) + kel) + alpg(1)*pg(1)*ql
           fr = alpg(2)*rhog(2)*qr*(eg(2) + ker) + alpg(2)*pg(2)*qr
           consl = alpg(1)*rhog(1)*(eg(1) + kel)
           consr = alpg(2)*rhog(2)*(eg(2) + ker)
           flux(5) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! mass fraction
           fl = alpg(1)*rhog(1)*ql*Ykl
           fr = alpg(2)*rhog(2)*qr*Ykr
           consl = alpg(1)*rhog(1)*ykl
           consr = alpg(2)*rhog(2)*ykr
           flux(6) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)


!--------------------------------------------------------------------

                     ! SOLID PHASE


               ql = us(1)*nx + vs(1)*ny + ws(1)*nz
               qr = us(2)*nx + vs(2)*ny + ws(2)*nz


             call compute_cs(alps(1),thetas(1),cl)
             call compute_cs(alps(2),thetas(2),cr)



           lamtil = max(abs(ql)+cl,abs(qr)+cr)


           ! continuity
           fl = alps(1)*rhos*ql
           fr = alps(2)*rhos*qr
           consl = alps(1)*rhos
           consr = alps(2)*rhos
           flux(7) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! x-momentum
           fl = alps(1)*rhos*ql*us(1) 
           fr = alps(2)*rhos*qr*us(2) 
           consl = alps(1)*rhos*us(1)
           consr = alps(2)*rhos*us(2)
           flux(8) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! y-momentum
           fl = alps(1)*rhos*ql*vs(1) 
           fr = alps(2)*rhos*qr*vs(2) 
           consl = alps(1)*rhos*vs(1)
           consr = alps(2)*rhos*vs(2)
           flux(9) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)      

           ! z-momentum
           fl = alps(1)*rhos*ql*ws(1) 
           fr = alps(2)*rhos*qr*ws(2) 
           consl = alps(1)*rhos*ws(1)
           consr = alps(2)*rhos*ws(2)
           flux(10) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! theta energy 
           fl = 3.0d0/2.0d0*alps(1)*rhos*ql*thetas(1) 
           fr = 3.0d0/2.0d0*alps(2)*rhos*qr*thetas(2) 
           consl = 3.0d0/2.0d0*alps(1)*rhos*thetas(1)
           consr = 3.0d0/2.0d0*alps(2)*rhos*thetas(2)
           flux(11) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)

           ! es energy
           fl = alps(1)*rhos*ql*es(1)
           fr = alps(2)*rhos*qr*es(2)
           consl = alps(1)*rhos*es(1)
           consr = alps(2)*rhos*es(2)
           flux(12) = 0.5d0*(fl+fr) - 0.5d0*lamtil*(consr - consl)


!--------------------------------------------------------------------

           return

          end subroutine

!---------------------------------------------------------------------------
