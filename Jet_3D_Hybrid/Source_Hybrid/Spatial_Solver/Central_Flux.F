                   subroutine central_flux(gamma,rhos,index,uu,flux)

          implicit none

          real(kind=8) :: gamma, rhos
          real(kind=8) :: uu(1:6,1:13), flux(1:12), ff(1:6,1:12) 
          integer :: i, index
          real(kind=8), dimension(1:6) :: rhog, ug, vg, wg, eg, pg, Yk
          real(kind=8), dimension(1:6) :: alpg, alps 
          real(kind=8), dimension(1:6) :: us, vs, ws, thetas, es

          real(kind=8), dimension(1:6) :: Et, qg, qs
          real(kind=8) :: ke
          real(kind=8) :: nx, ny, nz



            if(index.eq.1) then
              nx = 1.0d0; ny = 0.0d0; nz = 0.0d0
            else if(index.eq.2) then
              nx = 0.0d0; ny = 1.0d0; nz = 0.0d0
            else if(index.eq.3) then
              nx = 0.0d0; ny = 0.0d0; nz = 1.0d0
            else
              write(*,*) 'wrong index ', index
            endif



           rhog(1:6) = uu(1:6,1)
           ug(1:6) = uu(1:6,2)
           vg(1:6) = uu(1:6,3)             
           wg(1:6) = uu(1:6,4)
           eg(1:6) = uu(1:6,5)
           pg(1:6) = uu(1:6,6)  
           Yk(1:6) = uu(1:6,7)
                        
           alps(1:6) = uu(1:6,8)  
           do i = 1, 6
            alpg(i) = 1.0d0 - alps(i)
           enddo
           us(1:6) = uu(1:6,9)  
           vs(1:6) = uu(1:6,10)  
           ws(1:6) = uu(1:6,11)  
           thetas(1:6) = uu(1:6,12)  
           es(1:6) = uu(1:6,13)  
           
!-------------------------------------------------

           do i = 1, 6
            if(eg(i).le.0.0d0) then
              write(*,*) 'eg < 0 in central flux ', eg
              write(*,*) 'pg ', pg
              write(*,*) 'rhog ', rhog
              !call ps_stop('Central flux error ')
            endif
           enddo


!-------------------------------------------------


           do i = 1, 6
             qg(i) = ug(i)*nx + vg(i)*ny + wg(i)*nz
             qs(i) = us(i)*nx + vs(i)*ny + ws(i)*nz
             ke = 0.5d0*(ug(i)**2.0d0 + vg(i)**2.0d0 + wg(i)**2.0d0)
             Et(i) = alpg(i)*rhog(i)*(eg(i) + ke)
           enddo    


           do i = 1, 6
             ff(i,1) = alpg(i)*rhog(i)*qg(i)
             ff(i,2) = alpg(i)*rhog(i)*qg(i)*ug(i) + alpg(i)*pg(i)*nx
             ff(i,3) = alpg(i)*rhog(i)*qg(i)*vg(i) + alpg(i)*pg(i)*ny
             ff(i,4) = alpg(i)*rhog(i)*qg(i)*wg(i) + alpg(i)*pg(i)*nz
             ff(i,5) = qg(i)*(Et(i) + alpg(i)*pg(i))
             ff(i,6) = alpg(i)*rhog(i)*qg(i)*Yk(i)

             ff(i,7) = alps(i)*rhos*qs(i)
             ff(i,8) = alps(i)*rhos*qs(i)*us(i) 
             ff(i,9) = alps(i)*rhos*qs(i)*vs(i) 
             ff(i,10) = alps(i)*rhos*qs(i)*ws(i) 
             ff(i,11) = 3.0d0/2.0d0*alps(i)*rhos*qs(i)*thetas(i) 
             ff(i,12) = alps(i)*rhos*qs(i)*es(i)
           enddo


           ! i-2 = 1 
           ! i-1 = 2 
           ! i = 3 
           ! i+1 = 4 
           ! i+2 = 5 
           ! i+3 = 6 

            i = 3 


              flux(:) = (-ff(i-1,:) + 7.0d0*ff(i,:) + & 
                    7.0d0*ff(i+1,:) - ff(i+2,:))/12.0d0  


                    return
                    end subroutine
