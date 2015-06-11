                 program post

        implicit none

        real(kind=8), parameter :: gamma = 1.4d0 
        real(kind=8), parameter :: rhos = 2700.0d0 
        integer, parameter :: nprim = 12
        integer, parameter :: nghost = 6

        real(kind=8) :: time, dt
        integer :: step, nx, ny, nz, nv
        real(kind=8), allocatable, dimension(:,:,:,:) :: Q
        real(kind=8), allocatable, dimension(:,:,:) :: prim, prim1
        real(kind=8), allocatable, dimension(:,:) :: alps, alpg
        real(kind=8), allocatable, dimension(:,:) :: rhog, ug, vg, wg, pg 
        real(kind=8), allocatable, dimension(:,:) :: us, vs, ws, thetas, es 
        integer :: i, j, k, n, np, l, k1 
        real(kind=8) :: xx, ke, ie         
        integer :: file_num, imax, jmax, kmax  
        real(kind=8) :: dx, dy, dz
        character ( len = 100 ) filename
        real(kind=8) :: uu(1:5,1:5,1:5,1:3) 
        real(kind=8) :: term1, term2  
        real(kind=8) :: t1, t2, t3 

        !--------------------------------------
              ! Variables for paraview
              integer, parameter :: output_unit = 29
              character (len=100) title
              real(kind=8), allocatable, dimension(:,:,:) :: xyz
        !--------------------------------------


        ! cell size 
        dx = 0.075d0
        dy = dx
        dz = dx
   
 

           file_num = 51823
           !write(filename,'("jet_Restart_",I2.2,".dat")'),file_num
           !write(filename,'("jet_Restart_",I3.3,".dat")'),file_num
           !write(filename,'("jet_Restart_",I4.4,".dat")'),file_num
           write(filename,'("jet_Restart_",I5.5,".dat")'),file_num
           !write(filename,'("jet_Restart_",I6.6,".dat")'),file_num


        open(5,file=filename,form='unformatted')
          read(5) step, time, dt
          read(5) nx, ny, nz, nv, np
          allocate(Q(1:nx,1:ny,1:nz,1:nv))
          allocate(alpg(1:nx,1:nz))
          allocate(alps(1:nx,1:nz))
          allocate(rhog(1:nx,1:nz))
          allocate(ug(1:nx,1:nz))
          allocate(vg(1:nx,1:nz))
          allocate(wg(1:nx,1:nz))
          allocate(pg(1:nx,1:nz))
          allocate(us(1:nx,1:nz))
          allocate(vs(1:nx,1:nz))
          allocate(ws(1:nx,1:nz))
          allocate(thetas(1:nx,1:nz))
          allocate(es(1:nx,1:nz))



          do n = 1, nv
            read(5) Q(:,:,:,n)
          enddo
        close(5)


              print*, 'read file '

          
          j = ny/2
          do i = 1, nx
            do k = 1, nz

              alps(i,k) = Q(i,j,k,7)/rhos 
              alpg(i,k) = 1.0d0 - alps(i,k)

              rhog(i,k) = Q(i,j,k,1)/alpg(i,k)
              ug(i,k) = Q(i,j,k,2)/Q(i,j,k,1)
              vg(i,k) = Q(i,j,k,3)/Q(i,j,k,1)
              wg(i,k) = Q(i,j,k,4)/Q(i,j,k,1)
              ke = 0.5d0*(ug(i,k)**2.0d0 + vg(i,k)**2.0d0 + wg(i,k)**2.0d0)
              ie = Q(i,j,k,5)/Q(i,j,k,1) - ke
              pg(i,k) = (gamma-1.0d0)*rhog(i,k)*ie

              if(alps(i,k).gt.1.0d-6) then
               us(i,k) = Q(i,j,k,8)/Q(i,j,k,7)
               vs(i,k) = Q(i,j,k,9)/Q(i,j,k,7)
               ws(i,k) = Q(i,j,k,10)/Q(i,j,k,7)
 
               thetas(i,k) = Q(i,j,k,11)/Q(i,j,k,7)*2.0d0/3.0d0
               es(i,k) = Q(i,j,k,12)/Q(i,j,k,7)
              else
               us(i,k) = 0.0d0
               vs(i,k) = 0.0d0
               ws(i,k) = 0.0d0
               thetas(i,k) = 0.0d0
               es(i,k) = 0.0d0
              endif

            enddo 
          enddo  



           print*, 'line1 ', step, time, dt
           print*, 'line2 ', nx, ny, nz, nv, np


 
!------------------------------         


           imax = nx - 2*nghost 
           jmax = ny - 2*nghost 
           kmax = nz - 2*nghost 


           print*, 'imax, jmax, kmax ', imax, jmax, kmax



           title = 'jet_vtk'

!--------------------------------------------------------------

           allocate(xyz(1:imax,1:kmax,1:3))            
           allocate(prim(1:imax,1:kmax,1:nprim)) 
           allocate(prim1(1:imax,1:kmax,1:nprim)) 


             do i = 1, imax
              xyz(i,:,1) = dble(i-imax/2)*dx
             enddo
             do k = 1, kmax
              xyz(:,k,2) = dble(k-kmax/2)*dz
             enddo
             xyz(:,:,3) = 0.0d0 ! no y-direction



            do i = 1, imax
            do k = 1, kmax
             prim(i,k,1) = alpg(nghost+i,nghost+k)
             prim(i,k,2) = alps(nghost+i,nghost+k)
             prim(i,k,3) = log(rhog(nghost+i,nghost+k))
             prim(i,k,4) = ug(nghost+i,nghost+k)
             prim(i,k,5) = vg(nghost+i,nghost+k)
             prim(i,k,6) = wg(nghost+i,nghost+k)
             prim(i,k,7) = log(pg(nghost+i,nghost+k))
             prim(i,k,8) = us(nghost+i,nghost+k)
             prim(i,k,9) = vs(nghost+i,nghost+k)
             prim(i,k,10) = ws(nghost+i,nghost+k)
             prim(i,k,11) = thetas(nghost+i,nghost+k)
             prim(i,k,12) = es(nghost+i,nghost+k)
            enddo
            enddo

           write(filename,'("output_",I6.6,".vtk")'),file_num


            do i = 1, imax
            do k = 1, kmax
             k1 = kmax+1-k
             prim1(i,k,:) = prim(i,k1,:)
            enddo
            enddo
            prim1(:,:,6) = -prim1(:,:,6) 
            prim1(:,:,10) = -prim1(:,:,10) 



              call vtk_write(output_unit,filename,title,imax,kmax,xyz,nprim    &
              ,prim1(1:imax,1:kmax,1:nprim))


           deallocate(xyz,prim,prim1)  
!--------------------------------------------------------------

                 print*, 'ug ', sum(ug)
                 print*, 'pg ', sum(pg)
                 print*, 'us ', sum(us)
                 print*, 'thetas ', sum(thetas)

                 print*, 'max pg ', maxval(pg(:,:)) 
                 print*, 'min pg ', minval(pg(:,:)) 
                 print*, 'max rhog ', maxval(rhog(:,:)) 
                 print*, 'min rhog ', minval(rhog(:,:)) 

!--------------------------------------------------------------


               deallocate(pg,ug,vg,wg,rhog,Q,us,vs,ws,thetas,es,alpg,alps)


                    end program

!--------------------------------------------------------------

