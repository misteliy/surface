        PROGRAM cubecruncher
        IMPLICIT NONE
        INCLUDE "fftw3.f"
        INTEGER, PARAMETER :: dbl=8
        REAL ( dbl ), ALLOCATABLE  :: rho_xyz ( :, :, : ), r ( :, : ), fftw_avg (:,:,:)
        REAL ( dbl ), ALLOCATABLE :: surface ( :, :, :),fftw_out_sqr ( :, :)
        REAL ( dbl ), ALLOCATABLE :: xi(:),yi(:),zi(:),gridx(:),gridy(:)
        COMPLEX ( dbl ),ALLOCATABLE  :: fftw_out ( :, :)
        INTEGER :: natoms, nx, ny, nz, nf, ix, iy, iz, iatom, j,k,n
        INTEGER :: L1, L2, L3, U1, U2, U3, I1, I2, I3, iframe, nframe, iskip
        INTEGER :: argcount, i, iarg, ii, jj, kk, k1, k2, k3, is, geq0
        INTEGER :: iargc, xmax, ymax, zmax, jx, jy, jz, xgrid, ygrid, zgrid
        INTEGER :: plan(dbl), counter,nfiles,ifile,md,iwk   
        INTEGER,ALLOCATABLE :: iskiplist(:)
        REAL ( dbl)  :: kx,ky,sigma,wk 
        REAL ( dbl ) :: h ( 3, 3 ), dum_real, arg, dr ( 3 ), dr_au ( 3 ),  dv  , isoval
        REAL ( dbl ) :: pi, sigsq, gauss, rsq, rgrid ( 3 ), rdiff ( 3 ), dist,rmax, tmp
        CHARACTER ( len = 80 ) ::  wq_char
        CHARACTER ( len = 80 ) ::  fid_top, fid_bot, fid
        CHARACTER ( len = 80 ) :: inpfile, fxyz, forcefile
        CHARACTER ( len = 2 ),ALLOCATABLE :: asymb(:) 
        CHARACTER ( len = 80 ),ALLOCATABLE :: filelist(:)
!-------------------------------------------------------------------------------------
        counter = 0
        pi = ACOS ( -1.0_dbl )
        isoval = 0.016_dbl
        argcount = IARGC()
        DO i = 1, argcount
         CALL getarg ( i, wq_char )
         READ ( wq_char, * ) inpfile 
         inpfile  = TRIM ( inpfile )
         WRITE ( *, * ) 'INPFILE ', inpfile
        ENDDO
! READ nframs and xyz file name
        OPEN ( 11, file = inpfile )
        READ ( 11, * ) h(1,1), h(2,2), h(3,3)
        READ ( 11, * ) nx, ny, nz
        READ ( 11, * ) sigma
        READ ( 11, * ) nfiles
        allocate(filelist(nfiles))
        allocate(iskiplist(nfiles))
        READ ( 11, * ) filelist(:)
        READ ( 11, * ) iskiplist(:)
        CLOSE ( 11 )
        !------------------------------------------------------
        sigsq = sigma**2          
        rmax = -LOG(1.0e-6*SQRT((2*pi*sigsq)**3))*2.0_dbl*sigsq
        dr (1) = h (1,1)/nx
        dr (2) = h (2,2)/ny
        dr (3) = h (3,3)/nz
        dv = dr ( 1 ) * dr ( 2 ) * dr ( 3 )
        nf = nx/2 + 1
        ALLOCATE ( fftw_out ( nf, ny) )
        ALLOCATE ( fftw_avg (2, nf, ny) )
        ALLOCATE ( surface (2, nx, ny))
        ALLOCATE ( fftw_out_sqr ( nf, ny) )
        ALLOCATE ( rho_xyz ( nx, ny, nz ) )
        ALLOCATE ( xi(nx*ny),yi(nx*ny),zi(nx*ny) )
        fftw_avg = 0.0_dbl
        zmax = CEILING ( SQRT ( rmax ) / dr ( 3 ) )
        ymax = CEILING ( SQRT ( rmax ) / dr ( 2 ) )
        xmax = CEILING ( SQRT ( rmax ) / dr ( 1 ) )    
        if (xmax > nx/2) xmax = nx/2 
        if (ymax > ny/2) ymax = ny/2
        if (zmax > nz/2) zmax = nz/2 
        WRITE ( fid_top, '(a3,a4)') "top",".xyz"
        fid_top = TRIM ( fid_top )
        WRITE ( fid_bot, '(a3,a4)') "bot",".xyz"
        fid_bot = TRIM ( fid_bot )
        OPEN ( unit = 100, file=fid_bot )
        OPEN ( unit = 200, file=fid_top )
        call dfftw_plan_dft_r2c_2d(plan, nx, ny, surface(1,:,:), fftw_out, FFTW_ESTIMATE)
        WRITE(*,*) 'rmax',rmax,'sigma',sigma
        DO ifile=1,nfiles 
        fxyz = filelist(ifile)
        fxyz = TRIM ( fxyz )
        iskip = iskiplist(ifile)
        CALL READ_FRAMES(fxyz,nframe,natoms)
        IF (.NOT. ALLOCATED( r )) ALLOCATE ( r ( natoms, 3 ) )
        IF (.NOT. ALLOCATED(asymb )) ALLOCATE ( asymb(natoms) )
        WRITE(*,'(a4,a20,a10,i3,a5,i3)') 'file',fxyz,'filenumber',ifile,'from ',nfiles
        WRITE( *, * ) "number of frames in file ", nframe
        WRITE(*,*) 'number of atoms in simulation ', natoms
        WRITE(*,*) 'current iskip',iskip
        OPEN ( 12, file = fxyz )
!-----------------------------------------------------------------------------------
       DO iframe = 1,nframe
         rho_xyz ( :, :, : ) = 0.0_dbl
         surface ( :, :, : ) = 0._dbl
         fftw_out ( :, : ) = 0._dbl
         fftw_out_sqr ( :, :) = 0._dbl
         READ ( 12, * )
         READ ( 12, * )
         SELECT CASE ( MOD ( iframe, iskip ) )
         CASE DEFAULT
           DO iatom = 1, natoms
             READ ( 12, * ) asymb ( iatom ), r ( iatom, : )
           END DO
         CASE ( 0 )
           WRITE ( *, * ) 'IFRAME = ', iframe
! write the surface to xyz
           WRITE ( 100, * ) nx*ny
           WRITE ( 100, * )  iframe
           WRITE ( 200, * ) nx*ny
           WRITE ( 200, * ) iframe 
           DO iatom = 1, natoms
             READ ( 12, * ) asymb ( iatom ), r ( iatom, : )
!           WRITE ( *, * ) iframe, nframe,  asymb ( iatom ), r ( iatom, : )
! put them back in the box  0 < r < L
! real pbc wrap needed HERE !!!!!!!!!!!!!!!!!!!!!!!!! TO DO
             call pbc_wrap(r(iatom,:),h(1,1),h(2,2),h(3,3))
             xgrid = FLOOR ( r ( iatom, 1 )/dr ( 1 ) )
             ygrid = FLOOR ( r ( iatom, 2 )/dr ( 2 ) )
             zgrid = FLOOR ( r ( iatom, 3 )/dr ( 3 ) )
!             IF ( asymb ( iatom )(1:1) /= "O" ) CYCLE
             IF ( asymb ( iatom )(1:1) == "H" ) CYCLE
             DO ix = xgrid-xmax, xgrid+xmax
               DO iy = ygrid-ymax, ygrid+ymax
                 DO iz = zgrid-zmax, zgrid+zmax
                     rgrid = (/REAL(ix,dbl)*dr(1),REAL(iy,dbl)*dr(2),REAL(iz,dbl)*dr(3)/)
                     rdiff(:) = rgrid(:)-r(iatom,:)
                     rsq = DOT_PRODUCT ( rdiff, rdiff )
                     gauss = exp (-rsq/(2.0_dbl*sigsq))
!define local indices 
                     ii = ix
                     jj = iy
                     kk = iz
!PBC local indices
                     IF ( ii > nx ) ii = ii - nx
                     IF ( ii < 1 ) ii = ii + nx
                     IF ( jj > ny ) jj = jj - ny
                     IF ( jj < 1 ) jj = jj + ny
                     IF ( kk > nz ) kk = kk - nz
                     IF ( kk < 1 ) kk = kk + nz
                     rho_xyz ( ii, jj, kk ) = rho_xyz ( ii, jj, kk ) + &
                                       gauss/SQRT((2.0_dbl*pi*sigsq)**3)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
           print *,rho_xyz(nx/2,ny/2,nz/2)
!--------------------------- this was rho calc - now surface and fftw-----------------------
           DO ix = 1, nx
             DO iy = 1, ny
               DO iz = 1, nz
                  IF ( rho_xyz ( ix, iy, iz ) < isoval ) THEN 
                    CYCLE
                  ELSE
                    WRITE ( 100, * ) 'X', ix * dr ( 1 ), iy * dr ( 2 ), iz * dr ( 3 )
                    surface (1,ix,iy) = iz * dr( 3 )
                    EXIT
                  END IF
               END DO
               DO iz = nz, 1, -1
                  IF ( rho_xyz ( ix, iy, iz ) < isoval ) THEN 
                    CYCLE
                  ELSE 
                    WRITE ( 200, * ) 'X', ix * dr ( 1 ), iy * dr ( 2 ), iz * dr ( 3 )
                     surface (2,ix,iy) = iz * dr( 3 )
                    EXIT
                  END IF
               END DO
             END DO 
           END DO
!-------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------           
           do i=1,2
              call dfftw_execute_dft_r2c(plan, surface(i,:,:),fftw_out)
              fftw_out_sqr = (dr(1)*dr(2)/(nx*ny))*(fftw_out*CONJG(fftw_out))
              fftw_avg(i,:,:) = fftw_avg(i,:,:) + fftw_out_sqr
              counter = counter + 1
           end do
         END SELECT
       ENDDO
!--------------------------------------------------------------------------------------
        close(12)
        ENDDO
       call dfftw_destroy_plan(plan)
!-------------------------------------------------------------------------------------
       fftw_avg = fftw_avg/real(counter/2,dbl)
       do ix = 1, nf
         k1 = ix - 1
         do iy = 1, ny
              k2 = iy - 1
              IF ( k2 > ny/2 ) k2 = k2 - ny
              kx = pi/(real(nf)*dr(1))*k1 
              ky = 2*pi/(real(ny)*dr(2))*k2
              if (ix .eq. 1 .and. iy .eq. 1) then
                  tmp = 0.0
              else
                 tmp = (fftw_avg(1,ix,iy)+fftw_avg(2,ix,iy))/2.0d0 
              endif
              write( 3 , * ) sqrt(kx**2 + ky**2), fftw_avg(1,ix,iy)
              write( 2 , * ) sqrt(kx**2 + ky**2), fftw_avg(2,ix,iy)
              write( 1 , * ) sqrt(kx**2 + ky**2), &
              tmp,kx,ky
         enddo
       enddo
       CLOSE ( 100 )
       CLOSE ( 200 )
! Average block
       DEALLOCATE ( rho_xyz )
       DEALLOCATE ( fftw_avg )
       DEALLOCATE ( fftw_out )
       DEALLOCATE ( r )
       DEALLOCATE ( asymb )
       DEALLOCATE ( surface )
       END PROGRAM cubecruncher

!---------------------------------------------------------------------
      subroutine pbc_wrap(r,x,y,z)
      implicit none
      integer, parameter :: dbl=8
      real(dbl),dimension(3),intent(in out) :: r
      real(dbl),intent(in) :: x,y,z
      do while  ( r ( 1 ) > x ) 
           r ( 1 ) = r ( 1 ) - x 
      enddo
      do while  ( r ( 1 ) < 0._dbl )
           r ( 1 ) = r ( 1 ) + x 
      enddo
      do while  ( r ( 2 ) > y ) 
           r ( 2 ) = r ( 2 ) - y  
      enddo
      do while  ( r ( 2 ) < 0._dbl ) 
           r ( 2 ) = r ( 2 ) + y 
      enddo 
      do while  ( r ( 3 ) > z ) 
           r ( 3 ) = r ( 3 ) - z 
      enddo
      do while  ( r ( 3 ) < 0._dbl ) 
           r ( 3 ) = r ( 3 ) + z 
      enddo
      end subroutine pbc_wrap 


      subroutine read_frames(xyzfile,nframes,natom)
      implicit none
      character(len=*) :: xyzfile
      integer, intent(out) :: nframes,natom
      integer :: ioerror,i,j
      integer, parameter :: nframes_max = 100000000
      character(len=4) :: dummyc
      xyzfile = trim(xyzfile) 
      open(10,file=xyzfile,iostat=ioerror)
      if ( ioerror /= 0 ) then
      nframes = 0
      print *, 'error opening file'
      return
      end if
      do i = 1, nframes_max
       read(10,*,end=100) natom
       read(10,*) 
       do j = 1,natom
          read(10,*)
       end do 
      end do
    
100   nframes=i-1
      close(10)
    
      return
      end subroutine read_frames
