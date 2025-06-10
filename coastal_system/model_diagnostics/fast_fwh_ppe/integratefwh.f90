subroutine ifwh(Nx,Ny,Nz,salt,zeta,h,sref,Cs_w,hc,s_w,fwh)

implicit none

integer, intent(in)                               :: Nx, Ny, Nz
double precision, intent(in), dimension(Nz,Ny,Nx) :: salt
double precision, intent(in), dimension(Ny,Nx)    :: zeta, h
double precision, intent(in), dimension(Nz+1)     :: Cs_w, s_w
double precision, intent(in)                      :: sref, hc
double precision, intent(out), dimension(Ny,Nx)   :: fwh

integer :: i,j,k
double precision :: S, F, dz(Nz), z(Nz+1)

! Run loop
do j=1,Ny
    do i=1,Nx
    
        ! Calculate z/dz values according to stretching and transform
        do k=1,Nz+1                
            S = (hc*s_w(k) + h(j,i)*Cs_w(k))/(hc + h(j,i))
            z(k) = zeta(j,i) + (zeta(j,i) + h(j,i))*S
            if (z(k) + zeta(j,i) < -10.0) then 
               z(k) = -10.0 -zeta(j,i)
            end if
        end do
        
        do k=1,Nz                
            dz(k) = z(k+1)-z(k)
        end do
        
        ! Do the integration
        F = 0.0
        
        do k=1,Nz
            F = max(sref-salt(k,j,i),0.0)*(dz(k)/sref) + F
        end do
        
        fwh(j,i) = F

    end do
end do

end subroutine ifwh

