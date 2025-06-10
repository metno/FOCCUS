subroutine ippe(Nx,Ny,Nz,rho,zeta,h,rhoref,Cs_w,Cs_r,hc,s_w,s_rho,ppe)

implicit none

integer, intent(in)                               :: Nx, Ny, Nz
double precision, intent(in), dimension(Nz,Ny,Nx) :: rho
double precision, intent(in), dimension(Ny,Nx)    :: zeta, h
double precision, intent(in), dimension(Nz+1)     :: Cs_w, s_w
double precision, intent(in), dimension(Nz)       :: Cs_r, s_rho
double precision, intent(in)                      :: rhoref, hc
double precision, intent(out), dimension(Ny,Nx)   :: ppe

integer :: i,j,k
double precision :: S, F, dz(Nz), zw(Nz+1), zr(Nz), g

! Set constants
g = 9.81

! Run loop
do j=1,Ny
    do i=1,Nx
    
        ! Calculate z at w-points and dz values according to stretching and transform
        do k=1,Nz+1                
            S = (hc*s_w(k) + h(j,i)*Cs_w(k))/(hc + h(j,i))
            zw(k) = zeta(j,i) + (zeta(j,i) + h(j,i))*S
        end do

        ! Calculate z at rho-points according to stretching and transform
        do k=1,Nz                
            S = (hc*s_rho(k) + h(j,i)*Cs_r(k))/(hc + h(j,i))
            zr(k) = zeta(j,i) + (zeta(j,i) + h(j,i))*S
        end do
        
        do k=1,Nz                
            dz(k) = zw(k+1)-zw(k)
        end do
        
        ! Do the integration
        F = 0.0
        
        do k=1,Nz
            F = -(g/rhoref)*max(rhoref-rho(k,j,i),0.0)*zr(k)*dz(k) + F
        end do
        
        ppe(j,i) = F

    end do
end do

end subroutine ippe
