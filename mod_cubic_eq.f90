module mod_cubic_eq
  contains
pure subroutine cubic_eq_other(a,b,c,phase,z)
    implicit none
!---   argument   ------------------------
    double precision, intent(out) :: z
    double precision, intent(in)  :: a, b, c
    integer, intent(in) :: phase
        
!---   original   ------------------------
    double precision, parameter :: pai = 4d0*atan(1d0)    
    double precision :: p,q,d,theta
    double precision, allocatable, dimension(:) :: x, y, LA

!---   Set Coefficients p, q, d   ---
        p = -(a**2)/9d0 + b/3d0
        q = 2d0*(a**3) /27d0 - a*b/3d0 + c
        d = q**2 + 4d0*(p**3)

!---   d>0 one real solution   ---
    if(d>0d0) then
        allocate(x(1),y(1),LA(2))
        d = sqrt(d)
        LA(1) = (-q + d)*0.5d0
        LA(2) = (-q - d)*0.5d0
        LA(:) = fLA( LA(:) )
            
        y(1) = sum( LA(:) )

!---   d=0 two real solutions   ---        
    elseif (d == 0d0) then
        allocate(x(2), y(2), LA(1))
        LA = -q*0.5d0
        LA = fLA( LA )

        y(1) = 2d0 * LA(1)
        y(2) = -LA(1)

!---   d<0 three real solutions   ---
    else
        allocate(x(3), y(3))
        theta = acos( -0.5d0 * q /  sqrt(-p**3)  )
        y(1) = 2d0 * sqrt(-p) * cos(  theta            / 3d0 )
        y(2) = 2d0 * sqrt(-p) * cos( (theta + 2d0*pai) / 3d0 )
        y(3) = 2d0 * sqrt(-p) * cos( (theta + 4d0*pai) / 3d0 )
    endif
        
    x(:) = y(:) -a/3d0
        
    if(phase >= 2) then
        z = minval(x(:))           !Liquid
    else
        z = maxval(x(:))           !Vapor
    endif
    
contains
    elemental pure double precision function fLA(LA)
        double precision, intent(in) :: LA
        
        if(LA>0d0) then
            fLA = LA**(1d0/3d0)
        else
            fLA = -( -LA**(1d0/3d0) )
        endif
    
    end function fLA
    
    
end subroutine cubic_eq_other
subroutine cubic_eq(A,B,C,initial,xk)
    implicit none
    real(8),intent(in)::A,B,C,initial
    real(8),intent(out)::xk
    real(8)::f,df,dxk,xkp1,error,epsi
    integer::k,iter
    
    iter=1000
    epsi=0.001d0

  !zfactorは理想気体のときである1を参考にしてxk=1とする。

    xk=initial

    do k=1,iter
      f=xk**3.0d0+A*xk**2.0d0+B*xk+C
      df=3.0d0*xk**2.0d0+2.0d0*A*xk+B
      dxk=-f/df
      xkp1=xk+dxk
      error=dabs(xkp1-xk)
      if(error<epsi)exit
      xk=xkp1
    end do

end subroutine cubic_eq
end module mod_cubic_eq
