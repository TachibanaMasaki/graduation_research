module mod_Liquid_vapor_like
  use mod_initial
  use mod_autodiff
  use mod_gibbs
  use mod_pvgauss
  contains
  subroutine calc_Liquid_vapor_like(Po0,Xoil,Kgo0,like,PP)
    implicit none
    real(8),intent(in),dimension(grid)::Po0
    real(8),intent(inout),dimension(n,grid)::Xoil,Kgo0
    integer,intent(in),dimension(grid)::like
    integer,intent(out),dimension(grid)::PP
    real(8),allocatable,dimension(:,:)::jacobian
    real(8),allocatable,dimension(:)::zz,alpha,WW,residual,MF,a

    integer::r,i,j,c
    real(8)::Wt,lam,error,epsilon=0.000000000001d0
    type(diffs),allocatable::fxs(:)
    
    allocate(zz(n),alpha(n),WW(n),jacobian(n,n),residual(n),MF(n))
    !like=1 Å® vapor like(ÉKÉXÇ©ÇÁè≠Çµñ˚Ç™åªÇÍÇÈ)
    !like=2 Å® liquid like(ñ˚Ç©ÇÁè≠ÇµÉKÉXÇ™åªÇÍÇÈ)

    do i=1,grid  
        if(like(i)==1)then
        !vapor like(ÉKÉXÇ©ÇÁè≠Çµñ˚Ç™åªÇÍÇÈ) 
            do c=1,n
                zz(c)=Xoil(c,i)
                WW(c)=zz(c)/Kgo0(c,i)
                alpha(c)=2.0d0*sqrt(WW(c))
            end do
        else
        !liquid like(ñ˚Ç©ÇÁè≠ÇµÉKÉXÇ™åªÇÍÇÈ)
            do c=1,n
                zz(c)=Xoil(c,i)
                WW(c)=zz(c)*Kgo0(c,i)
                alpha(c)=2.0d0*sqrt(WW(c))
            end do
        end if

        do r=1,roopmax
            call calc_gibbs(alpha,Po0(i),zz,like(i),fxs)
            call jacobian_mat(fxs,jacobian)
            call outxs(fxs,residual)
            residual=-residual
            
            do c=1,n
                if(zz(c) <= 0.0d0)then
                    do j=1,n
                        jacobian(c,j)=0.0d0
                        jacobian(j,c)=0.0d0
                    end do
                    jacobian(c,c)=1.0d0
                    residual(c)=0.0d0
                end if
            end do
            
            call pvgauss(n,jacobian,residual)
            
            do c=1,n
                alpha(c)=alpha(c)+residual(c)
                WW(c)=(alpha(c)/2.0d0)**2.0d0
            end do

            error=maxval(residual)
            write(*,*)"error",error
            if(error < epsilon)then
                exit
            end if   
        end do
        
        Wt=0.0d0
        do c=1,n
            Wt=Wt+WW(c)
        end do
        
        do c=1,n
            MF(c)=WW(c)/Wt
        end do

        lam=1.0d0-log(Wt)
    
        if(lam>=1.0d0)then
            PP(i)=2
            write(*,*)2,"ëä"
        else
            PP(i)=3
            write(*,*)3,"ëä"
            if(like(i)==1)then 
                do c=1,n
                    Kgo0(c,i)=Xoil(c,i)/MF(c)
                end do
            else
                do c=1,n
                    Kgo0(c,i)=MF(c)/Xoil(c,i)
                end do
            end if
        end if
    end do   
    
    end subroutine calc_Liquid_vapor_like
end module mod_Liquid_vapor_like