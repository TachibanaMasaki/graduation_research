module mod_flash_2p
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_flash_2p(GZ,Kwo0,V0,L0,W0)
    implicit none
    real(8),intent(in),dimension(n,grid)::GZ,Kwo0
    real(8),intent(inout),dimension(grid)::V0,L0,W0
    integer,dimension(n)::num
    real(8)::F,dF,dW0,epsilon=0.00000000000001d0
    integer::i,r,c,roopmax=100
    
    V0=0.0d0
    L0=0.5d0
    W0=0.5d0
    

    do i=1,grid
        do r=1,roopmax
            F=0.0d0
            dF=0.0d0
            do c=1,n
                F=F+GZ(c,i)*(1.0d0-Kwo0(c,i))/(1.0d0-W0(i)+W0(i)*Kwo0(c,i))
                dF=dF+GZ(c,i)*(1.0d0-Kwo0(c,i))**2.0d0/(1.0d0-W0(i)+W0(i)*Kwo0(c,i))**2.0d0
            end do
            
            dW0=-F/dF
            W0(i)=W0(i)+dW0
            if(dW0 < epsilon)then
                exit
            end if
        end do
        L0(i)=1.0d0-W0(i)
    end do
      
        
    end subroutine calc_flash_2p
end module mod_flash_2p