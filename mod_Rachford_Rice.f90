module mod_Rachford_Rice
    use mod_initial
    use mod_autodiff
    use mod_flash
    use mod_pvgauss
    implicit none
    contains 
    subroutine calc_Rachford_Rice(n0,Kgo,Kwo,V0,L0,W0)
        real(8),intent(in),dimension(n)::n0,Kgo,Kwo
        real(8),intent(inout)::V0,L0,W0
        integer::k,roopmax=50
        real(8),allocatable,dimension(:)::residual
        real(8)::dL0,dW0,error,epsilon=10.0d0**(-10.0d0)
        real(8),allocatable,dimension(:,:)::jacobian
        type(diffs),allocatable::fxs(:)

        do k=1,roopmax
            call calc_flash(n0,L0,W0,Kgo,Kwo,fxs)
            call jacobian_mat(fxs,jacobian)
            call outxs(fxs,residual)
            residual=-residual
            call pvgauss(m-1,jacobian,residual)

            dL0=residual(1)
            dW0=residual(2)
    
            L0=L0+dL0
            W0=W0+dW0
            V0=1.0d0-L0-W0

            error=dot_product(residual,residual)
            error=sqrt(error)

            if(error<epsilon)then         
                exit  
            end if
        end do
    end subroutine calc_Rachford_Rice
end module mod_Rachford_Rice