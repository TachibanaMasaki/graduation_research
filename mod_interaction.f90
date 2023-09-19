module mod_interaction
    use mod_initial
    implicit none
    contains
subroutine calc_interaction(Vc,Theta)
    implicit none
    real(8),intent(in),dimension(n)::Vc
    real(8),intent(out),dimension(n,n)::Theta
    integer::i,j

    do i=1,n
        do j=1,n
            Theta(i,j)=1.0d0-(2.0d0*sqrt(Vc(i)**(1.0d0/3.0d0)*Vc(j)**(1.0d0/3.0d0))  &
                       /(Vc(i)**(1.0d0/3.0d0)+Vc(j)**(1.0d0/3.0d0)))**(1.2d0)
        end do
    end do
   
end subroutine calc_interaction
end module