module mod_Saturation
    use mod_initial
    use mod_autodiff
    implicit none
    contains
    subroutine calc_saturation(unknown,V,F,L,W,Nmol,Vp,fxs_sat)
        integer,intent(in)::unknown
        type(diffs),intent(in)::V,F,L,W
        type(diffs),intent(in),dimension(t)::Nmol
        type(diffs),intent(in),dimension(s)::Vp
        type(diffs),intent(out)::fxs_sat
        integer::com
        type(diffs)::ntotal

        call residualvectorset3(unknown*grid+BHin+BHout,ntotal)
        do com=1,t
            ntotal=ntotal+Nmol(com)
        end do

        fxs_sat=ntotal*V*Vp(1)+ntotal*F*Vp(2)+ntotal*L*Vp(3)+ntotal*W*Vp(4)-1.0d0
        
        
    end subroutine calc_Saturation
end module mod_Saturation