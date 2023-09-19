module mod_density_W_inj
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_density_W_inj(Po1,rho_inj)
        type(diffs),intent(in)::Po1
        type(diffs),intent(out)::rho_inj
        real(8),allocatable,dimension(:)::Cc_SC
        real(8)::W_SDS,W_NaCl,W_salt,T00,M_brine,M_salt,rho_inj_real

        type(diffs)::P00,rhow,rho_brine,vl_brine
        allocate(Cc_SC(t))

        Cc_SC(1)=Cwc_SC1
        Cc_SC(2)=Cwc_SC2
        Cc_SC(3)=Cwc_SC3
        Cc_SC(4)=Cwc_SC4
        Cc_SC(5)=Cwc_SC5
        Cc_SC(6)=Cwc_SC6
        Cc_SC(7)=Cwc_SC7
        Cc_SC(8)=Cwc_SC8
        Cc_SC(9)=Cwc_SC9
        Cc_SC(10)=Cwc_SC10
        Cc_SC(11)=Cwc_SC11

        W_SDS=Cc_SC(10)*MW10/(Cc_SC(9)*MW9+Cc_SC(10)*MW10+Cc_SC(11)*MW11)
        W_NaCl=Cc_SC(11)*MW11/(Cc_SC(9)*MW9+Cc_SC(10)*MW10+Cc_SC(11)*MW11)
        
        W_salt=W_SDS+W_NaCl

        P00=Po1*10.0d0**(-6.0d0)  !MPa
        T00=T0-273.15d0  !℃
        rhow=1.0d0+10.0d0**(-6.0d0)*(-80.0d0*T00-3.3d0*T00**2.0d0+0.00175d0*T00**3.0d0+489.0d0*P00-2.0d0*T00*P00  &
             +0.016d0*T00**2.0d0*P00-1.3d0*10.0d0**(-5.0d0)*T00**3.0d0*P00-0.333d0*P00**2.0d0-0.002d0*T00*P00**2.0d0)

        rho_brine=rhow+W_salt*(0.668d0+0.44d0*W_salt+10.0d0**(-6.0d0)*(300.0d0*P00-2400.0d0*P00*W_salt  &
                  +T00*(80.0d0+3.0d0*T00-3300.0d0*W_salt-13.0d0*P00+47.0d0*P00*W_salt)))


        if (Cc_SC(10) == 0 .AND. Cc_SC(11) == 0) then
            M_brine=Mw9*10.0d0**(-3.0d0)
        else
            M_salt=(Cc_SC(10)*MW10+Cc_SC(11)*MW11)/(Cc_SC(10)+Cc_SC(11))
            M_brine=Mw9*M_salt/((1.0d0-W_salt)*M_salt+W_salt*Mw9)
            M_brine=M_brine*10.0d0**(-3.0d0)
        end if

        vl_brine=M_brine/(rho_brine*10.0d0**3.0d0)

        rho_inj=1.0d0/vl_brine

        call out_diffsx(rho_inj,rho_inj_real)
        
    end subroutine calc_density_W_inj
end module