module mod_viscosity
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_viscosity(unknown,Ptotal,Xpc,Vp,FM1,Sw,So,mu)
        implicit none
        integer,intent(in)::unknown
        type(diffs),intent(in)::Ptotal
        type(diffs),intent(in),dimension(s,t)::Xpc
        type(diffs),intent(in),dimension(s)::Vp
        type(diffs),intent(in)::FM1,Sw,So
        type(diffs),intent(out),dimension(s)::mu
        real(8),allocatable,dimension(:)::Pc,Tc,Vc,MW,Tr,muc,rhoc,omega
        real(8)::av0,av1,av2,av3,av4,betaw,betae,Cs,TT,muwSC
        integer::c,com,phase
        type(diffs)::dummy,xi_V,xi_L,rho_V,rho_L,mu_V,mu_L,C_NaCl,betav,betar,Av,Bv,muw  &
                    ,FM2,FM7,FM
        allocate(Pc(n),Tc(n),Vc(n),MW(t),Tr(n),muc(n),rhoc(n),omega(n))
!!気相、油相の相粘度 (cP)  (cP)→(Pa*s)にした
!1:CO2, 2:BENZEN, 3:CH4, 4:H2O, 5:surfactant, 6:NaCl
        !atm
        Pc(1)=Pc1/101325.0d0  
        Pc(2)=Pc2/101325.0d0  
        Pc(3)=Pc3/101325.0d0  
        Pc(4)=Pc4/101325.0d0
        Pc(5)=Pc5/101325.0d0  
        Pc(6)=Pc6/101325.0d0  
        Pc(7)=Pc7/101325.0d0  
        Pc(8)=Pc8/101325.0d0
        Pc(9)=Pc9/101325.0d0

        !K
        Tc(1)=Tc1  
        Tc(2)=Tc2  
        Tc(3)=Tc3
        Tc(4)=Tc4 
        Tc(5)=Tc5 
        Tc(6)=Tc6  
        Tc(7)=Tc7
        Tc(8)=Tc8
        Tc(9)=Tc9

        !m^3/mol
        Vc(1)=Vc1
        Vc(2)=Vc2 
        Vc(3)=Vc3
        Vc(4)=Vc4
        Vc(5)=Vc5
        Vc(6)=Vc6 
        Vc(7)=Vc7
        Vc(8)=Vc8
        Vc(9)=Vc9

        !g/mol
        MW(1)=MW1
        MW(2)=MW2
        MW(3)=MW3
        MW(4)=MW4
        MW(5)=MW5
        MW(6)=MW6
        MW(7)=MW7
        MW(8)=MW8
        MW(9)=MW9
        MW(10)=MW10
        MW(11)=MW11

        !c成分の偏心因子 (-)
        omega(1)=omega1
        omega(2)=omega2
        omega(3)=omega3
        omega(4)=omega4
        omega(5)=omega5
        omega(6)=omega6
        omega(7)=omega7
        omega(8)=omega8
        omega(9)=omega9

        !c成分の還元温度 (-)
        do c=1,n
            Tr(c)=T0/Tc(c)   
        end do

    !gas相におけるξ
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*Tc(c)
        end do
        dummy=dummy**(1.0d0/6.0d0)
        xi_V=dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*MW(c)
        end do
        dummy=dummy**(-1.0d0/2.0d0)
        xi_V=xi_V*dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*Pc(c)
        end do
        dummy=dummy**(-2.0d0/3.0d0)
        xi_V=xi_V*dummy

    !oil相におけるξ
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*Tc(c)
        end do
        dummy=dummy**(1.0d0/6.0d0)
        xi_L=dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*MW(c)
        end do
        dummy=dummy**(-1.0d0/2.0d0)
        xi_L=xi_L*dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*Pc(c)
        end do
        dummy=dummy**(-2.0d0/3.0d0)
        xi_L=xi_L*dummy


    !gas相におけるρ
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*Vc(c)
        end do
        rho_V=dummy/Vp(1)

    !oil相におけるρ
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*Vc(c)
        end do
        rho_L=dummy/Vp(3)

    !相粘度 
    !gas相の粘度
        !成分粘度(気体)
        do c=1,n
            muc(c)=MW(c)**(1.0d0/2.0d0)*Pc(c)**(2.0d0/3.0d0)/Tc(c)**(1.0d0/6.0d0)  &
                   *(4.610d0*Tr(c)**0.618d0-2.04d0*exp(-0.449d0*Tr(c))+1.94d0*exp(-4.058d0*Tr(c))+0.1d0)*10**(-4.0d0)
        end do
        
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*muc(c)*MW(c)**0.5d0
        end do
        mu_V=dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(1,c)*MW(c)**0.5d0
        end do
        mu_V=mu_V/dummy

    !oil相の粘度
        !成分粘度(液体)
        do c=1,n
            muc(c)=MW(c)**(1.0d0/2.0d0)*Pc(c)**(2.0d0/3.0d0)/Tc(c)**(1.0d0/6.0d0)  &
                   *(4.610d0*Tr(c)**0.618d0-2.04d0*exp(-0.449d0*Tr(c))+1.94d0*exp(-4.058d0*Tr(c))+0.1d0)*10**(-4.0d0)
        end do
      
        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*muc(c)*MW(c)**0.5d0
        end do
        mu_L=dummy

        call residualvectorset3(unknown*grid+BHin+BHout,dummy)
        do c=1,n
            dummy=dummy+Xpc(3,c)*MW(c)**0.5d0
        end do
        mu_L=mu_L/dummy

    !Jossi, Stiel and Thodosの粘度定数
        av0=0.10230d0
        av1=0.023364d0
        av2=0.058533d0
        av3=-0.040758d0
        av4=0.0093324d0

        mu(1)=mu_V+((av0+av1*rho_V+av2*rho_V**2.0d0+av3*rho_V**3.0d0+av4*rho_V**4.0d0)**4.0d0-10**(-4.0d0))/xi_V
        mu(3)=mu_L+((av0+av1*rho_L+av2*rho_L**2.0d0+av3*rho_L**3.0d0+av4*rho_L**4.0d0)**4.0d0-10**(-4.0d0))/xi_L

    !1~7までがoil, 8:CO2, 9:H2O, 10:Surfactant, 11:NaCl
!!水相の相粘度 (cP)
        call residualvectorset3(unknown*grid+BHin+BHout,C_NaCl)
        do com=1,t
            C_NaCl=C_NaCl+Xpc(4,com)*MW(com)
        end do
        C_NaCl=Xpc(4,11)*1000.0d0/C_NaCl

        !水とNaClだけ
       !C_NaCl=Xpc(4,11)*1000.0d0/(Xpc(4,9)/(Xpc(4,9)+Xpc(4,11))*MW(9)+Xpc(4,11)/(Xpc(4,9)+Xpc(4,11))*MW(11))

        TT=T0-273.15d0
        betaw=-1.297d0+0.0574d0*TT-0.697d0*10.0d0**(-3.0d0)*TT**(2.0d0)  &
              +0.447d0*10.0d0**(-5.0d0)*TT**3.0d0-0.105d0*10.0d0**(-7.0d0)*TT**(4.0d0)

        betae=0.545d0+0.28*10.d0**(-2.0d0)*TT-betaw

        Cs=6.044d0+0.28*10.0d0**(-2.0d0)*TT+0.36*10.0d0**(-4.0d0)*TT**(2.0d0)
        betar=2.5d0*C_NaCl/Cs-2.0d0*(C_NaCl/Cs)**2.0d0+0.5d0*(C_NaCl/Cs)**3.0d0

        betav=betae*betar+betaw

        Av=0.3324d0*10.0d0**(-1.0d0)*C_NaCl+0.3624d0*10.0d0**(-2.0d0)*C_NaCl**2.0d0-0.1879d0*10.0d0**(-3.0d0)*C_NaCl**3.0d0
        Bv=0.102d0*10.0d0**(-1.0d0)*C_NaCl**2.0d0-0.702d0*10.0d0**(-3.0d0)*C_NaCl**3.0d0-0.396d0*10.0d0**(-1.0d0)*C_NaCl

        TT=20.0d0-TT
        muwSC=1002.0d0*10**((TT*(1.2378d0-1.303d0*10.0d0**(-3.0d0)*TT+3.06d0*10.0d0**(-6.0d0)*TT**2.0d0  &
                            +2.55d0*10.0d0**(-8.0d0)*TT**3.0d0))/(96.0d0+T0-273.15d0))

        muw=muwSC*10.0d0**(Av+Bv*log10(muwSC/1002.0d0))

        mu(4)=muw*10.0d0**(-6.0d0)*(1.0d0+betav*ptotal*10.0d0**(-9.0d0))
        mu(4)=mu(4)*1000.0d0



!!foam相の相粘度
        FM2=((fmoil-So)/(fmoil-Sor))**epoil
        FM7=0.5d0+atan(epdry*(Sw-fmdry))/(4.0d0*atan(1.0d0))

        !FM=1.0d0/(1.0d0+fmmob*FM1*FM7)  !F2未導入
        FM=1.0d0/(1.0d0+fmmob*FM1)
        mu(2)=mu(1)*10000.0d0!/FM**epvis

!!単位換算 (cP)→(Pa*s)

        do phase=1,s
            mu(phase)=mu(phase)*10**(-3.0d0) !(Pa・s)
        end do
        
    end subroutine calc_viscosity
end module