module mod_density
    use mod_initial
    use mod_autodiff
    use mod_interaction
    use mod_cubic_eq
    implicit none
    contains 
    subroutine calc_density(unknown,Ptotal,Xpc,rho,Vp)
        implicit none
        integer,intent(in)::unknown
        type(diffs),intent(in)::Ptotal
        type(diffs),intent(in),dimension(s,t)::Xpc
        type(diffs),intent(out),dimension(s)::rho,Vp
        real(8),allocatable,dimension(:,:)::Theta
        real(8),allocatable,dimension(:)::Pc,Tc,Vc,omega,Tr,acc,bcc,alf,ac,bc  &
                                        ,Aco,Bco,Cco,zfactor
        real(8)::R,Arc,Brc,Crc,Drc,Erc,Frc,Grc,Hrc,Vco2,Xpc_surf,Xpc_NaCl
        real(8)::T00
        integer::p,c,cc,i
        type(diffs)::X(m*n),Pp(m),Pr(m,n),ap(m),bp(m),LA(m),LB(m),A_coef(m),B_coef(m),C_coef(m),zfac(m)  &
                    ,MWw,Vw,Vl,W_s,W_NaCl,W_salt,M_salt,P00,rho_brine,rhow,M_brine,vl_brine,Xpc_gas,Vf

        allocate(Pc(n),Tc(n),Vc(n),omega(n),Tr(n),acc(n),bcc(n),alf(n),ac(n),bc(n),Theta(n,n)  &
                ,Aco(m),Bco(m),Cco(m),zfactor(m))

!気相, 油相の相モル密度 rho(mol/m^3)
    !Zfactor
        do c=1,n
            X(c)=Xpc(1,c)
        end do

        do p=2,m
            do c=1,n
                X(n*(p-1)+c)=Xpc(p+1,c)
            end do
        end do

     !fugacity
        Pc(1)=Pc1  
        Pc(2)=Pc2  
        Pc(3)=Pc3
        Pc(4)=Pc4
        Pc(5)=Pc5  
        Pc(6)=Pc6  
        Pc(7)=Pc7
        Pc(8)=Pc8
        Pc(9)=Pc9
  
        Tc(1)=Tc1  
        Tc(2)=Tc2  
        Tc(3)=Tc3
        Tc(4)=Tc4 
        Tc(5)=Tc5 
        Tc(6)=Tc6  
        Tc(7)=Tc7
        Tc(8)=Tc8
        Tc(9)=Tc9
        
        Vc(1)=Vc1
        Vc(2)=Vc2 
        Vc(3)=Vc3 
        Vc(4)=Vc4
        Vc(5)=Vc5
        Vc(6)=Vc6 
        Vc(7)=Vc7 
        Vc(8)=Vc8
        Vc(9)=Vc9
  
        omega(1)=omega1  
        omega(2)=omega2  
        omega(3)=omega3 
        omega(4)=omega4 
        omega(5)=omega5  
        omega(6)=omega6  
        omega(7)=omega7 
        omega(8)=omega8
        omega(9)=omega9
      
  ! Pp(p)=p相における圧力 (Pa)
        Pp(1)=Ptotal  
        Pp(2)=Ptotal   
        Pp(3)=Ptotal  
  
  !気体定数 (Pa*m^3/K)
        R=8.3144598d0 !(Pa*m^3/K) 
  
  ! (p,c)=(phase,component)=(m,n)
  
  ! Pr(p,c) = p相におけるc成分の還元圧力 (-)
        do p=1,m
              do c=1,n
                  Pr(p,c)=Pp(p)/Pc(c)    
              end do
          end do
      
      ! Tr(c) = p相におけるc成分の還元温度 (-)
          do c=1,n
              Tr(c)=T0/Tc(c)   
          end do
  
  !acc(c) = 臨界点における、c成分の未定係数a
  !bcc(c) = 臨界点における、c成分の未定係数b
          do c=1,n
              acc(c)=0.45724d0*(R**2.0d0*Tc(c)**2.0d0)/Pc(c)
              bcc(c)=0.07780d0*R*Tc(c)/Pc(c)                 
          end do
      
          do c=1,n
              alf(c)=(1.0d0+(0.37464d0+1.54226d0*omega(c)-0.26992d0*omega(c)**2.0d0)*(1.0d0-dsqrt(Tr(c))))**2.0d0      
          end do
      
          do c=1,n
              ac(c)=acc(c)*alf(c)  
              bc(c)=bcc(c)         
          end do
  
          
  !i成分とj成分における相互作用係数(-)
          call calc_interaction(Vc,Theta)
      
          Theta(1,8)=Theta18
          Theta(8,1)=Theta(1,8)
          Theta(1,9)=Theta19
          Theta(9,1)=Theta(1,9)
          Theta(2,8)=Theta28
          Theta(8,2)=Theta(2,8)
          Theta(2,9)=Theta29
          Theta(9,2)=Theta(2,9)
          Theta(3,8)=Theta38
          Theta(8,3)=Theta(3,8)
          Theta(3,9)=Theta39
          Theta(9,3)=Theta(3,9)
          Theta(4,8)=Theta48
          Theta(8,4)=Theta(4,8)
          Theta(4,9)=Theta49
          Theta(9,4)=Theta(4,9)
          Theta(5,8)=Theta58
          Theta(8,5)=Theta(5,8)
          Theta(5,9)=Theta59
          Theta(9,5)=Theta(5,9)
          Theta(6,8)=Theta68
          Theta(8,6)=Theta(6,8)
          Theta(6,9)=Theta69
          Theta(9,6)=Theta(6,9)
          Theta(7,8)=Theta78
          Theta(8,7)=Theta(7,8)
          Theta(7,9)=Theta79 
          Theta(9,7)=Theta(7,9) 
  
          Theta(8,9)=Theta89
          Theta(9,8)=Theta(8,9)
  
          do p=1,m
              call residualvectorset3(unknown*grid+BHin+BHout,ap(p))
              call residualvectorset3(unknown*grid+BHin+BHout,bp(p))
          end do
  
          do p=1,m
              do c=1,n
                  do cc=1,n
                      ap(p)=ap(p)+X(n*(p-1)+c)*X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))               
                  end do
              end do
          end do
      
          do p=1,m
              do c=1,n
                bp(p)=bp(p)+X(n*(p-1)+c)*bc(c)    
              end do
          end do
  
        do p=1,m
            LA(p)=ap(p)*Pp(p)/(R*T0)**2.0d0    
            LB(p)=bp(p)*Pp(p)/(R*T0)           
        end do
  
        do p=1,m
            A_coef(p)=-1.0d0+LB(p)                             
            B_coef(p)=LA(p)-3.0d0*LB(p)**2.0d0-2.0d0*LB(p)      
            C_coef(p)=LB(p)**3.0d0-LA(p)*LB(p)+LB(p)**2.0d0    
        end do
  
        call outxs(A_coef,Aco)
        call outxs(B_coef,Bco)
        call outxs(C_coef,Cco)
  
        call cubic_eq_other(Aco(1),Bco(1),Cco(1),1,zfactor(1))
        call cubic_eq_other(Aco(2),Bco(2),Cco(2),2,zfactor(2))
        call cubic_eq_other(Aco(3),Bco(3),Cco(3),3,zfactor(3))
  
        do p=1,m
            call residualvectorset4(zfactor(p),unknown*grid+BHin+BHout,zfac(p))
        end do
  
        rho(1)=Ptotal/(zfactor(1)*R*T0)

        rho(3)=Ptotal/(zfactor(2)*R*T0)

!水相の相モル密度 rho(mol/m^3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!玲音さんver
        !Rowe and Chouの式
        Arc=5.916365d0-0.01035794d0*T0+0.9270048d0*10.0d0**(-5.0d0)*T0**2.0d0-1127.522d0/T0+100674.1d0/(T0**2.0d0) 
        Brc=0.5204914d0*10.0d0**(-2.0d0)-0.10482101d0*10.0d0**(-4.0d0)*T0+0.8328532d0*10.0d0**(-8.0d0)*T0**(2.0d0)  &
           -1.1702939d0/T0+102.2783d0/T0**(2.0d0)
        Crc=0.118547d0*10.0d0**(-7.0d0)-0.6599143d0*10.0d0**(-10.0d0)*T0
        Drc=-2.5166d0+0.0111766d0*T0-0.170533d0*10.0d0**(-4.0d0)*T0**2.0d0
        Erc=2.84851d0-0.0154305d0*T0+0.228982d0*10.0d0**(-4.0d0)*T0**2.0d0
        Frc=-0.0014814d0-0.82969d0*10.0d0**(-5.0d0)*T0-0.12469d0*10.0d0**(-7.0d0)*T0**2.0d0
        Grc=0.0027141d0-0.15391d0*10.0d0**(-4.0d0)*T0+0.22655d0*10**(-7.0d0)*T0**2.0d0
        Hrc=0.62158d0*10.0d0**(-6.0d0)-0.40075d0*10.0d0**(-8.0d0)*T0+0.65972d0*10.0d0**(-11.0d0)*T0**2.0d0


        MWw=MW9*10.0d0**(-3.0d0)*MW11*10.0d0**(-3.0d0)/(Xpc(4,11)*MW9*10.0d0**(-3.0d0)+(1.0d0-Xpc(4,11))*MW11*10.0d0**(-3.0d0))

        Vw=MWw*10**(-3.0d0)*(Arc-Brc*(Ptotal/(9.8d0*10.0d0**4.0d0))-Crc*(Ptotal/(9.8d0*10**4.0d0))**2.0d0  &
          +Xpc(4,11)*Drc+Xpc(4,11)**2.0d0*Erc-Xpc(4,11)*Frc*(Ptotal/(9.8d0*10.0d0**4.0d0))  &
          -Xpc(4,11)**2.0d0*Grc*(Ptotal/(9.8d0*10.0d0**4.0d0))-Xpc(4,11)*Hrc/2.0d0*(Ptotal/(9.8d0*10**4.0d0))**2.0d0)

        !Duan and Sunの式
        Vco2=(-47.7518d0+4.336154d0*10.0d0**(-1.0d0)*T0-5.945771d0*10**(-4.0d0)*T0**2.0d0)/10**(6.0d0)

        Vl=Xpc(4,8)*Vco2+(1.0d0-Xpc(4,8))*Vw

        !rho(4)=1.0d0/Vl

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!吉田さんver
        W_s=Xpc(4,10)*MW10/(Xpc(4,9)*MW9+Xpc(4,10)*MW10+Xpc(4,11)*MW11)
        W_NaCl=Xpc(4,11)*MW11/(Xpc(4,9)*MW9+Xpc(4,10)*MW10+Xpc(4,11)*MW11)

        W_salt=W_s+W_NaCl

        P00=Ptotal*10.0d0**(-6.0d0)  !MPa
        T00=T0-273.15d0  !℃
        rhow=1.0d0+10.0d0**(-6.0d0)*(-80.0d0*T00-3.3d0*T00**2.0d0+0.00175d0*T00**3.0d0+489.0d0*P00-2.0d0*T00*P00  &
             +0.016d0*T00**2.0d0*P00-1.3d0*10.0d0**(-5.0d0)*T00**3.0d0*P00-0.333d0*P00**2.0d0-0.002d0*T00*P00**2.0d0)

        rho_brine=rhow+W_salt*(0.668d0+0.44d0*W_salt+10.0d0**(-6.0d0)*(300.0d0*P00-2400.0d0*P00*W_salt  &
                  +T00*(80.0d0+3.0d0*T00-3300.0d0*W_salt-13.0d0*P00+47.0d0*P00*W_salt)))

        call out_diffsx(Xpc(4,10),Xpc_surf)
        call out_diffsx(Xpc(4,11),Xpc_NaCl)

        if (Xpc_surf == 0 .AND. Xpc_NaCl == 0) then
            call residualvectorset4(Mw9*10.0d0**(-3.0d0),unknown*grid+BHin+BHout,M_brine)
        else
            M_salt=(Xpc(4,10)*MW10+Xpc(4,11)*MW11)/(Xpc(4,10)+Xpc(4,11))
            M_brine=Mw9*M_salt/((1.0d0-W_salt)*M_salt+W_salt*Mw9)
            M_brine=M_brine*10.0d0**(-3.0d0)
        end if

        vl_brine=M_brine/(rho_brine*10.0d0**3.0d0)

        !Duan and Sunの式
        Vco2=(-47.7518d0+4.336154d0*10.0d0**(-1.0d0)*T0-5.945771d0*10**(-4.0d0)*T0**2.0d0)/10**(6.0d0)

        Vl=Xpc(4,8)*Vco2+(Xpc(4,9)+Xpc(4,10)+Xpc(4,11))*vl_brine

        rho(4)=1.0d0/Vl

!!foam相の相モル密度 rho(mol/m^3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!立花ver
        !!CO2+HC のモル密度

        do i=1,m*n
            call residualvectorset3(unknown*grid+BHin+BHout,X(i))
        end do

        call residualvectorset3(unknown*grid+BHin+BHout,Xpc_gas)
        do c=1,8
            Xpc_gas=Xpc_gas+Xpc(2,c)
        end do
        
        do c=1,8
            X(c)=Xpc(2,c)/Xpc_gas
        end do
        
      
        !fugacity
        Pc(1)=Pc1  
        Pc(2)=Pc2  
        Pc(3)=Pc3
        Pc(4)=Pc4
        Pc(5)=Pc5  
        Pc(6)=Pc6  
        Pc(7)=Pc7
        Pc(8)=Pc8
        Pc(9)=Pc9
        
  
        Tc(1)=Tc1  
        Tc(2)=Tc2  
        Tc(3)=Tc3
        Tc(4)=Tc4 
        Tc(5)=Tc5 
        Tc(6)=Tc6  
        Tc(7)=Tc7
        Tc(8)=Tc8
        Tc(9)=Tc9
        
        Vc(1)=Vc1
        Vc(2)=Vc2 
        Vc(3)=Vc3 
        Vc(4)=Vc4
        Vc(5)=Vc5
        Vc(6)=Vc6 
        Vc(7)=Vc7 
        Vc(8)=Vc8
        Vc(9)=Vc9
  
        omega(1)=omega1  
        omega(2)=omega2  
        omega(3)=omega3 
        omega(4)=omega4 
        omega(5)=omega5  
        omega(6)=omega6  
        omega(7)=omega7 
        omega(8)=omega8
        omega(9)=omega9
      
  ! Pp(p)=p相における圧力 (Pa)
        Pp(1)=Ptotal  
        Pp(2)=Ptotal   
        Pp(3)=Ptotal  
  
  !気体定数 (Pa*m^3/K)
        R=8.3144598d0 !(Pa*m^3/K) 
  
  ! (p,c)=(phase,component)=(m,n)
  
  ! Pr(p,c) = p相におけるc成分の還元圧力 (-)
        do p=1,m
              do c=1,n
                  Pr(p,c)=Pp(p)/Pc(c)    
              end do
          end do
      
      ! Tr(c) = p相におけるc成分の還元温度 (-)
          do c=1,n
              Tr(c)=T0/Tc(c)   
          end do
  
  !acc(c) = 臨界点における、c成分の未定係数a
  !bcc(c) = 臨界点における、c成分の未定係数b
          do c=1,n
              acc(c)=0.45724d0*(R**2.0d0*Tc(c)**2.0d0)/Pc(c)
              bcc(c)=0.07780d0*R*Tc(c)/Pc(c)                 
          end do
      
          do c=1,n
              alf(c)=(1.0d0+(0.37464d0+1.54226d0*omega(c)-0.26992d0*omega(c)**2.0d0)*(1.0d0-dsqrt(Tr(c))))**2.0d0      
          end do
      
          do c=1,n
              ac(c)=acc(c)*alf(c)  
              bc(c)=bcc(c)         
          end do
  
          
  !i成分とj成分における相互作用係数(-)
          call calc_interaction(Vc,Theta)
      
          Theta(1,8)=Theta18
          Theta(8,1)=Theta(1,8)
          Theta(1,9)=Theta19
          Theta(9,1)=Theta(1,9)
          Theta(2,8)=Theta28
          Theta(8,2)=Theta(2,8)
          Theta(2,9)=Theta29
          Theta(9,2)=Theta(2,9)
          Theta(3,8)=Theta38
          Theta(8,3)=Theta(3,8)
          Theta(3,9)=Theta39
          Theta(9,3)=Theta(3,9)
          Theta(4,8)=Theta48
          Theta(8,4)=Theta(4,8)
          Theta(4,9)=Theta49
          Theta(9,4)=Theta(4,9)
          Theta(5,8)=Theta58
          Theta(8,5)=Theta(5,8)
          Theta(5,9)=Theta59
          Theta(9,5)=Theta(5,9)
          Theta(6,8)=Theta68
          Theta(8,6)=Theta(6,8)
          Theta(6,9)=Theta69
          Theta(9,6)=Theta(6,9)
          Theta(7,8)=Theta78
          Theta(8,7)=Theta(7,8)
          Theta(7,9)=Theta79 
          Theta(9,7)=Theta(7,9) 
  
          Theta(8,9)=Theta89
          Theta(9,8)=Theta(8,9)
  
          do p=1,m
              call residualvectorset3(unknown*grid+BHin+BHout,ap(p))
              call residualvectorset3(unknown*grid+BHin+BHout,bp(p))
          end do
  
          do p=1,m
              do c=1,n
                  do cc=1,n
                      ap(p)=ap(p)+X(n*(p-1)+c)*X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))               
                  end do
              end do
          end do
      
          do p=1,m
              do c=1,n
                bp(p)=bp(p)+X(n*(p-1)+c)*bc(c)    
              end do
          end do
  
        do p=1,m
            LA(p)=ap(p)*Pp(p)/(R*T0)**2.0d0    
            LB(p)=bp(p)*Pp(p)/(R*T0)           
        end do
  
        do p=1,m
            A_coef(p)=-1.0d0+LB(p)                             
            B_coef(p)=LA(p)-3.0d0*LB(p)**2.0d0-2.0d0*LB(p)      
            C_coef(p)=LB(p)**3.0d0-LA(p)*LB(p)+LB(p)**2.0d0    
        end do
  
        call outxs(A_coef,Aco)
        call outxs(B_coef,Bco)
        call outxs(C_coef,Cco)
  
        call cubic_eq_other(Aco(1),Bco(1),Cco(1),1,zfactor(1))

        call residualvectorset4(zfactor(1),unknown*grid+BHin+BHout,zfac(1))
   
        Vf=zfac(1)*R*T0/Ptotal

        !!水とSurfactantのみのモル密度
        W_s=Xpc(4,10)*MW10/(Xpc(4,9)*MW9+Xpc(4,10)*MW10)
        call residualvectorset3(unknown*grid+BHin+BHout,W_NaCl)

        W_salt=W_s+W_NaCl

        P00=Ptotal*10.0d0**(-6.0d0)  !MPa
        T00=T0-273.15d0  !℃
        rhow=1.0d0+10.0d0**(-6.0d0)*(-80.0d0*T00-3.3d0*T00**2.0d0+0.00175d0*T00**3.0d0+489.0d0*P00-2.0d0*T00*P00  &
             +0.016d0*T00**2.0d0*P00-1.3d0*10.0d0**(-5.0d0)*T00**3.0d0*P00-0.333d0*P00**2.0d0-0.002d0*T00*P00**2.0d0)

        rho_brine=rhow+W_salt*(0.668d0+0.44d0*W_salt+10.0d0**(-6.0d0)*(300.0d0*P00-2400.0d0*P00*W_salt  &
                  +T00*(80.0d0+3.0d0*T00-3300.0d0*W_salt-13.0d0*P00+47.0d0*W_salt)))

        call residualvectorset4(MW10,unknown*grid+BHin+BHout,M_salt)
        
        M_brine=Mw9*M_salt/((1.0d0-W_salt)*M_salt+W_salt*Mw9)
        M_brine=M_brine*10.0d0**(-3.0d0)

        vl_brine=M_brine/(rho_brine*10.0d0**3.0d0)

        !!上記二つよりfoam相における相モル密度を以下に示す。

        rho(2)=1.0d0/(Xpc_gas*Vf+(1.0d0-Xpc_gas)*vl_brine)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ideal mixing ver
        
        Vp(1)=1.0d0/rho(1)
        Vp(2)=1.0d0/rho(2)
        Vp(3)=1.0d0/rho(3)
        Vp(4)=1.0d0/rho(4)
        
    end subroutine calc_density
end module