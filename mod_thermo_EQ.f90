module mod_Thermo_EQ
    use mod_initial
    use mod_autodiff
    use mod_interaction
    use mod_cubic_eq
    implicit none
    contains 
    subroutine calc_Thermo_EQ(unknown,Ptotal,n0,V0,L0,W0,lnKgo,lnKwo,g)
      implicit none
      integer,intent(in)::unknown
      real(8),intent(in)::V0,L0,W0
      type(diffs),intent(in),dimension(n)::n0,lnKgo,lnKwo
      type(diffs),intent(in)::Ptotal
      type(diffs),intent(out),dimension((m-1)*n)::g(:)

      real(8),allocatable,dimension(:,:)::Theta
      real(8),allocatable,dimension(:)::lnK0,Tb,Pc,Tc,Vc,omega,Tr,ac,bc,acc,bcc,alf  &
                                        ,Aco,Bco,Cco,zfactor
      integer::p,c,cc
      real(8)::R,Arc,Brc,Crc,Drc,Erc,Frc,Grc,Hrc,MWw,Aant,Bant,Cant,Ahar,Bhar,Char  &
               ,Vco2,Psw,fsw,lnHsco2,Hsco2

      type(diffs)::ntotal,z0(n),lnK((m-1)*n),X(m*n),Pp(m),Pr(m,n),ap(m),bp(m),LA(m),LB(m)  &
                  ,A_coef(m),B_coef(m),C_coef(m),zfac(m),middle(m,n),lnf(m,n)  &
                  ,Vw,lnHco2,Hco2,lnf31,lnf33

      allocate(lnK0((m-1)*n),Tb(n),Pc(n),Tc(n),Vc(n),omega(n),Tr(n),ac(n),bc(n)  &
              ,acc(n),bcc(n),alf(n),Theta(n,n),Aco(m),Bco(m),Cco(m),zfactor(m))

      call residualvectorset3(unknown*grid+BHin+BHout,ntotal)
      do c=1,n
        ntotal=ntotal+n0(c)
      end do

      do c=1,n
        z0(c)=n0(c)/ntotal
      end do

      do c=1,n 
        lnK(c)=lnKgo(c)
        lnK(c+n)=lnKwo(c)
      end do


      do c=1,n 
          X(c)=exp(lnK(c))*z0(c)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
          X(c+n)=z0(c)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
          X(c+2*n)=exp(lnK(c+n))*z0(c)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
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

      call cubic_eq(Aco(1),Bco(1),Cco(1),0.7322d0,zfactor(1))
      call cubic_eq(Aco(2),Bco(2),Cco(2),0.4640d0,zfactor(2))
      call cubic_eq(Aco(3),Bco(3),Cco(3),0.0868d0,zfactor(3))

      do p=1,m
          call residualvectorset4(zfactor(p),unknown*grid+BHin+BHout,zfac(p))
      end do

      !水相のzfactorの計算の上書き
!Rowe and Chouが提唱した値
      Arc=5.916365d0-0.01035794d0*T0+0.9270048d0*10.0d0**(-5.0d0)*T0**2.0d0-1127.522d0/T0+100674.1d0/(T0**2.0d0) !-0.001035794説もあり
      Brc=0.5204914d0*10.0d0**(-2.0d0)-0.10482101d0*10.0d0**(-4.0d0)*T0+0.8328532d0*10.0d0**(-8.0d0)*T0**(2.0d0)  &
         -1.1702939d0/T0+102.2783d0/T0**(2.0d0)
      Crc=0.118547d0*10.0d0**(-7.0d0)-0.6599143d0*10.0d0**(-10.0d0)*T0
      Drc=-2.5166d0+0.0111766d0*T0-0.170533d0*10.0d0**(-4.0d0)*T0**2.0d0
      Erc=2.84851d0-0.0154305d0*T0+0.228982d0*10.0d0**(-4.0d0)*T0**2.0d0
      Frc=-0.0014814d0-0.82969d0*10.0d0**(-5.0d0)*T0-0.12469d0*10.0d0**(-7.0d0)*T0**2.0d0
      Grc=0.0027141d0-0.15391d0*10.0d0**(-4.0d0)*T0+0.22655d0*10**(-7.0d0)*T0**2.0d0
      Hrc=0.62158d0*10.0d0**(-6.0d0)-0.40075d0*10.0d0**(-8.0d0)*T0+0.65972d0*10.0d0**(-11.0d0)*T0**2.0d0


      MWw=MW3*10**(-3.0d0) 

      Vw=MWw*10**(-3.0d0)*(Arc-Brc*(Ptotal/(9.8d0*10.0d0**4.0d0))-Crc*(Ptotal/(9.8d0*10**4.0d0))**2.0d0)

      
!!↓外すとペンロビンソンバージョン
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      zfac(3)=Ptotal*Vw/(R*T0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!フガシティ計算 (気相、油相、水相)
      do p=1,m
            do c=1,n
                call residualvectorset3(unknown*grid+BHin+BHout,middle(p,c))
            end do
      end do
        do p=1,m
            do c=1,n
                do cc=1,n
                    middle(p,c)=middle(p,c)+X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))  
                end do
            end do
        end do

        do p=1,m
            do c=1,n
                lnf(p,c)=bc(c)/bp(p)*(zfac(p)-1.0d0)-log(zfac(p)-LB(p))  &
                         -LA(p)/(2.0d0*dsqrt(2.0d0)*LB(p))*(2.0d0*middle(p,c)/ap(p)-bc(c)/bp(p))  &
                         *log((zfac(p)+LB(p)*(1.0d0+dsqrt(2.0d0)))/(zfac(p)+LB(p)*(1.0d0-dsqrt(2.0d0))))
            end do
        end do

!!フガシティ計算 (水相)
!水相の二酸化炭素のフガシティ係数 (lnHco2はヘンリー定数の対数)
!HarveyによるCO2成分に関するヘンリー定数の補正
        Ahar=-9.4234d0
        Bhar=4.0087d0
        Char=10.3199d0
        
    !Antnieの式における係数
        Aant=8.02754d0  !(-)
        Bant=1705.616d0 !(K)
        Cant=-41.745d0  !(K)
    
    !水の飽和蒸気圧 (Pa)
        Psw=133.322368d0*10**(Aant-(Bant/(T0+Cant)))

    !CO2の液相モル体積Vco2はDuan and Sunの式を用いる ヘンリー定数の単位はMPa
        Vco2=(-47.7518d0+4.336154d0*10**(-1.0d0)*T0-5.945771d0*10**(-4.0d0)*T0**2.0d0)/(10**(6.0d0))
        lnHsco2=log(Psw)+Ahar/Tr(3)+(Bhar*(1.0d0-Tr(3))**(0.355d0))/Tr(3)+(Char*dexp(1.0d0-Tr(3)))/(Tr(3)**(0.41d0))
        Hsco2=dexp(lnHsco2)
        lnHco2=lnHsco2+(Vco2/(R*T0))*(Ptotal-Psw) 
        Hco2=exp(lnHco2)
        lnf31=log(Hco2/(Ptotal))
    
!!↓外すとペンロビンソンバージョン
!        call residualvectorset4(lnf31,unknown*grid+BHin+BHout,lnf(3,1))
    
    !飽和純水のフガシティ係数 (-)
        if (T0>303.7) then
            fsw=1.0d0
        else
            fsw=0.9958d0+9.68330d0*(1.8d0*T0-459.67d0)*10**(-5.0d0)  &
                -6.1750d0*((1.8d0*T0-459.67d0)**2.0d0)*10**(-7.0d0)  &
                -3.08333d0*((1.8d0*T0-459.67d0)**3.0d0)*10**(-10.0d0)
        end if
      
     
        lnf33=log((Psw*fsw/(Ptotal))*exp(MWw/(10**(3.0d0)*R*T0)*(Arc*(Ptotal-Psw)-0.5d0/(9.8d0*10**(4.0d0))*Brc  &
                  *((Ptotal)**2.0d0-Psw**2.0d0)-1.0d0/3.0d0*(1.0d0/(9.8d0*10**(4.0d0)))**2.0d0*Crc*((Ptotal)**3.0d0-Psw**3.0d0))))
!!↓外すとペンロビンソンバージョン
!      call residualvectorset4(lnf33,unknown*grid+BHin+BHout,lnf(3,3))  

      do c=1,n
        g(c)=lnK(c)+lnf(1,c)-lnf(2,c)
        g(c+n)=lnK(c+n)+lnf(3,c)-lnf(2,c)
      end do

    end subroutine calc_Thermo_EQ
  end module mod_Thermo_EQ