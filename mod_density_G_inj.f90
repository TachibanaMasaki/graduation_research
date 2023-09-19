module mod_density_G_inj
    use mod_initial
    use mod_autodiff
    use mod_interaction
    use mod_cubic_eq
    implicit none
    contains 
    subroutine calc_density_G_inj(unknown,Po1,rho_inj)
        integer,intent(in)::unknown
        type(diffs),intent(in)::Po1
        type(diffs),intent(out)::rho_inj
        real(8),allocatable,dimension(:,:)::Theta
        real(8),allocatable,dimension(:)::Pc,Tc,Vc,omega,Tr,acc,bcc,alf,ac,bc  &
                                        ,Aco,Bco,Cco,zfactor
        real(8)::R,rho_inj_real
        integer::p,c,cc,i
        type(diffs)::X(m*n),Pp(m),Pr(m,n),ap(m),bp(m),LA(m),LB(m),A_coef(m),B_coef(m),C_coef(m),zfac(m),Vg

        allocate(Pc(n),Tc(n),Vc(n),omega(n),Tr(n),acc(n),bcc(n),alf(n),ac(n),bc(n),Theta(n,n)  &
                ,Aco(m),Bco(m),Cco(m),zfactor(m))


            do i=1,m*n
                call residualvectorset3(unknown*grid+BHin+BHout,X(i))
            end do
        
            call residualvectorset4(1.0d0,unknown*grid+BHin+BHout,X(8))
        
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
            Pp(1)=Po1  
            Pp(2)=Po1   
            Pp(3)=Po1  
          
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
    
   
            Vg=zfac(1)*R*T0/Po1 
            rho_inj=1.0d0/Vg

            call out_diffsx(rho_inj,rho_inj_real)

    end subroutine calc_density_G_inj
end module