module mod_ORR
    use mod_initial
    use mod_autodiff
    use mod_interaction
    use mod_cubic_eq
    implicit none
    contains
    subroutine calc_ORR(unknown,Cgc,Cfc,Coc,Cwc,rho,qout,Vgsc0,Vosc0)
        integer,intent(in)::unknown
        type(diffs),intent(in),dimension(t,grid)::Cgc,Cfc,Coc,Cwc
        type(diffs),intent(in),dimension(s,grid)::rho
        type(diffs),intent(in),dimension(s)::qout
        real(8),intent(out)::Vgsc0,Vosc0
        integer::com,phase,c,cc
        real(8)::R
        real(8),allocatable,dimension(:,:)::Pr,Theta
        real(8),allocatable,dimension(:)::Pp,Pc,Tc,Vc,omega,Tr,ac,bc,acc,bcc,alf  &
                                        ,Aco,Bco,Cco,zfactor
        type(diffs)::Cc_pro(s,t),dNc(t),dNp(s),Cgctotal,Cfctotal,Coctotal,Cwctotal
        type(diffs)::X(s*n),ap(s),bp(s),LA(s),LB(s),A_coef(s),B_coef(s),C_coef(s),Vsc(s),Vgsc,Vosc

        allocate(Pp(s),Pr(s,n),Pc(n),Tc(n),Vc(n),omega(n),Tr(n),ac(n),bc(n)  &
            ,acc(n),bcc(n),alf(n),Theta(n,n),Aco(s),Bco(s),Cco(s),zfactor(s))

        do com=1,t
            Cc_pro(1,com)=Cgc(com,grid)
            Cc_pro(2,com)=Cfc(com,grid)
            Cc_pro(3,com)=Coc(com,grid)
            Cc_pro(4,com)=Cwc(com,grid)
        end do

        do phase=1,s
            dNp(phase)=rho(phase,grid)*(-1.0d0)*qout(phase)
        end do

        do com=1,t
            call residualvectorset3(unknown*grid+BHin+BHout,dNc(com))
        end do
        do com=1,t
            do phase=1,s
                dNc(com)=dNc(com)+dNp(phase)*Cc_pro(phase,com)
            end do
        end do

        call residualvectorset3(unknown*grid+BHin+BHout,Cgctotal)
        call residualvectorset3(unknown*grid+BHin+BHout,Cfctotal)
        call residualvectorset3(unknown*grid+BHin+BHout,Coctotal)
        call residualvectorset3(unknown*grid+BHin+BHout,Cwctotal)   
        do c=1,n
            Cgctotal=Cgctotal+Cgc(c,grid)
            Cfctotal=Cfctotal+Cfc(c,grid)
            Coctotal=Coctotal+Coc(c,grid)
            Cwctotal=Cwctotal+Cwc(c,grid)
        end do

        do c=1,n 
            X(c)=Cgc(c,grid)/Cgctotal
            X(c+n)=Cfc(c,grid)/Cfctotal
            X(c+2*n)=Coc(c,grid)/Coctotal
            X(c+3*n)=Cwc(c,grid)/Cwctotal
        end do

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
        Pp(1)=Psc  
        Pp(2)=Psc   
        Pp(3)=Psc 
        Pp(4)=Psc
  
  !気体定数 (Pa*m^3/K)
        R=8.3144598d0 !(Pa*m^3/(K*mol)) 
  
  
  ! Pr(p,c) = p相におけるc成分の還元圧力 (-)
        do phase=1,s
              do c=1,n
                  Pr(phase,c)=Pp(phase)/Pc(c)    
              end do
        end do
      
      ! Tr(c) = p相におけるc成分の還元温度 (-)
          do c=1,n
              Tr(c)=Tsc/Tc(c)   
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
  
          do phase=1,s
              call residualvectorset3(unknown*grid+BHin+BHout,ap(phase))
              call residualvectorset3(unknown*grid+BHin+BHout,bp(phase))
          end do
  
          do phase=1,s
              do c=1,n
                  do cc=1,n
                      ap(phase)=ap(phase)+X(n*(phase-1)+c)*X(n*(phase-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))               
                  end do
              end do
          end do
      
          do phase=1,s
              do c=1,n
                bp(phase)=bp(phase)+X(n*(phase-1)+c)*bc(c)    
              end do
          end do
  
        do phase=1,s
            LA(phase)=ap(phase)*Pp(phase)/(R*Tsc)**2.0d0    
            LB(phase)=bp(phase)*Pp(phase)/(R*Tsc)           
        end do
  
        do phase=1,s
            A_coef(phase)=-1.0d0+LB(phase)                             
            B_coef(phase)=LA(phase)-3.0d0*LB(phase)**2.0d0-2.0d0*LB(phase)      
            C_coef(phase)=LB(phase)**3.0d0-LA(phase)*LB(phase)+LB(phase)**2.0d0    
        end do
  
        call outxs(A_coef,Aco)
        call outxs(B_coef,Bco)
        call outxs(C_coef,Cco)
  
        call cubic_eq(Aco(1),Bco(1),Cco(1),0.9951d0,zfactor(1))
        call cubic_eq(Aco(2),Bco(2),Cco(2),0.0098d0,zfactor(2))
        call cubic_eq(Aco(3),Bco(3),Cco(3),0.0098d0,zfactor(3))
        call cubic_eq(Aco(4),Bco(4),Cco(4),0.0008d0,zfactor(4))
 
        write(*,*)"zfacter",Zfactor(1),zfactor(3)
        do phase=1,s
            Vsc(phase)=(zfactor(phase)*dNp(phase)*R*Tsc)/Psc
        end do
        call out_diffsx(Vsc(1),Vgsc0)
        call out_diffsx(Vsc(3),Vosc0)

        call residualvectorset3(unknown*grid+BHin+BHout,Vgsc)
        call residualvectorset3(unknown*grid+BHin+BHout,Vosc)
    !#TODOガスはモルと体積は比例関係にあるが、油はそうではない。それを考慮できていない。
        do c=1,7
            Vgsc=Vgsc+Cc_pro(1,c)*Vsc(1)
            Vosc=Vosc+Cc_pro(3,c)*Vsc(3)
        end do
        
        call out_diffsx(Vgsc,Vgsc0)
        call out_diffsx(Vosc,Vosc0)
        write(*,*)"gas",Vgsc0
        write(*,*)"oil",Vosc0
    end subroutine calc_ORR
end module mod_ORR