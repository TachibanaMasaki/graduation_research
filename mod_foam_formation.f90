module mod_foam_formation
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_foam_formation(unknown,Xf0pc,Xfpc,Xpc,V,F,L,W,FM1)
      implicit none
      integer,intent(in)::unknown
      type(diffs),intent(in),dimension(m,t)::Xf0pc
      type(diffs),intent(out),dimension(s,t)::Xfpc,Xpc
      type(diffs),intent(inout)::V,F,L,W
      type(diffs),intent(out)::FM1
      real(8)::Xfpc_H2O,Xfpc_surf,Xfpc_gas,Cs,FM1_real,alpha_real,beta_real,V0,f0,L0,W0
      integer::phase,com,c
      type(diffs)::Mg_f,alpha,beta,Xfpc_gas_diff
    
!1~7までがoil, 8:CO2, 9:H2O, 10:Surfactant, 11:NaCl
    
      do phase=1,s
          do com=1,t
              call residualvectorset3(unknown*grid+BHin+BHout,Xfpc(phase,com))
          end do
      end do

      do com=1,t
        Xfpc(1,com)=Xf0pc(1,com)  !gas Phase
                                  !foam Phase
        Xfpc(3,com)=Xf0pc(2,com)  !oil Phase
        Xfpc(4,com)=Xf0pc(3,com)  !water Pahse
      end do

!一部のHCもfoamにする。
      call residualvectorset3(unknown*grid+BHin+BHout,Xfpc_gas_diff)
      do c=1,7
          Xfpc_gas_diff=Xfpc_gas_diff+fmHC*Xfpc(1,c)
      end do
      Xfpc_gas_diff=Xfpc_gas_diff+Xfpc(1,8)
      
      call out_diffsx(Xfpc_gas_diff,Xfpc_gas)
      write(*,*)Xfpc_gas

      call out_diffsx(Xfpc(4,9),Xfpc_H2O)
      call out_diffsx(Xfpc(4,10),Xfpc_surf)

      
      if (Xfpc_gas /= 0 .AND. Xfpc_surf /= 0) then
          call out_diffsx((Xfpc(4,10)/Xfpc_gas),Cs)

          if (Cs >= fmsurf) then
              call residualvectorset4(1.0d0,unknown*grid+BHin+BHout,FM1)
          else
              FM1=log(FM_foam)+epsurf*log((Xfpc(4,10)/Xfpc_gas)/fmsurf)
              FM1=exp(FM1)
          end if
      else
          call residualvectorset3(unknown*grid+BHin+BHout,FM1)
      end if

      call out_diffsx(FM1,FM1_real)
      write(41,*)"FM :: ",FM1_real
      write(41,*)"CO2+HC@gas :: ",Xfpc_gas
      write(41,*)"SDS@water :: ",Xfpc_surf

      Mg_f=FM1*(fmHC*(Xfpc(1,1)*MW1+Xfpc(1,2)*MW2+Xfpc(1,3)*MW3+Xfpc(1,4)*MW4+Xfpc(1,5)*MW5  &
               +Xfpc(1,6)*MW6+Xfpc(1,7)*MW7)+Xfpc(1,8)*MW8)


      if(FM1_real /= 0.0d0) then
          beta=Mg_f*(1.0d0-Wtg)/(Xfpc_gas*FM1*MW10)*(Wts/Wtg)
          alpha=(1.0d0-Wts)*beta*MW10/(MW9*Wts)
          call out_diffsx(alpha,alpha_real)
          call out_diffsx(beta,beta_real)

          if (Xfpc_H2O > Xfpc_gas*FM1_real*alpha_real) then
              if (Xfpc_surf > Xfpc_gas*FM1_real*beta_real) then
                  write(*,*)"FM→foam"
                  call residualvectorset3(unknown*grid+BHin+BHout,Xfpc_gas_diff)
                  do c=1,7
                      Xfpc_gas_diff=Xfpc_gas_diff+fmHC*Xfpc(1,c)
                      Xfpc(2,c)=FM1*fmHC*Xfpc(1,c)
                      Xfpc(1,c)=(1.0d0-FM1*fmHC)*Xfpc(1,c)
                  end do
                  Xfpc_gas_diff=Xfpc_gas_diff+Xfpc(1,8)
                  Xfpc(2,8)=FM1*Xfpc(1,8)
                  Xfpc(1,8)=(1.0d0-FM1)*Xfpc(1,8)
                  
                  
                  Xfpc(2,9)=Xfpc_gas_diff*FM1*alpha
                  Xfpc(2,10)=Xfpc_gas_diff*FM1*beta
                  Xfpc(4,9)=Xfpc(4,9)-Xfpc(2,9)
                  Xfpc(4,10)=Xfpc(4,10)-Xfpc(2,10)

                  call out_diffsx(Xfpc(4,9),Xfpc_H2O)
                  call out_diffsx(Xfpc(4,10),Xfpc_surf)

              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SDSをすべてfoamにしてる。
                  write(*,*)"SDS→foam"
                  Xfpc(2,10)=Xfpc(4,10) 
                  call residualvectorset3(unknown*grid+BHin+BHout,Xfpc(4,10))
  
                  do c=1,7
                      Xfpc(2,c)=Xfpc(1,c)/Xfpc_gas_diff*1.0d0/beta*Xfpc(2,10)*fmHC
                      Xfpc(1,c)=Xfpc(1,c)-Xfpc(2,c)
                  end do
  
                  Xfpc(2,8)=Xfpc(1,8)/Xfpc_gas_diff*1.0d0/beta*Xfpc(2,10)
                  Xfpc(1,8)=Xfpc(1,8)-Xfpc(2,8)
                      
                  Xfpc(2,9)=alpha/beta*Xfpc(2,10)
                  Xfpc(4,9)=Xfpc(4,9)-Xfpc(2,9)
              end if
          else
              if (Xfpc_surf > Xfpc_gas*FM1_real*beta_real) then
                  write(*,*)"water→foam"
                  Xfpc(2,9)=Xfpc(4,9) 
                  call residualvectorset3(unknown*grid+BHin+BHout,Xfpc(4,9))

                  do c=1,7
                      Xfpc(2,c)=Xfpc(1,c)/Xfpc_gas_diff*1.0d0/alpha*Xfpc(2,9)*fmHC
                      Xfpc(1,c)=Xfpc(1,c)-Xfpc(2,c)
                  end do
                  
                  Xfpc(2,8)=Xfpc(1,8)/Xfpc_gas_diff*1.0d0/alpha*Xfpc(2,9)
                  Xfpc(1,8)=Xfpc(1,8)-Xfpc(2,8)
               
                                    

                  Xfpc(2,10)=beta/alpha*Xfpc(2,9)
                  Xfpc(4,10)=Xfpc(4,10)-Xfpc(2,10)
              else
                  if (1.0d0/alpha_real*Xfpc_H2O > 1.0d0/beta_real*Xfpc_surf) then
                      write(*,*)"SDS→foam"
                      Xfpc(2,10)=Xfpc(4,10) 
                      call residualvectorset3(unknown*grid+BHin+BHout,Xfpc(4,10))
    
                      do c=1,7
                          Xfpc(2,c)=Xfpc(1,c)/Xfpc_gas_diff*1.0d0/beta*Xfpc(2,10)*fmHC
                          Xfpc(1,c)=Xfpc(1,c)-Xfpc(2,c)
                      end do
    
                      Xfpc(2,8)=Xfpc(1,8)/Xfpc_gas_diff*1.0d0/beta*Xfpc(2,10)
                      Xfpc(1,8)=Xfpc(1,8)-Xfpc(2,8)
                          
                      
                      Xfpc(2,9)=alpha/beta*Xfpc(2,10)
                      Xfpc(4,9)=Xfpc(4,9)-Xfpc(2,9)
                  else
                    write(*,*)"water→foam"
                    Xfpc(2,9)=Xfpc(4,9) 
                    call residualvectorset3(unknown*grid+BHin+BHout,Xfpc(4,9))
  
                    do c=1,7
                        Xfpc(2,c)=Xfpc(1,c)/Xfpc_gas_diff*1.0d0/alpha*Xfpc(2,9)*fmHC
                        Xfpc(1,c)=Xfpc(1,c)-Xfpc(2,c)
                    end do
                    
                    Xfpc(2,8)=Xfpc(1,8)/Xfpc_gas_diff*1.0d0/alpha*Xfpc(2,8)
                    Xfpc(1,8)=Xfpc(1,8)-Xfpc(2,8)
                        
  
                    Xfpc(2,10)=beta/alpha*Xfpc(2,9)
                    Xfpc(4,10)=Xfpc(4,10)-Xfpc(2,10) 
                  end if
              end if
          end if
      else
          call residualvectorset3(unknown*grid+BHin+BHout,alpha)
          call residualvectorset3(unknown*grid+BHin+BHout,beta)
      end if

      call residualvectorset3(unknown*grid+BHin+BHout,V)
      call residualvectorset3(unknown*grid+BHin+BHout,F)
      call residualvectorset3(unknown*grid+BHin+BHout,L)
      call residualvectorset3(unknown*grid+BHin+BHout,W)

      do com=1,t
          V=V+Xfpc(1,com)
          F=F+Xfpc(2,com)
          L=L+Xfpc(3,com)
          W=W+Xfpc(4,com)
      end do

      call out_diffsx(V,V0)
      call out_diffsx(F,F0)
      call out_diffsx(L,L0)
      call out_diffsx(W,W0)
      
      
      if(V0 /= 0.0d0)then
          do com=1,t
              Xpc(1,com)=Xfpc(1,com)/V
          end do
      else
          call residualvectorset3(unknown*grid+BHin+BHout,V)
          do com=1,t
              call residualvectorset3(unknown*grid+BHin+BHout,Xpc(2,com))
          end do
      end if
              
      do com=1,t
        Xpc(1,com)=Xfpc(1,com)/V
        Xpc(3,com)=Xfpc(3,com)/L
        Xpc(4,com)=Xfpc(4,com)/W

          if (FM1_real /= 0) then
              Xpc(2,com)=Xfpc(2,com)/F
          else 
              call residualvectorset3(unknown*grid+BHin+BHout,Xpc(2,com))
              call residualvectorset3(unknown*grid+BHin+BHout,F)
          end if
      end do
 
    end subroutine calc_foam_formation
end module