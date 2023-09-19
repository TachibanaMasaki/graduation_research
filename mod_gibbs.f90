module mod_gibbs
  use mod_initial
  use mod_autodiff
  use mod_interaction
  use mod_cubic_eq
  contains
  subroutine calc_gibbs(alpha0,Ptotal,X_exist,like,fxs)
    implicit none
    real(8),intent(in),dimension(n)::alpha0
    real(8),intent(in),dimension(n)::X_exist
    real(8),intent(in)::Ptotal
    integer,intent(in)::like
    type(diffs),allocatable,intent(out)::fxs(:)
    real(8),allocatable,dimension(:)::x0,lnf_exist,aaa
    integer::i,c,cc
    type(diffs),allocatable,target::xd(:)
    type(diffs),pointer::alpha(:)
    type(diffs)::WW(n),X_appear(n),sumWW,lnf_appear(n)
        
    allocate(x0(n),lnf_exist(n),aaa(n))

    
    do c=1,n
        x0(c)=alpha0(c)
    end do
        
    call diffsset1(x0,xd)
    call sizeset(x0,fxs)
    
    alpha => xd
    
    do c=1,n
        WW(c)=(alpha(c)/2.0d0)**2.0d0
    end do
        
    
    call residualvectorset3(n,sumWW)
    do c=1,n
        sumWW=sumWW+WW(c)
    end do
    do c=1,n
        X_appear(c)=WW(c)/sumWW
    end do
    
    call calc_fugacity_appear(Ptotal,X_appear,like,lnf_appear)
    call calc_fugacity_exist(Ptotal,X_exist,like,lnf_exist)

    do c=1,n
        fxs(c)=sqrt(WW(c))*(log(WW(c))+lnf_appear(c)-(log(X_exist(c))+lnf_exist(c)))
    end do

  end subroutine calc_gibbs
  
  subroutine calc_fugacity_appear(Ptotal,X,like,lnf)
    implicit none
    real(8),intent(in)::Ptotal
    integer,intent(in)::like
    type(diffs),intent(in),dimension(n)::X
    type(diffs),intent(out),dimension(n)::lnf
    real(8),allocatable,dimension(:,:)::Pr,Theta
    real(8),allocatable,dimension(:)::Tb,Pp,Pc,Tc,Vc,omega,Tr,ac,bc,acc,bcc,alf,Aco,Bco,Cco,zfactor
    real(8)::R
    integer::p,c,cc    
    type(diffs)::ap(m_1),bp(m_1),LA(m_1),LB(m_1)  &
                  ,A_coef(m_1),B_coef(m_1),C_coef(m_1),zfac(m_1),middle(m_1,n)

    allocate(Tb(n),Pp(m_1),Pr(m_1,n),Pc(n),Tc(n),Vc(n),omega(n),Tr(n),ac(n),bc(n)  &
              ,acc(n),bcc(n),alf(n),Theta(n,n),Aco(m_1),Bco(m_1),Cco(m_1),zfactor(m_1))
    
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
      Pp=Ptotal  
   

!気体定数 (Pa*m^3/K)
      R=8.3144598d0 !(Pa*m^3/K) 

! (p,c)=(phase,component)=(m,n)

! Pr(p,c) = p相におけるc成分の還元圧力 (-)
      do p=1,m_1
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

        do p=1,m_1
            call residualvectorset3(n,ap(p))
            call residualvectorset3(n,bp(p))
        end do

        do p=1,m_1
            do c=1,n
                do cc=1,n
                    ap(p)=ap(p)+X(n*(p-1)+c)*X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))               
                end do
            end do
        end do
    
        do p=1,m_1
            do c=1,n
              bp(p)=bp(p)+X(n*(p-1)+c)*bc(c)    
            end do
        end do

      do p=1,m_1
          LA(p)=ap(p)*Pp(p)/(R*T0)**2.0d0    
          LB(p)=bp(p)*Pp(p)/(R*T0)           
      end do

      do p=1,m_1
          A_coef(p)=-1.0d0+LB(p)                             
          B_coef(p)=LA(p)-3.0d0*LB(p)**2.0d0-2.0d0*LB(p)      
          C_coef(p)=LB(p)**3.0d0-LA(p)*LB(p)+LB(p)**2.0d0    
      end do

      call outxs(A_coef,Aco)
      call outxs(B_coef,Bco)
      call outxs(C_coef,Cco)

      if (like == 1)then
          !call cubic_eq_other(Aco(1),Bco(1),Cco(1),2,zfactor(1))
          call cubic_eq(Aco(1),Bco(1),Cco(1),0.5d0,zfactor(1))
          write(*,*)"1",zfactor
          !call cubic_eq(Aco(1),Bco(1),Cco(1),1,zfactor(1))
          !call cubic_eq(Aco(1),Bco(1),Cco(1),1.0d0,zfactor(1))
          !write(*,*)"2",zfactor
          
      else
          call cubic_eq_other(Aco(1),Bco(1),Cco(1),1,zfactor(1))
      end if
      
      !write(*,*)"zfactor",zfactor

      do p=1,m_1
          call residualvectorset4(zfactor(p),n,zfac(p))
      end do


!フガシティ計算 (気相、油相、水相)
      do p=1,m_1
            do c=1,n
                call residualvectorset3(n,middle(p,c))
            end do
      end do
        do p=1,m_1
            do c=1,n
                do cc=1,n
                    middle(p,c)=middle(p,c)+X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))  
                end do
            end do
        end do

    do p=1,m_1
        do c=1,n
            lnf(c)=bc(c)/bp(p)*(zfac(p)-1.0d0)-log(zfac(p)-LB(p))  &
                    -LA(p)/(2.0d0*dsqrt(2.0d0)*LB(p))*(2.0d0*middle(p,c)/ap(p)-bc(c)/bp(p))  &
                    *log((zfac(p)+LB(p)*(1.0d0+dsqrt(2.0d0)))/(zfac(p)+LB(p)*(1.0d0-dsqrt(2.0d0))))
        end do
    end do
    
  end subroutine calc_fugacity_appear
  
  
  
  subroutine calc_fugacity_exist(Ptotal,X,like,lnf)
    implicit none
    real(8),intent(in)::Ptotal
    integer,intent(in)::like
    real(8),intent(in),dimension(n)::X
    real(8),intent(out),dimension(n)::lnf
    integer::p,c,cc 
    real(8)::R
    real(8),allocatable,dimension(:,:)::Pr,Theta,middle
    real(8),allocatable,dimension(:)::Pp,ap,bp,LA,LB,A_coef,B_coef,C_coef,zfac,Tb,Pc,Tc,Vc,omega,Tr,ac,bc,acc,bcc,alf  &
                                      ,Aco,Bco,Cco,zfactor

    allocate(Tb(n),Pc(n),Pp(m_1),Pr(m_1,n),Tc(n),Vc(n),omega(n),Tr(n),ac(n),bc(n),  &
            ap(m_1),bp(m_1),LA(m_1),LB(m_1),A_coef(m_1),B_coef(m_1),C_coef(m_1),zfac(m_1),  &
            acc(n),bcc(n),alf(n),Theta(n,n),Aco(m_1),Bco(m_1),Cco(m_1),zfactor(m_1),middle(m_1,n))

    write(*,*)X

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
      Pp=Ptotal   

!気体定数 (Pa*m^3/K)
      R=8.3144598d0 !(Pa*m^3/K) 

! (p,c)=(phase,component)=(m,n)

! Pr(p,c) = p相におけるc成分の還元圧力 (-)
      do p=1,m_1
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
           
        ap=0.0d0
        bp=0.0d0

        do p=1,m_1
            do c=1,n
                do cc=1,n
                    ap(p)=ap(p)+X(n*(p-1)+c)*X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))               
                end do
            end do
        end do
    
        do p=1,m_1
            do c=1,n
              bp(p)=bp(p)+X(n*(p-1)+c)*bc(c)    
            end do
        end do

      do p=1,m_1
          LA(p)=ap(p)*Pp(p)/(R*T0)**2.0d0    
          LB(p)=bp(p)*Pp(p)/(R*T0)           
      end do

      do p=1,m_1
          A_coef(p)=-1.0d0+LB(p)                             
          B_coef(p)=LA(p)-3.0d0*LB(p)**2.0d0-2.0d0*LB(p)      
          C_coef(p)=LB(p)**3.0d0-LA(p)*LB(p)+LB(p)**2.0d0    
      end do

      write(*,*)A_coef,B_coef,C_coef
      if (like == 1)then
          call cubic_eq_other(A_coef(1),B_coef(1),C_coef(1),1,zfactor(1))
      else
          call cubic_eq_other(A_coef(1),B_coef(1),C_coef(1),2,zfactor(1))
      end if
   
      zfac=zfactor
      write(*,*)"zfactor",zfactor

!フガシティ計算 (気相、油相、水相)
      middle=0.0d0
        do p=1,m_1
            do c=1,n
                do cc=1,n
                    middle(p,c)=middle(p,c)+X(n*(p-1)+cc)*(1.0d0-Theta(c,cc))*dsqrt(ac(c)*ac(cc))  
                end do
            end do
        end do

    do p=1,m_1
        do c=1,n
            lnf(c)=bc(c)/bp(p)*(zfac(p)-1.0d0)-log(zfac(p)-LB(p))  &
                    -LA(p)/(2.0d0*dsqrt(2.0d0)*LB(p))*(2.0d0*middle(p,c)/ap(p)-bc(c)/bp(p))  &
                    *log((zfac(p)+LB(p)*(1.0d0+dsqrt(2.0d0)))/(zfac(p)+LB(p)*(1.0d0-dsqrt(2.0d0))))
        end do
    end do

            
  end subroutine calc_fugacity_exist
  
end module mod_gibbs