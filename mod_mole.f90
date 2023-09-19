module mod_mole
    use mod_initial
    use mod_autodiff
    use mod_flash
    use mod_Kvalue
    use mod_pvgauss
    use mod_MolarVol
    contains 
    subroutine calc_mole(Press,Temp,Zcin,Kgo,Kwo,V0,F0,L0,W0,Sg,Sf,So,Sw,Vp,Nc)
    implicit none
    real(8),intent(in)::Press,Temp
    real(8),intent(in),dimension(t)::Zcin
    real(8),intent(inout),dimension(n)::Kgo,Kwo
    real(8),intent(out),dimension(t)::Nc
    real(8),intent(out),dimension(m)::Vp
    real(8),intent(out)::V0,F0,L0,W0,Sg,Sf,So,Sw
    integer::com,c,r,k,conv,i,j,roopmax_ini=100
    real(8),allocatable,dimension(:)::z0,K0,lnK0,residual1,residual2
    real(8)::ntotal,dL0,dW0,error,epsilon_ini=10.0d0**(-10.0d0)
    real(8),allocatable,dimension(:,:)::jacobian1,jacobian2
    type(diffs),allocatable::fxs1(:),fxs2(:)

    allocate(z0(n),K0((m-1)*n),jacobian1(m-1,m-1),residual1(m-1)  &
            ,jacobian2((m-1)*n,(m-1)*n),residual2((m-1)*n),lnK0((m-1)*n))

    V0=0.5d0
    L0=0.4d0
    W0=0.1d0

!成分のモル分率の入力

    ntotal=0.0d0
    do c=1,n
        ntotal=ntotal+Zcin(c)
    end do
    do c=1,n
        z0(c)=Zcin(c)/ntotal
    end do

    do c=1,n
        K0(c)=Kgo(c)
        K0(c+n)=Kwo(c)
    end do

    do r=1,roopmax_ini
        do k=1,roopmax
            call calc_flash(z0,L0,W0,Kgo,Kwo,fxs1)
            call jacobian_mat(fxs1,jacobian1)
            call outxs(fxs1,residual1)
            residual1=-residual1
            call pvgauss(m-1,jacobian1,residual1)

            dL0=residual1(1)
            dW0=residual1(2)
    
            L0=L0+dL0
            W0=W0+dW0
            V0=1.0d0-L0-W0
            F0=0.0d0
            error=dot_product(residual1,residual1)
            error=sqrt(error)
            if(error<epsilon_ini)then
                conv=k
                exit  
            end if
        end do

    
        call calc_Kvalue(Press,Temp,z0,V0,L0,W0,Kgo,Kwo,fxs2)
        call jacobian_mat(fxs2,jacobian2)
        call outxs(fxs2,residual2)

        do i=1,(m-1)*n
          if(K0(i)==0.0d0)then
            do j=1,(m-1)*n
              jacobian2(i,j)=0.0d0
              jacobian2(j,i)=0.0d0
            end do
            jacobian2(i,i)=1.0d0
            residual2(i)=0.0d0
          end if
        end do
 
        residual2=-residual2


        call pvgauss((m-1)*n,jacobian2,residual2)
        
        do c=1,n 
            lnK0(c)=log(Kgo(c))+residual2(c)
            lnK0(c+n)=log(Kwo(c))+residual2(c+n)
        end do

        do c=1,n 
            Kgo(c)=exp(lnK0(c))
            Kwo(c)=exp(lnK0(c+n))
        end do

        error=dot_product(residual2,residual2)
        error=sqrt(error)

        write(*,*)error
        write(*,*)Kgo,Kwo

        if(error<epsilon_ini)then
            conv=r
            exit  
        end if
    end do

    call calc_MolarVol(Press,Temp,Zcin,V0,L0,W0,Kgo,Kwo,Vp)


    write(*,*)"Vp",Vp
    !Ntotal 1 m^3あたりのmol数
    Ntotal=1.0d0/(Vp(1)*V0+Vp(2)*L0+Vp(3)*W0)

    do com=1,t
        Nc(com)=Zcin(com)*ntotal
    end do
    
    
    Sg=Vp(1)*V0/(Vp(1)*V0+Vp(2)*L0+Vp(3)*W0)
    Sf=0.0d0
    So=Vp(2)*L0/(Vp(1)*V0+Vp(2)*L0+Vp(3)*W0)
    Sw=Vp(3)*W0/(Vp(1)*V0+Vp(2)*L0+Vp(3)*W0)
   
            
end subroutine calc_mole
end module 