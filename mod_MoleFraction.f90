module mod_MoleFraction
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_MoleFraction(unknown,Nmol,lnKgo,lnKwo,V0,L0,W0,V,L,W,Xf0pc)
      implicit none
      real(8),intent(in)::V0,L0,W0
      type(diffs),intent(in),dimension(t)::Nmol
      type(diffs),intent(in),dimension(n)::lnKgo,lnKwo
      integer,intent(in)::unknown
      type(diffs),intent(out),dimension(m,t)::Xf0pc
      type(diffs),intent(out)::V,L,W
      integer::com,p,c,i
      real(8)::dammy
      type(diffs)::ntotal0,ntotal,lnK((m-1)*n),X0(m*n),Xf0(m*n),N0pc(m,t)

      call residualvectorset3(unknown*grid+BHin+BHout,ntotal0)
      do c=1,n
        ntotal0=ntotal0+Nmol(c)
      end do

      call residualvectorset3(unknown*grid+BHin+BHout,ntotal)
      do com=1,t
        ntotal=ntotal+Nmol(com)
      end do

      do c=1,n 
        lnK(c)=lnKgo(c)
        lnK(c+n)=lnKwo(c)
      end do

      do c=1,n 
          X0(c)=exp(lnK(c))*(Nmol(c)/ntotal0)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
          X0(c+n)=(Nmol(c)/ntotal0)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
          X0(c+2*n)=exp(lnK(c+n))*(Nmol(c)/ntotal0)/(V0*exp(lnK(c))+L0+W0*exp(lnK(c+n)))
      end do
            

    do c=1,n
        Xf0(c)=X0(c)*V0
        Xf0(c+n)=X0(c+n)*L0
        Xf0(c+2*n)=X0(c+2*n)*W0
    end do

    do p=1,m
        do com=1,t
            call residualvectorset3(unknown*grid+BHin+BHout,N0pc(p,com))
        end do
    end do

    do p=1,m
        do c=1,n
            N0pc(p,c)=Xf0(c+(p-1)*n)*ntotal0
        end do
    end do

!平衡状態にない物質はすべて水相にあるとしている。
    do i=1,t-n
        N0pc(3,n+i)=Nmol(n+i)
    end do

    do p=1,m
        do com=1,t
            Xf0pc(p,com)=N0pc(p,com)/ntotal
        end do
    end do

    do p=1,m
        do com=1,t
            call out_diffsx(Xf0pc(p,com),dammy)
            if(dammy < 0.0d0)then
                call residualvectorset3(unknown*grid+BHin+BHout,Xf0pc(p,com))
                write(41,*)"Xf0pc(p,com)<0",p,com,dammy
            end if
        end do
    end do
    
    call residualvectorset3(unknown*grid+BHin+BHout,V)
    call residualvectorset3(unknown*grid+BHin+BHout,L)
    call residualvectorset3(unknown*grid+BHin+BHout,W)

    do com=1,t
        V=V+Xf0pc(1,com)
        L=L+Xf0pc(2,com)
        W=W+Xf0pc(3,com)
    end do

end subroutine calc_MoleFraction
end module