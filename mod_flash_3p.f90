module mod_flash_3p
    use mod_initial
    use mod_autodiff
    implicit none
    contains 
    subroutine calc_flash_3p(n0,L0,W0,Kgo,Kwo,g)
      implicit none
      real(8),intent(in)::L0,W0
      real(8),intent(in),dimension(n)::n0,Kgo,Kwo
      type(diffs),allocatable,intent(out)::g(:)
      real(8),allocatable,dimension(:)::x0,z0
      real(8)::ntotal
      type(diffs),allocatable,target::xd(:)
      type(diffs),pointer::L,W
      type(diffs)::V,X(m*n),gg(m-1)
      integer::p,c

      allocate(x0(m-1),z0(n))

      x0(1)=L0
      x0(2)=W0
  
      call diffsset1(x0,xd)
      call sizeset(x0,g)
  
      L   => xd(1)
      W   => xd(2)
    

      V=1.0d0-L-W

      ntotal=0.0d0
      do c=1,n
          ntotal=ntotal+n0(c)
      end do

      do c=1,n
          z0(c)=n0(c)/ntotal
      end do
    
      do c=1,n
          X(c)=Kgo(c)*z0(c)/(V*Kgo(c)+L+W*Kwo(c))
          X(c+n)=z0(c)/(V*Kgo(c)+L+W*Kwo(c))
          X(c+2*n)=Kwo(c)*z0(c)/(V*Kgo(c)+L+W*Kwo(c))
      end do

      do p=1,m-1
          call residualvectorset3(m-1,gg(p))
      end do

      do p=1,m-1
          do c=1,n
              gg(p)=gg(p)+X(c+n*(p-1))-X(c+n*p)
          end do
      end do
    
      g(1)=gg(1)
      g(2)=gg(2)

    end subroutine calc_flash_3p
end module mod_flash_3p