module mod_Phase_stability
    use mod_initial
    use mod_flash_2p
    use mod_Liquid_vapor_like
    use mod_pvgauss
    contains
    subroutine calc_Phase_stability(Po0,Ncn0,Kgo0,Kwo0,V0,L0,W0,Pnum)
    implicit none
    real(8),intent(in),dimension(n,grid)::Ncn0,Kgo0,Kwo0
    real(8),intent(in),dimension(grid)::Po0,V0,L0,W0
    integer,intent(out),dimension(grid)::Pnum
    real(8),allocatable,dimension(:,:)::GZ,Xoil,KKgo0
    real(8),allocatable,dimension(:)::ntotal,VV0,LL0,WW0
    integer,allocatable,dimension(:)::like
    integer::i,c

    allocate(ntotal(grid),GZ(n,grid),Xoil(n,grid),KKgo0(n,grid),like(grid),VV0(grid),LL0(grid),WW0(grid))
    
    VV0=V0
    LL0=L0
    WW0=W0
    
    ntotal=0.0d0
    do i=1,grid
        do c=1,n
            ntotal(i)=ntotal(i)+Ncn0(c,i)
        end do
    end do

    do i=1,grid
        do c=1,n
            GZ(c,i)=Ncn0(c,i)/ntotal(i)
        end do
    end do
    

    
        do i=1,grid
        if(VV0(i)>=LL0(i))then
            like(i)=1
        else
            like(i)=2
        end if
    end do

    !気相を0としたときの2相フラッシュ
    call calc_flash_2p(GZ,Kwo0,VV0,LL0,WW0)

!2相で考える。1相:油相  2相:出現した相(気相)
    do i=1,grid
        do c=1,n
            if(Kwo0(c,i)==0.0d0)then
                Xoil(c,i)=GZ(c,i)/LL0(i)
            else
                Xoil(c,i)=GZ(c,i)/(LL0(i)+WW0(i)*Kwo0(c,i))
            end if
        end do
    end do
 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　　　　確かめ用　　　　!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !do i=1,grid
    !    Xoil(1,i)=0.45d0
    !    Xoil(2,i)=0.0d0
    !    Xoil(3,i)=0.0d0
    !    Xoil(4,i)=0.0d0
    !    Xoil(5,i)=0.0d0
    !    Xoil(6,i)=0.0d0
    !    Xoil(7,i)=0.01d0
    !    Xoil(8,i)=0.29d0
    !    Xoil(9,i)=0.25d0
    !    !自分とwinprop一致
    !    !Xoil(7,i)を0.1→0.0
    !    !Xoil(8,i)を0.29→0.3
    !    !!にするとWVの2phase 
    !end do

    !like=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    KKgo0=Kgo0
    call calc_Liquid_vapor_like(Po0,Xoil,KKgo0,like,Pnum)
    
!!いったん消した後で使うかも
    
    !do i=1,grid
    !    do c=1,n
    !        if(Xoil(c,i)<=0.0d0)then
    !            Kgo0(c,i)=0.0d0
    !        end if
    !    end do
    !end do

     
    do i=1,grid
        write(*,*)"phase number",Pnum(i)
    end do
    
    end subroutine calc_Phase_stability
end module