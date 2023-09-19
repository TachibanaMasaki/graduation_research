module mod_initial_Value
    use mod_initial
    use mod_mole
    contains 
    subroutine calc_initial_Value(Nc0,Nc0old,Kgo0,Kwo0,V0,F0,L0,W0,Pbh0in,Pbh0out,Po0,Poold,faiold,Vgsc_ini,Vosc_ini,K0)
        implicit none
        real(8),intent(out),dimension(t,grid)::Nc0,Nc0old
        real(8),intent(out),dimension(n,grid)::Kgo0,Kwo0
        real(8),intent(out),dimension(grid)::V0,F0,L0,W0,Po0,Poold,faiold
        real(8),intent(out),dimension((m-1)*n)::K0
        real(8),intent(out)::Pbh0in,Pbh0out,Vgsc_ini,Vosc_ini
        integer::com,c,i
        real(8),allocatable,dimension(:)::Zcin,Ncin,Kgoin,Kwoin,Vp_SC,Vp
        real(8)::Press,Temp,V_SC,F_SC,L_SC,W_SC,V,F,L,W,Sg,Sf,So,Sw,Vg_ini,Vo_ini,ntotal
        allocate(Zcin(t),Ncin(t),Kgoin(n),Kwoin(n),Vp_SC(m),Vp(m))
 
        do com=1,t
            read(1,*)Zcin(com)
        end do 

        do c=1,n
            read(1,*)Kgoin(c)
        end do

        do c=1,n
            read(1,*)Kwoin(c)
        end do

        Press=Psc
        Temp=Tsc
        call calc_mole(Press,Temp,Zcin,Kgoin,Kwoin,V_SC,F_SC,L_SC,W_SC,Sg,Sf,So,Sw,Vp_SC,Ncin)


        Press=Pinitial
        Temp=T0
        call calc_mole(Press,Temp,Zcin,Kgoin,Kwoin,V,F,L,W,Sg,Sf,So,Sw,Vp,Ncin)
        write(*,*)Ncin
!!!!! NcinÇÕ1 m^3Ç†ÇΩÇËÇÃÇªÇÍÇºÇÍÇÃÉÇÉãêî!!!!!!!!!!!!!!
        ntotal=0.0d0
        do com=1,t
            ntotal=ntotal+Ncin(com)
        end do
        ntotal=ntotal*A*Long*poro
        
        Vgsc_ini=ntotal*V_SC*Vp_SC(1)
        Vosc_ini=ntotal*L_SC*Vp_SC(2)
        
        do com=1,t
            do i=1,grid
                Nc0(com,i)=Ncin(com)
            end do
        end do

        do com=1,t
            do i=1,grid
                Nc0old(com,i)=Nc0(com,i)
            end do
        end do

        do c=1,n
            do i=1,grid
                Kgo0(c,i)=Kgoin(c)
            end do
        end do
        
        do c=1,n
            do i=1,grid
                Kwo0(c,i)=Kwoin(c)
            end do
        end do

        do i=1,grid
            V0(i)=V
            F0(i)=F
            L0(i)=L
            W0(i)=W
        enddo


        write(4,'(f10.7,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',Sg,',',Sg,',',Sg
        write(5,'(f10.7,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',Sf,',',Sf,',',Sf
        write(6,'(f10.7,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',So,',',So,',',So
        write(7,'(f10.7,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',Sw,',',Sw,',',Sw

        
        do i=1,grid
            write(8,*)Long/grid*(i-0.5),',',Sg
            write(9,*)Long/grid*(i-0.5),',',Sf
            write(10,*)Long/grid*(i-0.5),',',So
            write(11,*)Long/grid*(i-0.5),',',Sw
        end do
        

        Po0=Pinitial
        Poold=Pinitial
        faiold=poro
        Pbh0in=P_bhin
        Pbh0out=P_bhout

        do c=1,n
            K0(c)=Kgoin(c)
            K0(c+n)=Kwoin(c)
        end do
    end subroutine calc_initial_Value
end module 