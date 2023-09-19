module mod_main
    use mod_initial
    use mod_autodiff
    use mod_Phase_stability
    use mod_Rachford_Rice
    use mod_Thermo_EQ
    use mod_MoleFraction
    use mod_foam_formation
    use mod_density
    use mod_permeability
    use mod_permeability_interpolation
    use mod_viscosity
    use mod_flow
    use mod_Saturation
    use mod_flow_control
    use mod_ORR
    implicit none
    contains 
    subroutine calc_main(time_roop,r,Nup,V0,L0,W0,Sg0,Sf0,So0,Sw0,Kgo0,Kwo0,Nc0,Nc0old,Po0,Pbh0in,Pbh0out,faid,Poold,Vgsc,Vosc,fxs)
        implicit none
        integer,intent(in)::time_roop,r
        integer,intent(out),dimension(t,grid)::Nup
        real(8),intent(in),dimension(n,grid)::Kgo0,Kwo0
        real(8),intent(in),dimension(t,grid)::Nc0,Nc0old
        real(8),intent(in)::Pbh0in,Pbh0out
        real(8),intent(inout),dimension(grid)::V0,L0,W0,Po0,faid,Poold
        real(8),intent(out),dimension(grid)::Sg0,Sf0,So0,Sw0
        type(diffs),allocatable,intent(out)::fxs(:)
        integer,allocatable,dimension(:)::Pnum
        real(8),allocatable,dimension(:,:)::lnKgo0,lnKwo0,Ncn0,S0_p
        real(8),allocatable,dimension(:)::x0,Ncn0_i,Kgo0_i,Kwo0_i
        real(8),allocatable,dimension(:)::mu_real,rho_real,kr_real
        integer::unknown,i,p,c,com,k,phase
        real(8)::F0_real,ntotal0,Vgsc,Vosc,dammy
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::lnKgo1(:),lnKgo2(:),lnKgo3(:),lnKgo4(:),lnKgo5(:),lnKgo6(:),lnKgo7(:),lnKgo8(:),lnKgo9(:)  &
                            ,lnKwo1(:),lnKwo2(:),lnKwo3(:),lnKwo4(:),lnKwo5(:),lnKwo6(:),lnKwo7(:),lnKwo8(:),lnKwo9(:)  &
                            ,Nc1(:),Nc2(:),Nc3(:),Nc4(:),Nc5(:),Nc6(:),Nc7(:),Nc8(:),Nc9(:),Nc10(:),Nc11(:),Po(:),Pbhin,Pbhout
        type(diffs)::g_thermo((m-1)*n),g_f(t),g_s
        type(diffs),allocatable,dimension(:)::g_fc
        type(diffs)::lnKgo(n,grid),lnKwo(n,grid),Nc(t,grid),Ncn(n,grid),Nc_i(t),Ncn_i(n),lnKgo_i(n),lnKwo_i(n)  &
                    ,V(grid),F(grid),L(grid),W(grid),Xf0pc(m,t),Xfpc(s,t),Xpc(s,t),FM1(grid)  &
                    ,rho_i(s),Vp_i(s),rho(s,grid),Vp(s,grid),krg(grid),krf(grid),kro(grid),krw(grid),kr_p(s,grid)  &
                    ,mu_i(s),mu(s,grid),Sg(grid),Sf(grid),So(grid),Sw(grid)  &
                    ,Cgc(t,grid),Cfc(t,grid),Coc(t,grid),Cwc(t,grid),qin(s),qout(s),qin0(s),qout0(s),ntotal

        unknown=(m-1)*n+t+1
        
        Nup=0

        allocate(Pnum(grid))
        allocate(x0(unknown*grid+BHin+BHout),lnKgo0(n,grid),lnKwo0(n,grid),Ncn0(n,grid),Ncn0_i(n),Kgo0_i(n),Kwo0_i(n))
        allocate(S0_p(s,grid))
        allocate(mu_real(s),rho_real(s),kr_real(s))

        lnKgo0=log(Kgo0)
        lnKwo0=log(Kwo0)

        do c=1,n
            do i=1,grid
                x0(unknown*(i-1)+c)=lnKgo0(c,i)
                x0(unknown*(i-1)+c+n)=lnKwo0(c,i)
            end do
        end do
       
        do com=1,t
            do i=1,grid
                x0(unknown*(i-1)+com+2*n)=Nc0(com,i)
            end do
        end do
        do i=1,grid
            x0(unknown*i)=Po0(i)
        end do
        
        if (BHin == 1) then
            x0(unknown*grid+1)=Pbh0in
        end if
    
        if (BHin == 1 .AND. BHout == 1) then
            x0(unknown*grid+2)=Pbh0out
        else if (BHin /= 1 .AND. BHout == 1) then
            x0(unknown*grid+1)=Pbh0out
        end if

        call diffsset1(x0,xd)
        call sizeset(x0,fxs)
        
        lnKgo1 => xd(1:unknown*(grid-1)+1:unknown)
        lnKgo2 => xd(2:unknown*(grid-1)+2:unknown)
        lnKgo3 => xd(3:unknown*(grid-1)+3:unknown)
        lnKgo4 => xd(4:unknown*(grid-1)+4:unknown)
        lnKgo5 => xd(5:unknown*(grid-1)+5:unknown)
        lnKgo6 => xd(6:unknown*(grid-1)+6:unknown)
        lnKgo7 => xd(7:unknown*(grid-1)+7:unknown)
        lnKgo8 => xd(8:unknown*(grid-1)+8:unknown)
        lnKgo9 => xd(9:unknown*(grid-1)+9:unknown)

        lnKwo1 => xd(10:unknown*(grid-1)+10:unknown)
        lnKwo2 => xd(11:unknown*(grid-1)+11:unknown)
        lnKwo3 => xd(12:unknown*(grid-1)+12:unknown)
        lnKwo4 => xd(13:unknown*(grid-1)+13:unknown)
        lnKwo5 => xd(14:unknown*(grid-1)+14:unknown)
        lnKwo6 => xd(15:unknown*(grid-1)+15:unknown)
        lnKwo7 => xd(16:unknown*(grid-1)+16:unknown)
        lnKwo8 => xd(17:unknown*(grid-1)+17:unknown)
        lnKwo9 => xd(18:unknown*(grid-1)+18:unknown)

        Nc1    => xd(19:unknown*(grid-1)+19:unknown)
        Nc2    => xd(20:unknown*(grid-1)+20:unknown)
        Nc3    => xd(21:unknown*(grid-1)+21:unknown)
        Nc4    => xd(22:unknown*(grid-1)+22:unknown)
        Nc5    => xd(23:unknown*(grid-1)+23:unknown)
        Nc6    => xd(24:unknown*(grid-1)+24:unknown)
        Nc7    => xd(25:unknown*(grid-1)+25:unknown)
        Nc8    => xd(26:unknown*(grid-1)+26:unknown)
        Nc9    => xd(27:unknown*(grid-1)+27:unknown)
        Nc10   => xd(28:unknown*(grid-1)+28:unknown)
        Nc11   => xd(29:unknown*(grid-1)+29:unknown)

        Po     => xd(30:unknown*(grid-1)+30:unknown)

        if (BHin == 1) then
            Pbhin  => xd(unknown*grid+1)
        end if

        if (BHin == 1 .AND. BHout == 1) then
            Pbhout  => xd(unknown*grid+2)
        else if (BHin /= 1 .AND. BHout == 1) then
            Pbhout  => xd(unknown*grid+1)
        end if
               
        do i=1,grid
            lnKgo(1,i)=lnKgo1(i)
            lnKgo(2,i)=lnKgo2(i)
            lnKgo(3,i)=lnKgo3(i)
            lnKgo(4,i)=lnKgo4(i)
            lnKgo(5,i)=lnKgo5(i)
            lnKgo(6,i)=lnKgo6(i)
            lnKgo(7,i)=lnKgo7(i)
            lnKgo(8,i)=lnKgo8(i)
            lnKgo(9,i)=lnKgo9(i)

            lnKwo(1,i)=lnKwo1(i)
            lnKwo(2,i)=lnKwo2(i)
            lnKwo(3,i)=lnKwo3(i)
            lnKwo(4,i)=lnKwo4(i)
            lnKwo(5,i)=lnKwo5(i)
            lnKwo(6,i)=lnKwo6(i)
            lnKwo(7,i)=lnKwo7(i)
            lnKwo(8,i)=lnKwo8(i)
            lnKwo(9,i)=lnKwo9(i)
        
            Nc(1,i)=Nc1(i)
            Nc(2,i)=Nc2(i)
            Nc(3,i)=Nc3(i)
            Nc(4,i)=Nc4(i)
            Nc(5,i)=Nc5(i)
            Nc(6,i)=Nc6(i)
            Nc(7,i)=Nc7(i)
            Nc(8,i)=Nc8(i)
            Nc(9,i)=Nc9(i)
            Nc(10,i)=Nc10(i)
            Nc(11,i)=Nc11(i)
            
            Ncn(1,i)=Nc1(i)
            Ncn(2,i)=Nc2(i)
            Ncn(3,i)=Nc3(i)
            Ncn(4,i)=Nc4(i)
            Ncn(5,i)=Nc5(i)
            Ncn(6,i)=Nc6(i)
            Ncn(7,i)=Nc7(i)
            Ncn(8,i)=Nc8(i)
            Ncn(9,i)=Nc9(i)
        end do

        
        do i=1,grid
            do c=1,n
                Ncn0(c,i)=Nc0(c,i)
            end do
        end do
        
!相安定性解析
        call calc_Phase_stability(Po0,Ncn0,Kgo0,Kwo0,V0,L0,W0,Pnum)
        if(r == roopmax)then
            write(42,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Pnum(1)
        end if
        
!?grid iにおけるPVTとPhaseProperty------------------------------------------------------------
        do i=1,grid

            write(*,*)'grid',i
            do c=1,n
                Ncn0_i(c)=Ncn0(c,i)
                Kgo0_i(c)=Kgo0(c,i)
                Kwo0_i(c)=Kwo0(c,i)
            end do
            call calc_Rachford_Rice(Ncn0_i,Kgo0_i,Kwo0_i,V0(i),L0(i),W0(i))

            do c=1,n
                Ncn_i(c)=Ncn(c,i)
                lnKgo_i(c)=lnKgo(c,i)
                lnKwo_i(c)=lnKwo(c,i)
            end do
    
            if(V0(i) < 0.0d0 .OR. L0(i) < 0.0d0 .OR. W0(i) < 0.0d0)then
                write(41,*)"mole fracion",V0(i),L0(i),W0(i)
                write(41,*)"phase",Pnum(1)
            end if
!!Thermodynamic equilibrium
            call calc_Thermo_EQ(unknown,Po(i),Ncn_i,V0(i),L0(i),W0(i),lnKgo_i,lnKwo_i,g_thermo)
            
            do k=1,(m-1)*n
                fxs(k+(i-1)*unknown)=g_thermo(k)
            end do
    
            do com=1,t
                Nc_i(com)=Nc(com,i)
            end do

!!PVT
            call calc_MoleFraction(unknown,Nc_i,lnKgo_i,lnKwo_i,V0(i),L0(i),W0(i),V(i),L(i),W(i),Xf0pc)
            
            
            do p=1,m
                do com=1,t
                    call out_diffsx(Xf0pc(p,com),dammy)
                    if(dammy < 0.0d0)then
                        write(41,*)
                        write(41,*)"Xf0pc < 0.0d0",p,com,dammy
                        write(41,*)
                    end if
                end do
            end do
                    
    
            
            call calc_foam_formation(unknown,Xf0pc,Xfpc,Xpc,V(i),F(i),L(i),W(i),FM1(i))
            do phase=1,s
                do com=1,t
                    call out_diffsx(Xpc(phase,com),dammy)
                    if(dammy < 0.0d0)then
                        write(41,*)
                        write(41,*)"Xpc < 0.0d0",phase,com,dammy
                        write(41,*)
                    end if
                end do
            end do
                        
!!Phase Property
            call out_diffsx(F(i),F0_real)
            call calc_density(unknown,Po(i),Xpc,rho_i,Vp_i)
            
            if (F0_real <= 0.0d0)then
                call residualvectorset3(unknown*grid+BHin+BHout,rho_i(2))
                call residualvectorset3(unknown*grid+BHin+BHout,Vp_i(2))
            end if

!相対浸透率の求め方の決定
            call calc_permeability(unknown,V(i),F(i),L(i),W(i),Vp_i,krg(i),krf(i),kro(i),krw(i),Sg(i),Sf(i),So(i),Sw(i))
            !call calc_permeability_interpolation(unknown,time_roop,r,i,FM1(i),V(i),F(i),L(i),W(i),Vp_i,krg(i),krf(i),kro(i),krw(i),Sg(i),Sf(i),So(i),Sw(i))
            
            
            call calc_viscosity(unknown,Po(i),Xpc,Vp_i,FM1(i),Sw(i),So(i),mu_i)
 
            call out_diffsx(Sg(i),Sg0(i))
            call out_diffsx(Sf(i),Sf0(i))
            call out_diffsx(So(i),So0(i))
            call out_diffsx(Sw(i),Sw0(i))


            call out_diffsx(V(i),V0(i))
            call out_diffsx(L(i),L0(i))
            call out_diffsx(W(i),W0(i))
            
            do phase=1,s
                rho(phase,i)=rho_i(phase)
                Vp(phase,i)=Vp_i(phase)
                mu(phase,i)=mu_i(phase)
            end do

            call outxs(mu_i,mu_real)
            call outxs(rho_i,rho_real)

            !total mol
            ntotal=(rho(1,i)*V(i)+rho(2,i)*F(i)+rho(3,i)*L(i)+rho(4,i)*W(i))*Long*A/10.0d0
            call out_diffsx(ntotal,ntotal0)
            write(*,*)'ntotal',ntotal0

            do com=1,t
                Cgc(com,i)=Xpc(1,com)
                Cfc(com,i)=Xpc(2,com)
                Coc(com,i)=Xpc(3,com)
                Cwc(com,i)=Xpc(4,com)

                call out_diffsx(F(i),F0_real)
                if (F0_real <= 0.0d0)then
                    call residualvectorset3(unknown*grid+BHin+BHout,Cfc(com,i))
                    call residualvectorset3(unknown*grid+BHin+BHout,F(i))
                end if
            end do
            
        end do


        do i=1,grid
            kr_p(1,i)=krg(i)
            kr_p(2,i)=krf(i)
            kr_p(3,i)=kro(i)
            kr_p(4,i)=krw(i)
        end do

        if(r == roopmax)then
            call outxs(mu_i,mu_real)
            write(25,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",mu_real(1)
            write(26,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",mu_real(2)
            write(27,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",mu_real(3)
            write(28,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",mu_real(4)
            
            call outxs(rho_i,rho_real)
            write(29,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",rho_real(1)
            write(30,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",rho_real(2)
            write(31,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",rho_real(3)
            write(32,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",rho_real(4)
         
            call outxs(krg,kr_real)
            write(33,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",kr_real(1)
            call outxs(krf,kr_real)
            write(34,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",kr_real(1)
            call outxs(kro,kr_real)
            write(35,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",kr_real(1)
            call outxs(krw,kr_real)
            write(36,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",kr_real(1)
        end if
                
!?material balance(flow)-----------------------------------------------------------------------------
    do i=1,grid
        call calc_flow(time_roop,r,unknown,i,Nup,Nc,Nc0old,Po,Pbhin,Pbhout,rho,mu,kr_p,Cgc,Cfc,Coc,Cwc,qin0,qout0,faid,Poold,g_f)

        do com=1,t
            fxs(com+(i-1)*unknown+(m-1)*n)=g_f(com)
        end do

        if (i == 1) then
            do phase=1,s
                qin(phase)=qin0(phase)
            end do
        else if (i == grid) then
            do phase=1,s
                qout(phase)=qout0(phase)
            end do
        end if
    end do
    
!?saturation-------------------------------------------------------------------------------------
        do i=1,grid
            do com=1,t
                Nc_i(com)=Nc(com,i)
            end do
            do phase=1,s
                Vp_i(phase)=Vp(phase,i)
            end do
            call calc_saturation(unknown,V(i),F(i),L(i),W(i),Nc_i,Vp_i,g_s)
            fxs(i*unknown)=g_s
        end do
        
!?流量制御-----------------------------------------------------------------------------------------
        if (BHin+BHout /= 0) then
            allocate(g_fc(BHin+BHout))
            call calc_flow_control(unknown,qin,qout,g_fc)
        

            if (BHin == 1) then
                fxs(unknown*grid+1)=g_fc(1)
            end if
      
            if (BHin == 1 .AND. BHout == 1) then
                fxs(unknown*grid+2)=g_fc(2)
            else if (BHin /= 1 .AND. BHout == 1) then
                fxs(unknown*grid+1)=g_fc(1)
            end if
        end if
      
!?累計油生産量

        call calc_ORR(unknown,Cgc,Cfc,Coc,Cwc,rho,qout,Vgsc,Vosc)

    end subroutine calc_main
end module mod_main