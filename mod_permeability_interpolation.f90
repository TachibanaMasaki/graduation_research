module mod_permeability_interpolation
    use mod_initial
    use mod_autodiff
    contains
    subroutine calc_permeability_interpolation(unknown,time_roop,r,i_grid,FM,V0,F0,L0,W0,Vp,krg,krf,kro,krw,Sg,Sf,So,Sw)
        implicit none
        integer,intent(in)::unknown,time_roop,r,i_grid
        type(diffs),intent(in)::FM,V0,F0,L0,W0
        type(diffs),intent(in),dimension(s)::Vp
        type(diffs),intent(out)::krg,krf,kro,krw,Sg,Sf,So,Sw
        real(8),dimension(data_ow)::swd,krwd,krowd
        real(8),dimension(data_go)::sgd,krgd,krogd
        integer::i
        real(8)::Sg_real,Sf_real,So_real,Sw_real,krocw
        type(diffs)::krow,krog

        call residualvectorset3(unknown*grid+BHin+BHout,krg)
        call residualvectorset3(unknown*grid+BHin+BHout,krf)
        call residualvectorset3(unknown*grid+BHin+BHout,kro)
        call residualvectorset3(unknown*grid+BHin+BHout,krw)
        
        if (time_roop == 1 .AND. r == 1 .AND. i_grid == 1) then
            do i=1,data_ow
                read(2,*)swd(i),krwd(i),krowd(i)
            end do
            do i=1,data_go
                read(3,*)sgd(i),krgd(i),krogd(i)
            end do
        end if
        
        Sg=Vp(1)*V0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
        Sf=Vp(2)*F0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
        So=Vp(3)*L0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
        Sw=Vp(4)*W0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)

        call out_diffsx(Sg,Sg_real)
        call out_diffsx(Sf,Sf_real)
        call out_diffsx(So,So_real)
        call out_diffsx(Sw,Sw_real)

!ÉKÉXëäÇÃëäëŒêZìßó¶ (ê¸å`ï‚ä‘)
        do i=1,data_go-1
            if(Sg_real <= sgd(i+1))then
                krg=(krgd(i+1)-krgd(i))/(sgd(i+1)-sgd(i))*(Sg-sgd(i))+krgd(i)
                krog=(krogd(i+1)-krogd(i))/(sgd(i+1)-sgd(i))*(Sg-sgd(i))+krogd(i)
                exit
            endif
        end do


!foamëäÇÃëäëŒêZìßó¶(ï‚ê≥)
        krf=FM*krg

!êÖëäÇÃëäëŒêZìßó¶ (ê¸å`ï‚ä‘)
      
        do i=1,data_ow-1
            if(Sw_real <= swd(i+1))then
                krw=(krwd(i+1)-krwd(i))/(swd(i+1)-swd(i))*(Sw-swd(i))+krwd(i)
                krow=(krowd(i+1)-krowd(i))/(swd(i+1)-swd(i))*(Sw-swd(i))+krowd(i)
                exit
            endif
        end do

!ñ˚ëäÇÃëäëŒêZìßó¶
        do i=1,data_ow-1
            if(Swc <= swd(i+1))then
                krocw=(krowd(i+1)-krowd(i))/(swd(i+1)-swd(i))*(Swc-swd(i))+krowd(i)
                exit
            endif
        end do

!stoneëÊàÍñ@ë•
        !kro=krow*krog*((So-Sor)/(1.0d0-Swc-Sor))/krocw/(1.0d0-(Sw-Swc)/(1.0d0-Swc-Sor))/(1.0d0-Sg/(1.0d0-Swc-Sor))
!stoneëÊìÒñ@ë•        
        kro=krocw*((krow/krocw+krw)*(krog/krocw+krg)-(krw+krg))
   
  end subroutine calc_permeability_interpolation
end module