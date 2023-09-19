module mod_permeability
    use mod_initial
    use mod_autodiff
    contains
    subroutine calc_permeability(unknown,V0,F0,L0,W0,Vp,krg,krf,kro,krw,Sg,Sf,So,Sw)
      implicit none
      integer,intent(in)::unknown
      type(diffs),intent(in)::V0,F0,L0,W0
      type(diffs),intent(in),dimension(s)::Vp
      type(diffs),intent(out)::krg,krf,kro,krw,Sg,Sf,So,Sw
      real(8)::Sg_real,Sf_real,So_real,Sw_real,krg_real,krf_real,kro_real,krw_real,Sgc_now,Sfc_now

      Sg=Vp(1)*V0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
      Sf=Vp(2)*F0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
      So=Vp(3)*L0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)
      Sw=Vp(4)*W0/(Vp(1)*V0+Vp(2)*F0+Vp(3)*L0+Vp(4)*W0)

      call out_diffsx(Sg,Sg_real)
      call out_diffsx(Sf,Sf_real)
      call out_diffsx(So,So_real)
      call out_diffsx(Sw,Sw_real)

      if (Sgc < Sg_real)then
          Sgc_now=Sgc
      else
          Sgc_now=Sg_real
      end if
      
      if (Sfc < Sf_real)then
          Sfc_now=Sfc
      else
          Sfc_now=Sf_real
      end if  

!ガス相の相対浸透率
      krg=cogperm*((Sg-Sgc_now)/(1.0d0-Sgc_now-Sfc_now-Sor-Swc))**epgperm
      call out_diffsx(krg,krg_real)
      if(krg_real < 0.0d0)then
          call residualvectorset3(unknown*grid+BHin+BHout,krg)
      end if
!foam相の相対浸透率
      krf=cofperm*((Sf-Sfc_now)/(1.0d0-Sgc_now-Sfc_now-Sor-Swc))**epfperm
      call out_diffsx(krf,krf_real)
      if(krf_real < 0.0d0)then
          call residualvectorset3(unknown*grid+BHin+BHout,krf)
      end if
!油相の相対浸透率
      kro=cooperm*((So-Sor)/(1.0d0-Sgc_now-Sfc_now-Sor-Swc))**epoperm
      call out_diffsx(kro,kro_real)
      if(kro_real < 0.0d0)then
          call residualvectorset3(unknown*grid+BHin+BHout,kro)
      end if
!水相の相対浸透率
      krw=cowperm*((Sw-Swc)/(1.0d0-Sgc_now-Sfc_now-Sor-Swc))**epwperm
      call out_diffsx(krw,krw_real)
      if(krw_real < 0.0d0)then
          call residualvectorset3(unknown*grid+BHin+BHout,krw)
      end if
     
      if (krg_real > 1.0d0)then
          write(35,*)"krg",krg_real
      end if
      if (krf_real > 1.0d0)then
          write(35,*)"krf",krf_real
      end if
      if (kro_real > 1.0d0)then
          write(35,*)"kro",kro_real
      end if
      if (krw_real > 1.0d0)then
          write(35,*)"krw",krw_real
      end if
      
  end subroutine calc_permeability
end module