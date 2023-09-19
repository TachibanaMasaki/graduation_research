module mod_flow
    use mod_initial
    use mod_autodiff
    use mod_density_W_inj
    use mod_density_G_inj
    contains 
  subroutine calc_flow(time_roop,r,unknown,i,Nup,Nc,Nc0old,Po,Pbhin,Pbhout,rho,mu,kr_p,Cgc,Cfc,Coc,Cwc,qin,qout,faid,Poold,g)
        implicit none
        integer,intent(in)::time_roop,r,unknown,i
        integer,intent(inout),dimension(t,grid)::Nup
        type(diffs),intent(in),dimension(grid)::Po
        type(diffs),intent(in),dimension(s,grid)::mu,rho,kr_p
        type(diffs),intent(in),dimension(t,grid)::Nc,Cgc,Cfc,Coc,Cwc
        type(diffs),intent(in)::Pbhin,Pbhout
        real(8),intent(in),dimension(t,grid)::Nc0old
        real(8),intent(in),dimension(grid)::faid,Poold
        real(8),allocatable,dimension(:,:)::Cc_ip1_r,Cc_i_r,Cc_im1_r
        real(8),allocatable,dimension(:)::Cc_SC,rhoc_SC,kr_ip1_r,kr_i_r,kr_im1_r
        type(diffs),intent(out),dimension(s)::qin,qout
        type(diffs),intent(out),dimension(t)::g(:)
        integer::com,phase,WFp,WFm,interval0,range0=1,num=1
        real(8)::dx,dy,dz,re,Wwell,Poim1,Poi,Poip1,dammy_real
        type(diffs)::transp(s,t),transm(s,t),fai(grid),kr_ip1(s),kr_i(s),kr_im1(s),Cc_ip1(s,t),Cc_i(s,t),Cc_im1(s,t)  &
                    ,dammy,lam,P1,Pgrid,rho_W_inj,rho_G_inj,Psc_diff,rho_W_SC,rho_G_SC

        allocate(Cc_SC(t),rhoc_SC(t),kr_ip1_r(s),kr_i_r(s),kr_im1_r(s),Cc_ip1_r(s,t),Cc_i_r(s,t),Cc_im1_r(s,t))


        dx=Long/grid
        dy=sqrt(A)  !(m)
        dz=sqrt(A)  !(m)

        Nup=0
 
        fai(i)=faid(i)*exp(cr*(Po(i)-Poold(i)))

        do com=1,t
            call residualvectorset3(unknown*grid+BHin+BHout,g(com))
        end do

!!!!!!!!!!!!!!!!!Middle conditions of grid i/=1,i/=grid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (i /= 1 .AND. i /= grid) then
            do phase=1,s
                kr_ip1(phase)=kr_p(phase,i+1)
            end do

            do phase=1,s
                kr_i(phase)=kr_p(phase,i)
            end do

            do phase=1,s
                kr_im1(phase)=kr_p(phase,i-1)
            end do

            do com=1,t
                Cc_ip1(1,com)=Cgc(com,i+1)
                Cc_ip1(2,com)=Cfc(com,i+1)
                Cc_ip1(3,com)=Coc(com,i+1)
                Cc_ip1(4,com)=Cwc(com,i+1)
            end do
            
            do com=1,t
                Cc_i(1,com)=Cgc(com,i)
                Cc_i(2,com)=Cfc(com,i)
                Cc_i(3,com)=Coc(com,i)
                Cc_i(4,com)=Cwc(com,i)
            end do

            do com=1,t
                Cc_im1(1,com)=Cgc(com,i-1)
                Cc_im1(2,com)=Cfc(com,i-1)
                Cc_im1(3,com)=Coc(com,i-1)
                Cc_im1(4,com)=Cwc(com,i-1)
            end do

            call out_diffsx(Po(i-1),Poim1)
            call out_diffsx(Po(i),Poi)
            call out_diffsx(Po(i+1),Poip1)

            if (Poi > Poip1) then
                WFp=1
            else
                WFp=0
            end if

            if (Poim1 > Poi) then
                WFm=1
            else 
                WFm=0
            end if

            do com=1,t
                do phase=1,s
                    transp(phase,com)=WFp*absk*kr_i(phase)/mu(phase,i)*rho(phase,i)*Cc_i(phase,com)  &
                                     +(1-WFp)*absk*kr_ip1(phase)/mu(phase,i+1)*rho(phase,i+1)*Cc_ip1(phase,com)
                end do
            end do

            do com=1,t
                do phase=1,s
                    transm(phase,com)=WFm*absk*kr_im1(phase)/mu(phase,i-1)*rho(phase,i-1)*Cc_im1(phase,com)  &
                                     +(1-WFm)*absk*kr_i(phase)/mu(phase,i)*rho(phase,i)*Cc_i(phase,com)
                end do
            end do


            do com=1,t
                do phase=1,s
                    g(com)=g(com)+1.0d0/dx**2.0d0*(transp(phase,com)*(Po(i+1)-Po(i))  &
                                                  -transm(phase,com)*(Po(i)-Po(i-1)))
                end do

                g(com)=g(com)-1.0d0/dt*(fai(i)*Nc(com,i)-faid(i)*Nc0old(com,i))
            end do

            do phase=1,s
                call out_diffsx(kr_ip1(phase),kr_ip1_r(phase))
                call out_diffsx(kr_i(phase),kr_i_r(phase))
                call out_diffsx(kr_im1(phase),kr_im1_r(phase))
            end do

            do com=1,t
                call out_diffsx(Cc_ip1(1,com),Cc_ip1_r(1,com))
                call out_diffsx(Cc_ip1(2,com),Cc_ip1_r(2,com))
                call out_diffsx(Cc_ip1(3,com),Cc_ip1_r(3,com))
                call out_diffsx(Cc_ip1(4,com),Cc_ip1_r(4,com))

                call out_diffsx(Cc_i(1,com),Cc_i_r(1,com))
                call out_diffsx(Cc_i(2,com),Cc_i_r(2,com))
                call out_diffsx(Cc_i(3,com),Cc_i_r(3,com))
                call out_diffsx(Cc_i(4,com),Cc_i_r(4,com))

                call out_diffsx(Cc_im1(1,com),Cc_im1_r(1,com))
                call out_diffsx(Cc_im1(2,com),Cc_im1_r(2,com))
                call out_diffsx(Cc_im1(3,com),Cc_im1_r(3,com))
                call out_diffsx(Cc_im1(4,com),Cc_im1_r(4,com))
            end do

            do com=1,t
                if((Poim1 > Poi .AND. Poi > Poip1) .AND.  &
                   (kr_im1_r(1)==0.0d0 .OR. Cc_im1_r(1,com)==0.0d0)  .AND.  &
                   (kr_im1_r(2)==0.0d0 .OR. Cc_im1_r(2,com)==0.0d0)  .AND.  &
                   (kr_im1_r(3)==0.0d0 .OR. Cc_im1_r(3,com)==0.0d0)  .AND.  &
                   (kr_im1_r(4)==0.0d0 .OR. Cc_im1_r(4,com)==0.0d0)  .AND.  &

                   (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0)      .AND.  &
                   (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0)      .AND.  &
                   (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0)      .AND.  &
                   (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                    !Nup(com,i)=1
                else if((Poim1 < Poi .AND. Poi < Poip1) .AND.  &
                        (kr_ip1_r(1)==0.0d0 .OR. Cc_ip1_r(1,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(2)==0.0d0 .OR. Cc_ip1_r(2,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(3)==0.0d0 .OR. Cc_ip1_r(3,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(4)==0.0d0 .OR. Cc_ip1_r(4,com)==0.0d0)  .AND.  &
 
                        (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0)      .AND.  &
                        (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0)      .AND.  &
                        (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0)      .AND.  &
                        (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                        !Nup(com,i)=1
                else if((Poim1 > Poi .AND. Poi < Poip1)                   .AND.  &
                        (kr_im1_r(1)==0.0d0 .OR. Cc_im1_r(1,com)==0.0d0)  .AND.  &
                        (kr_im1_r(2)==0.0d0 .OR. Cc_im1_r(2,com)==0.0d0)  .AND.  &
                        (kr_im1_r(3)==0.0d0 .OR. Cc_im1_r(3,com)==0.0d0)  .AND.  &
                        (kr_im1_r(4)==0.0d0 .OR. Cc_im1_r(4,com)==0.0d0)  .AND.  &
     
                        (kr_ip1_r(1)==0.0d0 .OR. Cc_ip1_r(1,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(2)==0.0d0 .OR. Cc_ip1_r(2,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(3)==0.0d0 .OR. Cc_ip1_r(3,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(4)==0.0d0 .OR. Cc_ip1_r(4,com)==0.0d0))then
                        !Nup(com,i)=1
                else if((Poim1 < Poi .AND. Poi > Poip1)               .AND.  &
                        (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0)  .AND.  &
                        (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0)  .AND.  &
                        (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0)  .AND.  &
                        (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                        !Nup(com,i)=1
                end if
            end do

!!!!!!!!!!!!!!!!!Boundary conditions of grid i=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (i == 1) then
            !坑井式
            do phase=1,s
                kr_ip1(phase)=kr_p(phase,i+1)
            end do

            do phase=1,s
                kr_i(phase)=kr_p(phase,i)
            end do

            do com=1,t
                Cc_ip1(1,com)=Cgc(com,i+1)
                Cc_ip1(2,com)=Cfc(com,i+1)
                Cc_ip1(3,com)=Coc(com,i+1)
                Cc_ip1(4,com)=Cwc(com,i+1)
            end do

            do com=1,t
                Cc_i(1,com)=Cgc(com,i)
                Cc_i(2,com)=Cfc(com,i)
                Cc_i(3,com)=Coc(com,i)
                Cc_i(4,com)=Cwc(com,i)
            end do

            call out_diffsx(Po(i),Poi)
            call out_diffsx(Po(i+1),Poip1)

            if (Poi > Poip1) then
                WFp=1
            else
                WFp=0
            end if

            do com=1,t
                do phase=1,s
                    transp(phase,com)=WFp*absk*kr_i(phase)/mu(phase,i)*rho(phase,i)*Cc_i(phase,com)  &
                                     +(1-WFp)*absk*kr_ip1(phase)/mu(phase,i+1)*rho(phase,i+1)*Cc_ip1(phase,com)
                end do
            end do
    
            re=0.14d0*sqrt(dx**2.0d0+dy**2.0d0)
            Wwell=2.0d0*4.0d0*atan(1.0d0)*absk*dz/(log(re/rin)+skinin)
            lam=kr_i(1)/mu(1,i)+kr_i(2)/mu(2,i)+kr_i(3)/mu(3,i)+kr_i(4)/mu(4,i)
            
            P1=Po(1)
            do Phase=1,s
                call residualvectorset3(unknown*grid+BHin+BHout,qin(Phase))
            end do
!water injection
            if (time_roop*dt <= change_time*range0)then
                if (BHin == 1) then
                    qin(4)=Wwell*lam*(Pbhin-Po(1))
                else
                    qin(4)=Wwell*lam*(P_bhin-Po(1))
                end if

                Cc_SC(1)=Cwc_SC1
                Cc_SC(2)=Cwc_SC2
                Cc_SC(3)=Cwc_SC3
                Cc_SC(4)=Cwc_SC4
                Cc_SC(5)=Cwc_SC5
                Cc_SC(6)=Cwc_SC6
                Cc_SC(7)=Cwc_SC7
                Cc_SC(8)=Cwc_SC8
                Cc_SC(9)=Cwc_SC9
                Cc_SC(10)=Cwc_SC10
                Cc_SC(11)=Cwc_SC11
            

                rhoc_SC(1)=rhoc_SC1
                rhoc_SC(2)=rhoc_SC2
                rhoc_SC(3)=rhoc_SC3
                rhoc_SC(4)=rhoc_SC4
                rhoc_SC(5)=rhoc_SC5
                rhoc_SC(6)=rhoc_SC6
                rhoc_SC(7)=rhoc_SC7
                rhoc_SC(8)=rhoc_SC8
                rhoc_SC(9)=rhoc_SC9
                rhoc_SC(10)=rhoc_SC10
                rhoc_SC(11)=rhoc_SC11

                call calc_density_W_inj(Po(i),rho_W_inj)

                do com=1,t
                    do phase=1,s
                        g(com)=g(com)+1.0d0/dx**2.0d0*transp(phase,com)*(Po(i+1)-Po(i))                     
                    end do
                    g(com)=g(com)+qin(4)/(dx*dy*dz)*Cc_SC(com)*rho_W_inj

                    g(com)=g(com)-1.0d0/dt*(fai(i)*Nc(com,i)-faid(i)*Nc0old(com,i))
                end do
                
!!!!!!!!!!!!!!Surface flow control!!!!!!!!!!!!!!!!!!!!!!!
                !call residualvectorset4(Psc,unknown*grid+BHin+BHout,Psc_diff)
                !call calc_density_W_inj(Psc_diff,rho_W_SC)
                !qin(4)=qin(4)*rho_W_inj/rho_W_SC
                
!CO2 injection
            else   
                if (BHin == 1) then
                    qin(1)=Wwell*lam*(Pbhin-P1)
                else
                    qin(1)=Wwell*lam*(P_bhin-P1)
                end if

                Cc_SC(1)=Cgc_SC1
                Cc_SC(2)=Cgc_SC2
                Cc_SC(3)=Cgc_SC3
                Cc_SC(4)=Cgc_SC4
                Cc_SC(5)=Cgc_SC5
                Cc_SC(6)=Cgc_SC6
                Cc_SC(7)=Cgc_SC7
                Cc_SC(8)=Cgc_SC8
                Cc_SC(9)=Cgc_SC9
                Cc_SC(10)=Cgc_SC10
                Cc_SC(11)=Cgc_SC11

                call calc_density_G_inj(unknown,Po(i),rho_G_inj)

                do com=1,t
                    do phase=1,s
                        g(com)=g(com)+1.0d0/dx**2.0d0*transp(phase,com)*(Po(i+1)-Po(i))                     
                    end do
                    g(com)=g(com)+qin(1)/(dx*dy*dz)*Cc_SC(com)*rho_G_inj
                    g(com)=g(com)-1.0d0/dt*(fai(i)*Nc(com,i)-faid(i)*Nc0old(com,i))
                end do           
                
!!!!!!!!!!!!!!Surface flow control!!!!!!!!!!!!!!!!!!!!!!!
                !call residualvectorset4(Psc,unknown*grid+BHin+BHout,Psc_diff)
                !call calc_density_G_inj(unknown,Psc_diff,rho_G_SC)
                !qin(1)=qin(1)*rho_G_inj/rho_G_SC
                
                interval0=time_roop*dt/change_time-2*num
                if(interval0 /= 0)then
                    range0=0
                else if(interval0 == 0 .AND. r == roopmax)then
                    range0=2*num+1
                    num=num+1
                end if
            end if

            
            do phase=1,s
                call out_diffsx(kr_ip1(phase),kr_ip1_r(phase))
                call out_diffsx(kr_i(phase),kr_i_r(phase))
            end do

            do com=1,t
                call out_diffsx(Cc_ip1(1,com),Cc_ip1_r(1,com))
                call out_diffsx(Cc_ip1(2,com),Cc_ip1_r(2,com))
                call out_diffsx(Cc_ip1(3,com),Cc_ip1_r(3,com))
                call out_diffsx(Cc_ip1(4,com),Cc_ip1_r(4,com))

                call out_diffsx(Cc_i(1,com),Cc_i_r(1,com))
                call out_diffsx(Cc_i(2,com),Cc_i_r(2,com))
                call out_diffsx(Cc_i(3,com),Cc_i_r(3,com))
                call out_diffsx(Cc_i(4,com),Cc_i_r(4,com))
            end do

            do com=1,t
                if((Poi > Poip1)                                     .AND.  &
                   (Cc_SC(com)==0.0d0)                               .AND.  &
                   (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0)      .AND.  &
                   (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0)      .AND.  &
                   (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0)      .AND.  &
                   (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                    !Nup(com,i)=1
                else if((Poi < Poip1)                                     .AND.  &
                        (Cc_SC(com)==0.0d0)                               .AND.  &
                        (kr_ip1_r(1)==0.0d0 .OR. Cc_ip1_r(1,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(2)==0.0d0 .OR. Cc_ip1_r(2,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(3)==0.0d0 .OR. Cc_ip1_r(3,com)==0.0d0)  .AND.  &
                        (kr_ip1_r(4)==0.0d0 .OR. Cc_ip1_r(4,com)==0.0d0))then
                        !Nup(com,i)=1
                end if
            end do
            
!!!!!!!!!!!!!!!!!Boundary conditions of grid i=grid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (i == grid) then
            do phase=1,s
                kr_i(phase)=kr_p(phase,i)
            end do

            do phase=1,s
                kr_im1(phase)=kr_p(phase,i-1)
            end do

            do com=1,t
                Cc_i(1,com)=Cgc(com,i)
                Cc_i(2,com)=Cfc(com,i)
                Cc_i(3,com)=Coc(com,i)
                Cc_i(4,com)=Cwc(com,i)
            end do

            do com=1,t
                Cc_im1(1,com)=Cgc(com,i-1)
                Cc_im1(2,com)=Cfc(com,i-1)
                Cc_im1(3,com)=Coc(com,i-1)
                Cc_im1(4,com)=Cwc(com,i-1)
            end do

            call out_diffsx(Po(i-1),Poim1)
            call out_diffsx(Po(i),Poi)

            if (Poim1 > Poi) then
                WFm=1
            else 
                WFm=0
            end if


            do com=1,t
                do phase=1,s
                    transm(phase,com)=WFm*absk*kr_im1(phase)/mu(phase,i-1)*rho(phase,i-1)*Cc_im1(phase,com)  &
                                     +(1-WFm)*absk*kr_i(phase)/mu(phase,i)*rho(phase,i)*Cc_i(phase,com)
                end do
            end do

            re=0.14d0*sqrt(dx**2.0d0+dy**2.0d0)
            Wwell=2.0d0*4.0d0*atan(1.0d0)*absk*dz/(log(re/rout)+skinout)

            Pgrid=Po(i)
            
            do Phase=1,s
                call residualvectorset3(unknown*grid+BHin+BHout,qout(Phase))
            end do                     
            
            if (BHout == 1) then
                do phase=1,s
                   qout(phase)=Wwell*kr_i(phase)/mu(phase,i)*(Pbhout-Pgrid)
                end do
            else
                do phase=1,s
                    qout(phase)=Wwell*kr_i(phase)/mu(phase,i)*(P_bhout-Pgrid)
                end do
            end if

            do com=1,t
                do phase=1,s      
                    g(com)=g(com)-1.0d0/dx**2.0d0*transm(phase,com)*(Po(i)-Po(i-1)) 
                end do
                    call residualvectorset3(unknown*grid+BHin+BHout,dammy)
                do phase=1,s
                    dammy=dammy+qout(phase)/(dx*dy*dz)*Cc_i(phase,com)*rho(phase,i)
                end do
                g(com)=g(com)+dammy
                g(com)=g(com)-1.0d0/dt*(fai(i)*Nc(com,i)-faid(i)*Nc0old(com,i))
            end do
        end if


        do phase=1,s
            call out_diffsx(kr_i(phase),kr_i_r(phase))
            call out_diffsx(kr_im1(phase),kr_im1_r(phase))
        end do

        do com=1,t
            call out_diffsx(Cc_i(1,com),Cc_i_r(1,com))
            call out_diffsx(Cc_i(2,com),Cc_i_r(2,com))
            call out_diffsx(Cc_i(3,com),Cc_i_r(3,com))
            call out_diffsx(Cc_i(4,com),Cc_i_r(4,com))

            call out_diffsx(Cc_im1(1,com),Cc_im1_r(1,com))
            call out_diffsx(Cc_im1(2,com),Cc_im1_r(2,com))
            call out_diffsx(Cc_im1(3,com),Cc_im1_r(3,com))
            call out_diffsx(Cc_im1(4,com),Cc_im1_r(4,com))
        end do

        do com=1,t
            if((Poim1 > Poi)                                     .AND.  &
               (kr_im1_r(1)==0.0d0 .OR. Cc_im1_r(1,com)==0.0d0)  .AND.  &
               (kr_im1_r(2)==0.0d0 .OR. Cc_im1_r(2,com)==0.0d0)  .AND.  &
               (kr_im1_r(3)==0.0d0 .OR. Cc_im1_r(3,com)==0.0d0)  .AND.  &
               (kr_im1_r(4)==0.0d0 .OR. Cc_im1_r(4,com)==0.0d0)  .AND.  &

               (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0)      .AND.  &
               (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0)      .AND.  &
               (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0)      .AND.  &
               (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                !Nup(com,i)=1
            else if((Poim1 < Poi)                                .AND.  &
                    (kr_i_r(1)==0.0d0 .OR. Cc_i_r(1,com)==0.0d0) .AND.  &
                    (kr_i_r(2)==0.0d0 .OR. Cc_i_r(2,com)==0.0d0) .AND.  &
                    (kr_i_r(3)==0.0d0 .OR. Cc_i_r(3,com)==0.0d0) .AND.  &
                    (kr_i_r(4)==0.0d0 .OR. Cc_i_r(4,com)==0.0d0))then
                    !Nup(com,i)=1
            end if
        end do

    end subroutine calc_flow
end module 