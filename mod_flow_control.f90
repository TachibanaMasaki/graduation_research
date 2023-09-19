module mod_flow_control
    use mod_initial
    use mod_autodiff
    implicit none
    contains
    subroutine calc_flow_control(unknown,qin,qout,fxs_fc)
        integer,intent(in)::unknown
        type(diffs),intent(in),dimension(s)::qin,qout
        type(diffs),allocatable,intent(out)::fxs_fc(:)
        integer::phase
        real(8),allocatable,dimension(:)::qin_r,qout_r

        allocate(qin_r(s),qout_r(s))

        call outxs(qin,qin_r)
        call outxs(qout,qout_r)

        write(*,*)qin_r
        write(*,*)qout_r

        if (BHin+BHout /= 0)then
            allocate(fxs_fc(BHin+BHout))

            if (BHin == 1) then
                call residualvectorset3(unknown*grid+BHin+BHout,fxs_fc(1))
                do phase=1,s
                    fxs_fc(1)=fxs_fc(1)+qin(phase)
                end do
                fxs_fc(1)=fxs_fc(1)-Qtotalin
            end if
      
            if (BHin == 1 .AND. BHout == 1) then
                call residualvectorset3(unknown*grid+BHin+BHout,fxs_fc(2))
                do phase=1,s
                    fxs_fc(2)=fxs_fc(2)+qout(phase)
                end do
                fxs_fc(2)=fxs_fc(2)+Qtotalout
            else if (BHin /= 1 .AND. BHout == 1) then
                call residualvectorset3(unknown*grid+BHin+BHout,fxs_fc(1))
                do phase=1,s
                    fxs_fc(1)=fxs_fc(1)+qout(phase)
                end do
                fxs_fc(1)=fxs_fc(1)+Qtotalout
            end if
       end if   
    end subroutine calc_flow_control
end module mod_flow_control