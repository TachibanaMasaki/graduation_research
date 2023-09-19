program main
    use mod_initial
    use mod_initial_Value
    use mod_main
    use mod_pvgauss
    implicit none
    integer::unknown,BH,com,c,r,k,time_roop,time,i,j,conv,interval,range=1,com_out
    integer,allocatable,dimension(:,:)::Nup
    real(8),allocatable,dimension(:,:)::Nc0,Nc0old,GZ
    real(8),allocatable,dimension(:,:)::Kgo0,Kwo0,lnKgo0,lnKwo0
    real(8),allocatable,dimension(:)::Ncin
    real(8),allocatable,dimension(:)::Kgoin,Kwoin,K0
    real(8),allocatable,dimension(:)::Po0,Poold,faiold
    real(8),allocatable,dimension(:)::V0,F0,L0,W0,Sg0,Sf0,So0,Sw0,ntotal
    real(8),allocatable,dimension(:,:)::jacobian
    real(8),allocatable,dimension(:)::residual
    real(8)::Qtotal,PV,ORR,GRR,Vgsc_ini,Vosc_ini,Vgsc_total=0.0d0,Vosc_total=0.0d0,Vgsc,Vosc,Pbh0in,Pbh0out,error
    type(diffs),allocatable::fxs(:)

    open(1,file="data.txt")
    open(2,file="permeability_o-w.txt")
    open(3,file="permeability_g-o.txt")
    
    open(4,file="result_saturation-gas.csv")
    open(5,file="result_saturation-foam.csv")
    open(6,file="result_saturation-oil.csv")
    open(7,file="result_saturation-water.csv")
    
    open(8,file="result_all_saturation-gas.csv")
    open(9,file="result_all_saturation-foam.csv")
    open(10,file="result_all_saturation-oil.csv")
    open(11,file="result_all_saturation-water.csv")
    
    open(12,file="result_pressure.csv")
    open(13,file="result_BH-pressure.csv")
    
    open(14,file="result_mole-C1.csv")
    open(15,file="result_mole-C2.csv")
    open(16,file="result_mole-C3.csv")
    open(17,file="result_mole-C4.csv")
    open(18,file="result_mole-C5.csv")
    open(19,file="result_mole-C6.csv")
    open(20,file="result_mole-C7.csv")
    open(21,file="result_mole-CO2.csv")
    open(22,file="result_mole-H2O.csv")
    open(23,file="result_mole-SDS.csv")
    open(24,file="result_mole-NaCl.csv")  
    
    open(25,file="result_viscosity_gas.csv")
    open(26,file="result_viscosity_foam.csv")
    open(27,file="result_viscosity_oil.csv")
    open(28,file="result_viscosity_water.csv")
    
    open(29,file="result_density_gas.csv")
    open(30,file="result_density_foam.csv")
    open(31,file="result_density_oil.csv")
    open(32,file="result_density_water.csv")
    
    open(33,file="result_permeability_gas.csv")
    open(34,file="result_permeability_foam.csv")
    open(35,file="result_permeability_oil.csv")
    open(36,file="result_permeability_water.csv")
    
    open(37,file="result_prod_GRR.csv")
    open(38,file="result_prod_ORR.csv")
    open(39,file="result_CumGas_SC.csv")
    open(40,file="result_CumOil_SC.csv")
    
    open(41,file="result_check.csv")
    
    open(42,file="result_Phase_check.csv")
    allocate(Nup(t,grid))
    allocate(Ncin(t),Kgoin(n),Kwoin(n),K0((m-1)*n),F0(grid),V0(grid),L0(grid),W0(grid))
    allocate(Nc0(t,grid),Nc0old(t,grid),GZ(t,grid),Kgo0(n,grid),Kwo0(n,grid),Po0(grid),Poold(grid),faiold(grid))
    allocate(Sg0(grid),Sf0(grid),So0(grid),Sw0(grid),ntotal(grid))
    allocate(lnKgo0(n,grid),lnKwo0(n,grid))
    unknown=(m-1)*n+t+1

    BH=BHin+BHout
    allocate(jacobian(unknown*grid+BH,unknown*grid+BH),residual(unknown*grid+BH))


    call calc_initial_Value(Nc0,Nc0old,Kgo0,Kwo0,V0,F0,L0,W0,Pbh0in,Pbh0out,Po0,Poold,faiold,Vgsc_ini,Vosc_ini,K0)

    
    ntotal=0.0d0
    do i=1,grid
        do c=1,n
            ntotal(i)=ntotal(i)+Nc0(c,i)
        end do
    end do

    do i=1,grid
        do c=1,n
            GZ(c,i)=Nc0(c,i)/ntotal(i)
        end do
    end do
 
    write(14,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(1,1),",",GZ(1,3),",",GZ(1,5)
    write(15,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(2,1),",",GZ(2,3),",",GZ(2,5)
    write(16,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(3,1),",",GZ(3,3),",",GZ(3,5)
    write(17,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(4,1),",",GZ(4,3),",",GZ(4,5)
    write(18,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(5,1),",",GZ(5,3),",",GZ(5,5)
    write(19,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(6,1),",",GZ(6,3),",",GZ(6,5)
    write(20,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(7,1),",",GZ(7,3),",",GZ(7,5)
    write(21,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(8,1),",",GZ(8,3),",",GZ(8,5)
    write(22,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',GZ(9,1),",",GZ(9,3),",",GZ(9,5)
    write(23,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',0.0d0,",",0.0d0,",",0.0d0
    write(24,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')0.0d0,',',0.0d0,",",0.0d0,",",0.0d0
    
    write(41,*)"Simulation start\(^_^)/"
    time=day/dt
    do time_roop=1,time
        write(*,*)
        write(*,*)"time roop",time_roop

    do r=1,roopmax
        write(*,*)"raphson roop=",r
        call calc_main(time_roop,r,Nup,V0,L0,W0,Sg0,Sf0,So0,Sw0,Kgo0,Kwo0,Nc0,Nc0old,Po0,Pbh0in,Pbh0out,faiold,Poold,Vgsc,Vosc,fxs)
        call jacobian_mat(fxs,jacobian)
        call outxs(fxs,residual)

        do k=1,grid
            do i=1,(m-1)*n
                if (K0(i) == 0.0d0) then
                  do j=1,unknown*grid+BH
                      jacobian(i+(k-1)*unknown,j)=0.0d0
                      jacobian(j,i+(k-1)*unknown)=0.0d0
                  end do
                  jacobian(i+(k-1)*unknown,i+(k-1)*unknown)=1.0d0
                  residual(i+(k-1)*unknown)=0.0d0
                end if
             end do
        end do

        do i=1,grid
            do com=1,t
                if (Nup(com,i) == 1) then
                    do j=1,unknown*grid+BH
                        jacobian(com+(i-1)*unknown+(m-1)*n,j)=0.0d0
                        jacobian(j,com+(i-1)*unknown+(m-1)*n)=0.0d0
                    end do
                      jacobian(com+(i-1)*unknown+(m-1)*n,com+(i-1)*unknown+(m-1)*n)=1.0d0
                      residual(com+(i-1)*unknown+(m-1)*n)=0.0d0
                end if
             end do
        end do

!!!!!!!!!!!!!!!!!!!!!!存在しない物質を除外!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (time_roop*dt <= change_time)then
            com_out=8
            do i=1,grid
                do j=1,unknown*grid+BH
                    jacobian(com_out+(i-1)*unknown+(m-1)*n,j)=0.0d0
                    jacobian(j,com_out+(i-1)*unknown+(m-1)*n)=0.0d0
                end do
                jacobian(com_out+(i-1)*unknown+(m-1)*n,com_out+(i-1)*unknown+(m-1)*n)=1.0d0
                residual(com_out+(i-1)*unknown+(m-1)*n)=0.0d0
            end do
        end if
        !com_out=10
        !do i=1,grid
        !    do j=1,unknown*grid+BH
        !        jacobian(com_out+(i-1)*unknown+(m-1)*n,j)=0.0d0
        !        jacobian(j,com_out+(i-1)*unknown+(m-1)*n)=0.0d0
        !    end do
        !    jacobian(com_out+(i-1)*unknown+(m-1)*n,com_out+(i-1)*unknown+(m-1)*n)=1.0d0
        !    residual(com_out+(i-1)*unknown+(m-1)*n)=0.0d0
        !end do
        
        !com_out=11
        !do i=1,grid
        !    do j=1,unknown*grid+BH
        !        jacobian(com_out+(i-1)*unknown+(m-1)*n,j)=0.0d0
        !        jacobian(j,com_out+(i-1)*unknown+(m-1)*n)=0.0d0
        !    end do
        !    jacobian(com_out+(i-1)*unknown+(m-1)*n,com_out+(i-1)*unknown+(m-1)*n)=1.0d0
        !    residual(com_out+(i-1)*unknown+(m-1)*n)=0.0d0
        !end do  
     
!!流動のみ--------------------------------------------------------------------------
        !do k=1,grid
        !    do i=1,(m-1)*n
        !          do j=1,unknown*grid+BH
        !              jacobian(i+(k-1)*unknown,j)=0.0d0
        !              jacobian(j,i+(k-1)*unknown)=0.0d0
        !          end do
        !          jacobian(i+(k-1)*unknown,i+(k-1)*unknown)=1.0d0
        !          residual(i+(k-1)*unknown)=0.0d0
        !     end do
        !end do
!!流動なし---------------------------------------------------------------------------
        !do k=1,grid
        !    do com=1,t+1
        !        do j=1,unknown*grid+BH
        !            jacobian(com+(k-1)*unknown+(m-1)*n,j)=0.0d0
        !            jacobian(j,com+(k-1)*unknown+(m-1)*n)=0.0d0
        !        end do
        !          jacobian(com+(k-1)*unknown+(m-1)*n,com+(k-1)*unknown+(m-1)*n)=1.0d0
        !          residual(com+(k-1)*unknown+(m-1)*n)=0.0d0
        !     end do
        !end do
!!圧力なし----------------------------------------------------------------------------
        !do k=1,grid
        !    do j=1,unknown*grid+BH
        !        jacobian(k*unknown,j)=0.0d0
        !        jacobian(j,k*unknown)=0.0d0
        !    end do
        !    jacobian(k*unknown,k*unknown)=1.0d0
        !    residual(k*unknown)=0.0d0
        !end do
!!-----------------------------------------------------------------------------------

        residual=-residual

        call pvgauss(unknown*grid+BH,jacobian,residual)

        do c=1,n
            do i=1,grid
                lnKgo0(c,i)=log(Kgo0(c,i))+residual(unknown*(i-1)+c)
                Kgo0(c,i)=exp(lnKgo0(c,i))

                lnKwo0(c,i)=log(Kwo0(c,i))+residual(unknown*(i-1)+c+n)
                Kwo0(c,i)=exp(lnKwo0(c,i))
            end do
        end do

        do com=1,t
            do i=1,grid
                Nc0(com,i)=Nc0(com,i)+residual(unknown*(i-1)+com+2*n)
            end do
        end do

        do com=1,t
            do i=1,grid
                if(Nc0(com,i) < 0.0d0)then
                    write(*,*)"Nc0<0",com,i,Nc0(com,i)
                end if
            end do
        end do
        
        do i=1,grid
            Po0(i)=Po0(i)+residual(unknown*i)
        end do

        
        if (BHin == 1) then
            Pbh0in=Pbh0in+residual(unknown*grid+1)
        end if

        if (BHin == 1 .AND. BHout == 1) then
            Pbh0out=Pbh0out+residual(unknown*grid+2)
        else if (BHin /= 1 .AND. BHout == 1) then
            Pbh0out=Pbh0out+residual(unknown*grid+1)
        end if

        do i=1,grid
            residual(unknown*i)=residual(unknown*i)
        end do

        if (BHin == 1) then
            residual(unknown*grid+1)=residual(unknown*grid+1)
        end if

        if (BHin == 1 .AND. BHout == 1) then
            residual(unknown*grid+2)=residual(unknown*grid+2)
        else if (BHin /= 1 .AND. BHout == 1) then
            residual(unknown*grid+1)=residual(unknown*grid+1)
        end if

        
        error=dot_product(residual,residual)
        error=sqrt(error)

        write(*,*)"error=",error

        if(error<epsilon)then
            conv=r
            write(*,*)"convergence",r
            exit  
        end if
    end do

!!!!!!!!!!!値の更新!!!!!!!!!!!!!!!!!!!!!(Poold=Po0の場所おかしい faioldに影響)
    Poold=Po0
    Nc0old=Nc0
    do i=1,grid
        faiold(i)=faiold(i)*exp(cr*(Po0(i)-Poold(i)))
    end do
    
    Vgsc_total=Vgsc_total+Vgsc*dt
    Vosc_total=Vosc_total+Vosc*dt
    GRR=Vgsc_total/Vgsc_ini*100.0d0
    ORR=Vosc_total/Vosc_ini*100.0d0
    
!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!!!
    
    interval=time_roop*dt/(24*60*60)-range

    if (interval==0)then
        range=range+1
!!PVの計算
        Qtotal=Qtotalin*time_roop*dt
        PV=Qtotal/(Long*A*poro)

        write(37,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",PV,",",GRR
        write(38,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",PV,",",ORR
        write(39,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",PV,",",Vgsc_total
        write(40,*)time_roop*dt/(24.0d0*60.0d0*60.0d0),",",PV,",",Vosc_total

        
        write(4,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Sg0(1),',',Sg0(3),',',Sg0(5)
        write(5,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Sf0(1),',',Sf0(3),',',Sf0(5)
        write(6,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',So0(1),',',So0(3),',',So0(5)
        write(7,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Sw0(1),',',Sw0(3),',',Sw0(5)
        
        do i=1,grid
            write(8,*)Long/grid*(i-0.5),',',Sg0(i)
            write(9,*)Long/grid*(i-0.5),',',Sf0(i)
            write(10,*)Long/grid*(i-0.5),',',So0(i)
            write(11,*)Long/grid*(i-0.5),',',Sw0(i)
        end do
        
        write(12,'(f10.5,A,f10.5,A,f10.5,A,f10.5)')  &
            time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Po0(1)/6894.757d0,',',Po0(3)/6894.757d0,",",Po0(5)/6894.757d0  ![Pa]→[psi]
        write(13,'(f10.5,A,f10.5,A,f10.5,A,f10.5)')  &
            time_roop*dt/(24.0d0*60.0d0*60.0d0),',',Pbh0in/6894.757d0,',',Pbh0out/6894.757d0
        
          
        ntotal=0.0d0
        do i=1,grid
            do c=1,n
                ntotal(i)=ntotal(i)+Nc0(c,i)
            end do
        end do

        do i=1,grid
            do c=1,n
                GZ(c,i)=Nc0(c,i)/ntotal(i)
            end do
        end do
        
        write(14,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(1,1),",",GZ(1,3),",",GZ(1,5)
        write(15,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(2,1),",",GZ(2,3),",",GZ(2,5)
        write(16,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(3,1),",",GZ(3,3),",",GZ(3,5)
        write(17,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(4,1),",",GZ(4,3),",",GZ(4,5)
        write(18,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(5,1),",",GZ(5,3),",",GZ(5,5)
        write(19,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(6,1),",",GZ(6,3),",",GZ(6,5)
        write(20,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(7,1),",",GZ(7,3),",",GZ(7,5)
        write(21,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(8,1),",",GZ(8,3),",",GZ(8,5)
        write(22,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(9,1),",",GZ(9,3),",",GZ(9,5)
        write(23,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(10,1),",",GZ(10,3),",",GZ(10,5)
        write(24,'(f10.5,A,f10.7,A,f10.7,A,f10.7)')time_roop*dt/(24.0d0*60.0d0*60.0d0),',',GZ(11,1),",",GZ(11,3),",",GZ(11,5)
        
       
    end if
  
  end do

end program main