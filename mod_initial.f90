module mod_initial
    implicit none
!平衡状態にある相の数 p
    integer,parameter::m=3
!平衡状態にある成分の数 c
    integer,parameter::n=9
!すべての相の数 phase
    integer,parameter::s=4
!すべての成分の数 com
    integer,parameter::t=11
    
    integer,parameter::m_1=1
!grid数
    integer,parameter::grid=5
    real(8),parameter::grid_real=5.0d0
!誤差
    real(8),parameter::epsilon=10.0d0**(-30.0d0)
!newtonroopMax
    integer,parameter::roopmax=50

!地上状態の圧力 (Pa)
    real(8),parameter::Psc=14.7d0*6894.757d0  ![psi]→[Pa]  14.7 psi
!地上状態の温度 (K)
    real(8),parameter::Tsc=15.5556d0+273.15d0  ![℃]→[K]  !60 F

!全圧 (Pa)
    real(8),parameter::Pinitial=2000.0d0*6894.757d0  ![psi]→[Pa]  2000 psi   
!温度 (K)
    real(8),parameter::T0=80.000d0+273.15d0  ![℃]→[K] !176 F

!!condition
!時間 (s)
    integer,parameter::day=2000*24*60*60  ![s]  !2000 day
!dt (s)
    !integer,parameter::dt=6*60*60  ![s]  !6 hour
    !real(8),parameter::dt_real=6.0d0*60.0d0*60.0d0  ![s]  !6 hour
    integer,parameter::dt=24*60*60  ![s]  !24 hour
    real(8),parameter::dt_real=24.0d0*60.0d0*60.0d0  ![s]  !24 hour
    
!!!!!!!!!!!!!!!!!!!!!!!!!!圧入井の条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!圧入井の坑井圧力 (Pa)
    !!流量制御をするときは下のBHinを1にし、坑底圧力制御をするときは0にする。
    integer,parameter::BHin=1
    real(8),parameter::Qtotalin=5000.0d0*2.831685d0*10**(-2.0d0)/(24*60*60)  ![ft^3/day]→[m^3/s]
    real(8),parameter::P_bhin=2000.0d0*6894.757d0  ![psi]→[Pa]  2000 psi
    
!WAGで切り替えるタイミング (s)
    integer,parameter::change_time=60*24*60*60  ![s]  60 day
    
!圧入井の坑井の半径 (m)
    real(8),parameter::rin=0.28d0*0.3048d0  ![ft]→[m]  0.28 ft

    real(8),parameter::skinin=0.0d0

!水圧入におけるモル分率
    real(8),parameter::Cwc_SC1=0.0d0
    real(8),parameter::Cwc_SC2=0.0d0
    real(8),parameter::Cwc_SC3=0.0d0
    real(8),parameter::Cwc_SC4=0.0d0
    real(8),parameter::Cwc_SC5=0.0d0
    real(8),parameter::Cwc_SC6=0.0d0
    real(8),parameter::Cwc_SC7=0.0d0
    real(8),parameter::Cwc_SC8=0.0d0
    real(8),parameter::Cwc_SC9=0.995192156d0
    real(8),parameter::Cwc_SC10=0.00012652d0
    real(8),parameter::Cwc_SC11=0.004681324d0

!ガス圧入におけるモル分率
    real(8),parameter::Cgc_SC1=0.0d0
    real(8),parameter::Cgc_SC2=0.0d0
    real(8),parameter::Cgc_SC3=0.0d0
    real(8),parameter::Cgc_SC4=0.0d0
    real(8),parameter::Cgc_SC5=0.0d0
    real(8),parameter::Cgc_SC6=0.0d0
    real(8),parameter::Cgc_SC7=0.0d0
    real(8),parameter::Cgc_SC8=1.0d0
    real(8),parameter::Cgc_SC9=0.0d0
    real(8),parameter::Cgc_SC10=0.0d0
    real(8),parameter::Cgc_SC11=0.0d0

    
!圧入量におけるモル密度(mol/m^3)
    real(8),parameter::rhoc_SC1=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC2=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC3=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC4=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC5=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC6=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC7=0.0d0!使わないから無視
    real(8),parameter::rhoc_SC8=1.0d0/(23.54d0*10.0d0**(-3.0d0))
    real(8),parameter::rhoc_SC9=1.0d0/(0.02112d0*10.0d0**(-3.0d0))
    real(8),parameter::rhoc_SC10=1.01d0/288.38d0*10**6.0d0
    real(8),parameter::rhoc_SC11=2.16d0/58.44d0*10**6.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!生産井の条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!生産井の坑井圧力 (Pa)
    !!流量制御をするときは下のBHoutを1にし、坑底圧力制御をするときは0にする。
    integer,parameter::BHout=0
    real(8),parameter::Qtotalout=5000.0d0*2.831685d0*10**(-2.0d0)/(24*60*60)  ![ft^3/day]→[m^3/s]
    real(8),parameter::P_bhout=2000.0d0*6894.757d0  ![psi]→[Pa]  2000 psi
        
!生産井の坑井の半径 (m)
    real(8),parameter::rout=0.28d0*0.3048d0  ![ft]→[m]  0.28 ft

    real(8),parameter::skinout=0.0d0
!!貯留層の特性
!Length of the system (m)
    real(8),parameter::Long=1000.0d0*0.3048d0  ![ft]→[m]  1000 ft
!Cross-sectional area (m^2)
    real(8),parameter::A=10000.0d0*0.3048d0**2.0d0  ![ft^2]→[m^2]  10000 ft2
!Porosity
    real(8),parameter::poro=0.2d0
    real(8),parameter::poro0=0.2d0
    real(8),parameter::cr=3.0d0*10**(-6.0d0)/(6.894757d0*10**3.0d0)  ![1/Pa]
!Permeability (m^2)
    real(8),parameter::absk=3000.0d0*9.86923d0*10**(-16.0d0)  ![mD]→[m^2]
    
!不動水飽和率    
    real(8),parameter::Sgc=0.1d0
    real(8),parameter::Sfc=0.3d0
    real(8),parameter::Sor=0.2d0
    real(8),parameter::Swc=0.15d0
!線形補間
    integer,parameter::data_ow=20
    integer,parameter::data_go=20
!!成分の特性
!MW(c) = c成分の物質量 (g/mol)
    real(8),parameter::MW1=16.043d0 
    real(8),parameter::MW2=30.07d0
    real(8),parameter::MW3=44.097d0
    real(8),parameter::MW4=58.124d0
    real(8),parameter::MW5=72.151d0
    real(8),parameter::MW6=86.0d0
    real(8),parameter::MW7=190.0d0
    real(8),parameter::MW8=44.01d0
    real(8),parameter::MW9=18.015d0
    real(8),parameter::MW10=288.31d0
    real(8),parameter::MW11=58.44d0

! Pc(c) = c成分の臨界圧力 (Pa)
    real(8),parameter::Pc1=45.4d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc2=48.2d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc3=41.9d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc4=37.5d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc5=33.3d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc6=32.46d0*101325.0d0 ![atm]→[Pa] 
    real(8),parameter::Pc7=19.33d0*101325.0d0 ![atm]→[Pa] 
    real(8),parameter::Pc8=72.8d0*101325.0d0  ![atm]→[Pa] 
    real(8),parameter::Pc9=217.6d0*101325.0d0 ![atm]→[Pa] 

!Tc(c) = c成分の臨界温度 (K)
    real(8),parameter::Tc1=190.6d0 
    real(8),parameter::Tc2=305.4d0 
    real(8),parameter::Tc3=369.8d0 
    real(8),parameter::Tc4=425.2d0 
    real(8),parameter::Tc5=469.6d0 
    real(8),parameter::Tc6=507.5d0 
    real(8),parameter::Tc7=700.7d0 
    real(8),parameter::Tc8=304.2d0
    real(8),parameter::Tc9=647.3d0

! Vc(c) = c成分の臨界体積 (m^3/mol)
    real(8),parameter::Vc1=0.000099d0
    real(8),parameter::Vc2=0.000148d0
    real(8),parameter::Vc3=0.000203d0
    real(8),parameter::Vc4=0.000255d0
    real(8),parameter::Vc5=0.000304d0
    real(8),parameter::Vc6=0.000344d0
    real(8),parameter::Vc7=0.000723d0
    real(8),parameter::Vc8=0.000094d0
    real(8),parameter::Vc9=0.000056d0
    
!i成分とj成分における相互作用係数(-)
    real(8)::Theta18=0.105d0
    real(8)::Theta19=0.4907d0
    real(8)::Theta28=0.13d0
    real(8)::Theta29=0.4911d0
    real(8)::Theta38=0.125d0
    real(8)::Theta39=0.5469d0
    real(8)::Theta48=0.115d0
    real(8)::Theta49=0.508d0
    real(8)::Theta58=0.115d0
    real(8)::Theta59=0.5d0
    real(8)::Theta68=0.115d0
    real(8)::Theta69=0.48d0
    real(8)::Theta78=0.115d0
    real(8)::Theta79=0.48d0
    real(8)::Theta89=0.2d0

! omega = 偏心因子
    real(8),parameter::omega1=0.008d0
    real(8),parameter::omega2=0.098d0
    real(8),parameter::omega3=0.152d0
    real(8),parameter::omega4=0.193d0
    real(8),parameter::omega5=0.251d0
    real(8),parameter::omega6=0.27504d0
    real(8),parameter::omega7=0.604823d0
    real(8),parameter::omega8=0.225d0
    real(8),parameter::omega9=0.344d0


!!!!! foam formation model !!!!!!!!!!!
    real(8),parameter::CMC=0.00014307d0
    real(8),parameter::fmsurf=0.0005d0 
    real(8),parameter::epsurf=1.0d0
    real(8),parameter::FM_foam=1.0d0
    real(8),parameter::fmoil=0.2d0
    real(8),parameter::epoil=1.0d0
    real(8),parameter::epdry=100.0d0
    real(8),parameter::fmdry=0.044d0
    !real(8),parameter::fmmob=33614.0d0
    real(8),parameter::fmmob=10000.0d0
    real(8),parameter::fmHC=0.01d0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!変えた

    real(8),parameter::cogperm=0.3d0!!!!!!!!立花ver
    real(8),parameter::epgperm=3.0d0!!!!!!!!立花ver
    
    real(8),parameter::cofperm=0.3d0!!!!!!!!立花ver
    real(8),parameter::epfperm=3.0d0!!!!!!!!立花ver
    
    real(8),parameter::cooperm=0.4d0!!!!!!!!立花ver
    real(8),parameter::epoperm=3.0d0!!!!!!!!立花ver
    
    real(8),parameter::cowperm=1.0d0!!!!!!!!立花ver
    real(8),parameter::epwperm=3.0d0!!!!!!!!立花ver

    real(8),parameter::epvis=1.0d0!!!!!!!!立花ver

    real(8),parameter::Wtg=0.8d0  !Gas Fraction in the Foam.wt％
    real(8),parameter::Wts=0.2d0  !surfactant Fraction in the water.wt％

    !real(8),parameter::WperS=3.0d0!foam中の水とSDSのmol比

end module mod_initial