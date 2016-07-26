
module hyster 

    use nml 

    implicit none 

    integer,  parameter :: dp  = kind(1.0d0)

    real(dp), parameter :: MV  = -9999.0_dp 

    type hyster_par_class 
        integer  :: ntot 
        real(dp) :: fac
        real(dp) :: f_init, df_sign 
        real(dp) :: dv_dt_max, df_dt_min, df_dt_max 
        
        real(dp) :: v_min, v_max 
        real(dp) :: f_min, f_max 
        
        real(dp) :: expfac
    end type 

    type hyster_class

        type(hyster_par_class) :: par 

        ! variables 
        integer :: n
        real(dp), allocatable :: time(:), v(:)
        real(dp) :: dv_dt, df_dt, f

        
        
    end type 

    private

contains

    function hyster_init(filename,time) result(hyst)

        character(len=*), intent(IN) :: filename 
        real(dp),         intent(IN) :: time 
        type(hyster_class)   :: hyst 
        
        integer :: ntot 

        ! Load parameters 
        call nml_read(filename,"hyster_par","ntot",hyst%par%ntot)
        call nml_read(filename,"hyster_par","fac",hyst%par%fac)
        call nml_read(filename,"hyster_par","f_init",hyst%par%f_init)
        call nml_read(filename,"hyster_par","f_init",hyst%par%f_init)
        call nml_read(filename,"hyster_par","df_sign",hyst%par%df_sign)
        call nml_read(filename,"hyster_par","dv_dt_max",hyst%par%dv_dt_max)
        call nml_read(filename,"hyster_par","df_dt_min",hyst%par%df_dt_min)
        call nml_read(filename,"hyster_par","df_dt_max",hyst%par%df_dt_max)
        
        call nml_read(filename,"hyster_par","v_min",hyst%par%v_min)
        call nml_read(filename,"hyster_par","v_max",hyst%par%v_max)
        call nml_read(filename,"hyster_par","f_min",hyst%par%f_min)
        call nml_read(filename,"hyster_par","f_max",hyst%par%f_max)

        hyst%par%expfac    = dexp(hyst%par%fac)
        hyst%par%df_dt_min = hyst%par%df_dt_min/1d6          ! deg/million years
        hyst%par%df_dt_max = hyst%par%df_dt_max/1d6          ! deg/million years



        ! (Re)initialize hyster vectors
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%v))    deallocate(hyst%v)
        allocate(hyst%time(hyst%par%ntot))
        allocate(hyst%v(hyst%par%ntot))

        ! Initialize variable values
        hyst%time = MV
        hyst%v    = MV 
        hyst%f    = hyst%par%f_init  
        hyst%dv_dt = 0.0_dp 
        hyst%df_dt = 0.0_dp

        return 

    end function hyster_init 

  
!     function hyster_calc_rate(hyst,time,fnow) result(dforc_dt)
!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         ! Subroutine :  d t T r a n s 1
!         ! Author     :  Alex Robinson
!         ! Purpose    :  Generate correct T_warming for gradual changes over
!         !               time (continuous stability diagram!)
!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

!         type(hyster_class) :: hyst 
!         real(dp) :: time, fnow
!         real(dp) :: dvar_dt 

!         if (n_step .eq. 0) then

!             hyst%f   = fnow
!             hyst%n = 0
!             hyst%dVdt = 0.0_dp

!             transT%dVdtm = 0.0_dp
!             transT%mstep = 0

!         else if ( n_step .eq. transT%nstep ) then

!             dTtrans1         = transT%T_warming

!         else

!             transT%mstep = transT%mstep+1
!             if (transT%mstep .gt. size(transT%dVdtm)) transT%mstep = 1
!             transT%dVdtm(transT%mstep) = transT%dVdt

!             ! Get the absolute value of the current change in vol (Gt/a)
!             ! (not greater than the max allowed)
!             !       dVdt_now = min(dabs(transT%dVdt),dVdt_max)
!             dVdtm_now = sum(transT%dVdtm,transT%dVdtm .ne. 0.d0)/max(1,count(transT%dVdtm .ne. 0.d0))
!             dVdt_now  = min(dabs(dVdtm_now),dVdt_max)

!             ! Calculate the current dTdt based on allowed values and rate of smb 
!             ! BASED ON COS (smoother transition)
!             !dTdt = (dTdt_max-dTdt_min)*0.5d0*(dcos(pi*dVdt_now/dVdt_max)+1.d0) + dTdt_min

!             ! BASED ON EXPONENTIAL (sharper transition, tuneable)
!             dTdt = (dTdt_max-dTdt_min)*(1.d0-(dexp(dVdt_now/dVdt_max*fac)-1.d0)/(expfac-1.d0)) + dTdt_min

!             ! Cooling for ice-free init.
!             if (slow_hyst .lt. 0) dTdt = -dTdt


!             ! Make sure warming should be applied...
!             if ( nstep .ge. T_warming_delay ) then
!                 ! Assign new T_warming values
!                 transT%T_warming = transT%T_warming + dTdt*(n_step-transT%nstep)
!             end if

!             ! Output current T_warming value
!             dTtrans1 = transT%T_warming

!             ! Activate kill switch if max/min temperature has been reached
!             if (slow_hyst .gt. 0) then  ! T_warming is increasing
!                 Tlim = T_warming + T_diff
!                 if ( dTtrans1 .gt. Tlim ) kill = 1
!             else                        ! T_warming is decreasing
!                 Tlim = T_warming - T_diff
!                 if ( dTtrans1 .lt. Tlim ) kill = 1
!             end if

!             write(stdout,"(a5,3f10.3,3g15.3)") "dTdt ", n_step*1d-3, dTtrans1, dTdt*1e6, &
!                               transT%dVdt, dVdt_now, dVdtm_now

!             transT%nstep = n_step
!             transT%dTdt  = dTdt*1e6

!         end if

!         return

!     end function hyster_calc_rate

  
end module hyster 
