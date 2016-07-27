
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
        
        logical  :: kill 
        real(dp) :: expfac
    end type 

    type hyster_class

        type(hyster_par_class) :: par 

        ! variables 
        integer :: n
        real(dp), allocatable :: time(:), var(:)
        real(dp) :: dv_dt, df_dt, f_now

        
        
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

        ! Additional parameter modification
        hyst%par%expfac    = dexp(hyst%par%fac)
!         hyst%par%df_dt_min = hyst%par%df_dt_min/1d6          ! deg/million years
!         hyst%par%df_dt_max = hyst%par%df_dt_max/1d6          ! deg/million years


        ! (Re)initialize hyster vectors
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%var))  deallocate(hyst%var)
        allocate(hyst%time(hyst%par%ntot))
        allocate(hyst%var(hyst%par%ntot))

        ! Initialize variable values
        hyst%time  = MV
        hyst%var   = MV 
        hyst%n     = 0 
        hyst%f_now = hyst%par%f_init  
        hyst%dv_dt = 0.0_dp 
        hyst%df_dt = 0.0_dp

        ! Initially set kill to false
        ! (to be activated when min/max forcing is reached) 
        hyst%par%kill = .FALSE. 

        return 

    end function hyster_init 

  
    subroutine hyster_calc_rate(hyst,time,var,verbose)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Subroutine :  d t T r a n s 1
        ! Author     :  Alex Robinson
        ! Purpose    :  Generate correct T_warming for gradual changes over
        !               time (continuous stability diagram!)
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

        type(hyster_class), intent(INOUT) :: hyst 
        real(dp),           intent(IN)    :: time
        real(dp),           intent(IN)    :: var 
        logical,            intent(IN), optional :: verbose 

        ! Local variables 
        real(dp) :: dv_dt(hyst%par%ntot-1)

        if ( hyst%n .lt. hyst%par%ntot ) then 
            ! Number of timesteps not reached yet, fill in the hyst vectors

            hyst%n = hyst%n + 1 
            hyst%time(hyst%n) = time 
            hyst%var(hyst%n)  = var 

        end if 

        if ( hyst%n .eq. hyst%par%ntot ) then  
            ! Maximum number of timesteps reached, update the forcing rate of change

            ! Calculate mean rate of change for ntot time steps 
            dv_dt = (hyst%var(2:hyst%par%ntot)-hyst%var(1:hyst%par%ntot-1)) / &
                      (hyst%time(2:hyst%par%ntot)-hyst%time(1:hyst%par%ntot-1))
            hyst%dv_dt = sum(dv_dt) / real(hyst%par%ntot,kind=dp)

            ! Limit the absolute dv_dt to the max value threshold
            hyst%dv_dt = min(dabs(hyst%dv_dt),hyst%par%dv_dt_max)

            ! Calculate the current df_dt based on allowed values and dv_dt 
            ! BASED ON COS (smoother transition)
!             hyst%df_dt = (hyst%par%df_dt_max-hyst%par%df_dt_min) * &
!                        0.5_dp*(dcos(pi*hyst%dv_dt/hyst%par%dv_dt_max)+1.0_dp) + hyst%par%df_dt_min

            ! BASED ON EXPONENTIAL (sharper transition, tuneable)
            hyst%df_dt = (hyst%par%df_dt_max-hyst%par%df_dt_min) * &
                         (1.0_dp-(dexp(hyst%dv_dt/hyst%par%dv_dt_max*hyst%par%fac)-1.0_dp)/(expfac-1.0_dp)) &
                          + hyst%par%df_dt_min

            ! Adjust sign of forcing rate of change
            hyst%df_dt = hyst%df_dt * hyst%par%df_sign

            ! Reset hyst vectors
            hyst%time  = MV
            hyst%var   = MV 
            hyst%n     = 0 

        end if 


        ! Update the forcing every time step 
        hyst%f_now = hyst%f_now + hyst%df_dt * (time-hyst%time(hyst%n-1))


        if (hyst%f_now .le. hyst%par%f_min .or. hyst%f_now .ge. hyst%par%f_max) then 
            ! Activate kill switch if max/min temperature has been reached
            hyst%par%kill = .TRUE. 

        end if 

        if (present(verbose)) then 
            if (verbose) then 
                write(*,"(a5,3f10.3,3g15.3)") "hyster ", &
                time*1d-3, var, dTtrans1, dTdt*1e6, &
                          transT%dVdt, dVdt_now, dVdtm_now
            end if
        end if 

        return

    end subroutine hyster_calc_rate

  
end module hyster 
