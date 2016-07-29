
module hyster 

    use nml 

    implicit none 

    integer,  parameter :: dp  = kind(1.0d0)

    real(dp), parameter :: MV  = -9999.0_dp 
    real(dp), parameter :: pi = 3.141592653589793_dp

    type hyster_par_class 
        logical  :: use_hyster 
        character(len=56) :: label 
        integer  :: ntot 
        character(len=3) :: func
        real(dp) :: fac
        real(dp) :: f_init, df_sign 
        real(dp) :: dv_dt_max, df_dt_min, df_dt_max 
        
        real(dp) :: v_min, v_max 
        real(dp) :: f_min, f_max 
        
        logical  :: active, kill 
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
    public :: hyster_class
    public :: hyster_init 
    public :: hyster_calc_rate

contains

    function hyster_init(filename,time,label) result(hyst)

        character(len=*), intent(IN) :: filename 
        real(dp),         intent(IN) :: time 
        character(len=*), intent(IN), optional :: label 

        type(hyster_class)   :: hyst 
        
        integer :: ntot 

        ! Load parameters 
        call nml_read(filename,"hyster_par","use_hyster",hyst%par%use_hyster)
        call nml_read(filename,"hyster_par","ntot",hyst%par%ntot)
        call nml_read(filename,"hyster_par","func",hyst%par%func)
        call nml_read(filename,"hyster_par","fac",hyst%par%fac)
        call nml_read(filename,"hyster_par","f_init",hyst%par%f_init)
        call nml_read(filename,"hyster_par","df_sign",hyst%par%df_sign)
        call nml_read(filename,"hyster_par","dv_dt_max",hyst%par%dv_dt_max)
        call nml_read(filename,"hyster_par","df_dt_min",hyst%par%df_dt_min)
        call nml_read(filename,"hyster_par","df_dt_max",hyst%par%df_dt_max)
        
!         call nml_read(filename,"hyster_par","v_min",hyst%par%v_min)
!         call nml_read(filename,"hyster_par","v_max",hyst%par%v_max)
        call nml_read(filename,"hyster_par","f_min",hyst%par%f_min)
        call nml_read(filename,"hyster_par","f_max",hyst%par%f_max)

        ! Additional parameter modification
        hyst%par%expfac    = dexp(hyst%par%fac)
!         hyst%par%df_dt_min = hyst%par%df_dt_min/1d6          ! deg/million years
!         hyst%par%df_dt_max = hyst%par%df_dt_max/1d6          ! deg/million years
    
        ! Make sure sign is only +1/-1 
        hyst%par%df_sign = sign(1.0_dp,hyst%par%df_sign)

        ! Define label for this hyster object 
        hyst%par%label = "hyster" 
        if (present(label)) hyst%par%label = trim(label)

        ! (Re)initialize hyster vectors
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%var))  deallocate(hyst%var)
        allocate(hyst%time(hyst%par%ntot))
        allocate(hyst%var(hyst%par%ntot))

        ! Initialize variable values
        hyst%time  = MV
        hyst%var   = MV 
        hyst%n     = 0   
        hyst%dv_dt = 0.0_dp 
        hyst%df_dt = 0.0_dp

        hyst%f_now = hyst%par%f_init
        if (hyst%par%f_init .eq. MV) then 
            if (hyst%par%df_sign .lt. 0.0_dp) then 
                hyst%f_now = hyst%par%f_max 
            else 
                hyst%f_now = hyst%par%f_min 
            end if 
        end if 

        ! Initially set kill to false
        ! (to be activated when min/max forcing is reached) 
        hyst%par%kill = .FALSE. 

        ! Intially set active to false, to populate array first
        hyst%par%active = .FALSE. 

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
        real(dp) :: dv_dt_now 

        if ( hyst%n .lt. hyst%par%ntot ) then 
            ! Number of timesteps not reached yet, fill in the hyst vectors

            hyst%n = hyst%n + 1 
            hyst%time(hyst%n) = time 
            hyst%var(hyst%n)  = var 

        else 
            ! Keep a running average removing oldest point and adding current one
            hyst%time = eoshift(hyst%time,1,boundary=time)
            hyst%var  = eoshift(hyst%var, 1,boundary=var)
            
!         if ( hyst%n .eq. hyst%par%ntot ) then  
!             ! Maximum number of timesteps reached, update the forcing rate of change

            ! Calculate mean rate of change for ntot time steps 
            dv_dt = (hyst%var(2:hyst%par%ntot)-hyst%var(1:hyst%par%ntot-1)) / &
                      (hyst%time(2:hyst%par%ntot)-hyst%time(1:hyst%par%ntot-1))
            hyst%dv_dt = sum(dv_dt) / real(hyst%par%ntot,kind=dp)

            ! Limit the absolute dv_dt to the max value threshold for calculating function
            dv_dt_now = min(dabs(hyst%dv_dt),hyst%par%dv_dt_max)

            ! Calculate the current df_dt based on allowed values and dv_dt
            select case(trim(hyst%par%func))

                case("cos")
                    ! BASED ON COS (smoother transition)
                    hyst%df_dt = (hyst%par%df_sign * (hyst%par%df_dt_max-hyst%par%df_dt_min) * &
                       0.5_dp*(dcos(pi*dv_dt_now/hyst%par%dv_dt_max)+1.0_dp) &
                       + hyst%par%df_dt_min) *1d-6 ! [f/1e6 a] => [f/a]

                case("exp")
                    ! BASED ON EXPONENTIAL (sharper transition, tuneable)
                    hyst%df_dt = (hyst%par%df_sign * (hyst%par%df_dt_max-hyst%par%df_dt_min) * &
                         (1.0_dp-(dexp(dv_dt_now/hyst%par%dv_dt_max*hyst%par%fac)-1.0_dp) / &
                            (hyst%par%expfac-1.0_dp)) + hyst%par%df_dt_min) *1d-6  ! [f/1e6 a] => [f/a]

                case DEFAULT 

                    write(*,*) "hyster_calc_rate:: error: function not recognized (cos,exp): ", &
                                trim(hyst%par%func)
                    write(*,*) "hyster_label = ", trim(hyst%par%label)
                    stop 

            end select 

!             ! Reset hyst vectors
!             hyst%time  = MV
!             hyst%var   = MV 
!             hyst%n     = 0 

        end if 


        ! Update the forcing every time step 
        hyst%f_now = hyst%f_now + hyst%df_dt * (time-hyst%time(hyst%n-1))


        ! When forcing reaches limits, stop changing forcing 
        ! and activate kill switch 
        if (hyst%f_now .le. hyst%par%f_min .and. hyst%par%df_sign .lt. 0.0) then 
            hyst%f_now = hyst%par%f_min
            hyst%par%kill = .TRUE. 
        else if (hyst%f_now .ge. hyst%par%f_max .and. hyst%par%df_sign .gt. 0.0) then 
            hyst%f_now = hyst%par%f_max
            hyst%par%kill = .TRUE. 
        end if 

        if (present(verbose)) then 
            if (verbose) then 
                write(*,"(a,1x,1f10.3,5g15.3)") trim(hyst%par%label), &
                time*1d-3, var, hyst%dv_dt, dv_dt_now, hyst%df_dt*1e6, hyst%f_now
            end if
        end if 

        return

    end subroutine hyster_calc_rate

  
end module hyster 
