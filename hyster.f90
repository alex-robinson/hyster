
module hyster 

    use nml 

    implicit none 

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: wp  = sp 

    real(wp), parameter :: MV  = -9999.0_wp 
    real(wp), parameter :: pi = 3.141592653589793_wp

    type hyster_par_class  
        character(len=56) :: label 
        integer  :: ntot 
        real(wp) :: fac
        real(wp) :: df_sign 
        real(wp) :: dv_dt_max
        real(wp) :: df_dt_min
        real(wp) :: df_dt_max 
    end type 

    type hyster_class

        type(hyster_par_class) :: par 

        ! variables 
        integer :: n
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: var(:)
        real(wp) :: dv_dt
        real(wp) :: df_dt
        
    end type 

    private
    public :: sp, dp, wp, pi 
    public :: hyster_class
    public :: hyster_init 
    public :: hyster_calc_rate

contains

    function hyster_init(filename,time,label) result(hyst)

        character(len=*), intent(IN) :: filename 
        real(wp),         intent(IN) :: time 
        character(len=*), intent(IN), optional :: label 

        type(hyster_class)   :: hyst 
        
        integer :: ntot 
        character(len=56) :: par_label 

        par_label = "hyster"
        if (present(label)) par_label = trim(par_label)//"_"//trim(label)

        ! Load parameters 
        call nml_read(filename,trim(par_label),"ntot",hyst%par%ntot)
        call nml_read(filename,trim(par_label),"fac",hyst%par%fac)
        call nml_read(filename,trim(par_label),"df_sign",hyst%par%df_sign)
        call nml_read(filename,trim(par_label),"dv_dt_max",hyst%par%dv_dt_max)
        call nml_read(filename,trim(par_label),"df_dt_min",hyst%par%df_dt_min)
        call nml_read(filename,trim(par_label),"df_dt_max",hyst%par%df_dt_max)
        
        ! Make sure sign is only +1/-1 
        hyst%par%df_sign = sign(1.0_wp,hyst%par%df_sign)

        ! Define label for this hyster object 
        hyst%par%label = "hyster" 
        if (present(label)) hyst%par%label = trim(hyst%par%label)//"_"//trim(label)

        ! (Re)initialize hyster vectors
        if (allocated(hyst%time)) deallocate(hyst%time)
        if (allocated(hyst%var))  deallocate(hyst%var)
        allocate(hyst%time(hyst%par%ntot))
        allocate(hyst%var(hyst%par%ntot))

        ! Initialize variable values
        hyst%time  = MV
        hyst%var   = MV 
        hyst%n     = 0   
        hyst%dv_dt = 0.0_wp 
        hyst%df_dt = 0.0_wp

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
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: var 
        logical,            intent(IN), optional :: verbose 

        ! Local variables 
        real(wp) :: dv_dt(hyst%par%ntot-1)
        real(wp) :: dv_dt_now 
        real(wp) :: f_scale 

        if ( hyst%n .lt. hyst%par%ntot ) then 
            ! Number of timesteps not reached yet, fill in the hyst vectors

            hyst%n = hyst%n + 1 
            hyst%time(hyst%n) = time 
            hyst%var(hyst%n)  = var 

        else 
            ! Keep a running average removing oldest point and adding current one
            hyst%time = eoshift(hyst%time,1,boundary=time)
            hyst%var  = eoshift(hyst%var, 1,boundary=var)
            
            ! Calculate mean rate of change for ntot time steps 
            dv_dt = (hyst%var(2:hyst%par%ntot)-hyst%var(1:hyst%par%ntot-1)) / &
                      (hyst%time(2:hyst%par%ntot)-hyst%time(1:hyst%par%ntot-1))
            hyst%dv_dt = sum(dv_dt) / real(hyst%par%ntot,wp)

            ! Limit the absolute value of dv_dt to the max value threshold for calculating function
            dv_dt_now = min(abs(hyst%dv_dt),hyst%par%dv_dt_max)

            ! Calculate the current df_dt based on allowed values and dv_dt

            ! BASED ON EXPONENTIAL (sharp transition, tuneable)
            f_scale = exp(-hyst%par%fac*dv_dt_now/hyst%par%dv_dt_max)   ! Returns scalar in range [0-1]

            ! Get forcing rate of change in [f/1e6 a]
            hyst%df_dt = hyst%par%df_sign * ( hyst%par%df_dt_min + f_scale*(hyst%par%df_dt_max-hyst%par%df_dt_min) )

            ! Convert [f/1e6 a] => [f/a]
            hyst%df_dt = hyst%df_dt *1e-6 

        end if 

        return

    end subroutine hyster_calc_rate

  
end module hyster 
