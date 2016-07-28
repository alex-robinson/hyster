program test

    use hyster 

    implicit none 

    real (8), parameter :: pi = 3.141592653589793d0

    type(hyster_class) :: hyst1

    real(8), allocatable :: time(:), var(:), dv_dt(:)
    integer :: nt, k  

    ! Initialize the hyster object
    hyst1 = hyster_init("Greenland.nml",time=0.d0,label="hyster")

    ! Define input
    nt = 10000
    allocate(time(nt),var(nt),dv_dt(nt))

    time(1)  = 0.0
    var(1)   = 0.0
    dv_dt(1) = (cos(pi*time(1)/1000.d0)+1.d0)/2.d0 * hyst1%par%dv_dt_max*1.1d0

    ! First check increasing forcing
    do k = 2, nt  
        time(k)  = time(k-1) + 1.d0
        dv_dt(k) = (cos(pi*time(k)/1000.d0)+1.d0)/2.d0 * hyst1%par%dv_dt_max*1.1d0 
        var(k)   = var(k-1) + dv_dt(k)
        
        call hyster_calc_rate(hyst1,time=time(k),var=var(k),verbose=.TRUE.)

    end do 


    ! Now check decreasing forcing
    hyst1%par%df_sign = -1.d0 
    do k = nt+1, 2*nt  
        time(k)  = time(k-1) + 1.d0
        dv_dt(k) = (cos(pi*time(k)/1000.d0)+1.d0)/2.d0 * hyst1%par%dv_dt_max*1.1d0 
        var(k)   = var(k-1) + dv_dt(k)
        
        call hyster_calc_rate(hyst1,time=time(k),var=var(k),verbose=.TRUE.)

    end do 









end program 

