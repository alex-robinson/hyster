program test

    use hyster 

    implicit none 

    type(hyster_class) :: hyst1

    real(wp), allocatable :: time(:), var(:), dv_dt(:), for(:) 
    real(wp) :: dt 
    integer  :: nt, k  
    
    ! Initialize the hyster object
    hyst1 = hyster_init("Greenland.nml",time=0.0,label="test1")

    ! Define input
    nt = 10000
    allocate(time(nt*2),var(nt*2),dv_dt(nt*2),for(nt*2))

    time(1)  = 0.0
    var(1)   = 0.0
    dv_dt(1) = (cos(pi*time(1)/1000.0)+1.0)/2.0 * hyst1%par%dv_dt_max*1.1
    for(1)   = -10.0   ! [K] Initial forcing value 

    dt = 1.0 

    ! First check increasing forcing
    do k = 2, nt  
        time(k)  = time(1) + (k-1)*dt
        dv_dt(k) = (cos(pi*time(k)/1000.0)+1.0)/2.0 * hyst1%par%dv_dt_max*1.1
        var(k)   = var(k-1) + dv_dt(k)*dt
        
        call hyster_calc_rate(hyst1,time=time(k),var=var(k),verbose=.TRUE.)

        for(k)   = for(k-1) + hyst1%df_dt*dt 
        
        write(*,"(a,1x,1f10.3,5g15.3)") trim(hyst1%par%label), &
                time(k)*1e-3, var(k), dv_dt(k), hyst1%df_dt*1e6, for(k)

    end do 


    ! Now check decreasing forcing

    hyst1%par%df_sign = -1.0 

    do k = nt+1, 2*nt  
        time(k)  = time(1) + (k-1)*dt
        dv_dt(k) = (cos(pi*time(k)/1000.0)+1.0)/2.0 * hyst1%par%dv_dt_max*1.1
        var(k)   = var(k-1) + dv_dt(k)*dt
        
        call hyster_calc_rate(hyst1,time=time(k),var=var(k),verbose=.TRUE.)

        for(k)   = for(k-1) + hyst1%df_dt*dt 

        write(*,"(a,1x,1f10.3,5g15.3)") trim(hyst1%par%label), &
                time(k)*1e-3, var(k), dv_dt(k), hyst1%df_dt*1e6, for(k)

    end do 



end program 

