


# Check hyst results 
if (TRUE) {

    dat = read.table("test.txt")
    names(dat) = c("label","time","var","dv_dt","df_dt","f_now")

    # Delete initial points (up to ntot)
    ii = c(51:length(dat$time))
    dat = dat[ii,]

    par(mfrow=c(4,1))
    par(plt=c(0.1,0.95,0.1,0.95))

    plot(dat$time,dat$var,type="l",ann=FALSE)
    title(ylab="var")

    plot(dat$time,dat$dv_dt,type="l",ann=FALSE)
    title(ylab="dv_dt")

    plot(dat$time,dat$df_dt,type="l",ann=FALSE)
    title(ylab="df_dt")

    plot(dat$time,dat$f_now,type="l",ann=FALSE)
    title(ylab="f_now")


    quartz()
    plot(dat$dv_dt,dat$df_dt)

}


hystrate = function(df_sign,dv_dt_now,dv_dt_max,dv_dt_min,df_dt_max,df_dt_min,fac) 
{
    expfac = exp(fac) 

    df_dt = (df_sign * (df_dt_max-df_dt_min) *
                         (1.0-(exp(dv_dt_now/dv_dt_max*fac)-1.0) / 
                            (expfac-1.0)) + df_dt_min) #*1e-6

    return(df_dt)
}

hystrate2 = function(df_sign,dv_dt_now,dv_dt_max,dv_dt_min,df_dt_max,df_dt_min,fac) 
{

    f_scale = exp(-fac*dv_dt_now/dv_dt_max)   # [0-1]

    df_dt = df_sign * ( df_dt_min + f_scale*(df_dt_max-df_dt_min) )

    #df_dt = df_dt * 1e-6  # [f/1e6 a] => [f/a]

    return(df_dt)
}

if (FALSE) {
    dv_dt = seq(0,30,by=0.5)
    df_dt_up  = hystrate2( 1,dv_dt,dv_dt_max=20.0,df_dt_max=20.0,df_dt_min=0.0,fac=10.0)
    df_dt_dwn = hystrate2(-1,dv_dt,dv_dt_max=20.0,df_dt_max=20.0,df_dt_min=0.0,fac=10.0)
    xlim = range(dv_dt)
    ylim = range(df_dt_up,df_dt_dwn)
    plot(xlim,ylim,type="n",ann=FALSE)
    abline(h=0,lwd=3,col="grey80")
    lines(dv_dt,df_dt_up,lwd=2)
    lines(dv_dt,df_dt_dwn,lwd=2)
}


