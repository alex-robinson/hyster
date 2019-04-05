library(myr)


# Check hyst results 
if (FALSE) {

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

    dv_dt_now = min(abs(dv_dt_now),dv_dt_max)

    df_dt = (df_sign * (df_dt_max-df_dt_min) *
                         (1.0-(exp(dv_dt_now/dv_dt_max*fac)-1.0) / 
                            (expfac-1.0)) + df_dt_min) #*1e-6

    return(df_dt)
}

hystrate2 = function(dv_dt_now,df_sign,df_dt_min,df_dt_max,dv_dt_scale) 
{

    f_scale = exp(-abs(dv_dt_now)/dv_dt_scale)   # [0-1], 0.6 at dv_dt_now==dv_dt_scale

    df_dt = df_sign * ( df_dt_min + f_scale*(df_dt_max-df_dt_min) )

    #df_dt = df_dt * 1e-6  # [f/1e6 a] => [f/a]

    return(df_dt)
}

if (TRUE) {

    dv_dt_scale = 10.0 
    
    dv_dt = seq(0,100,by=0.5)
    df_dt_up  = hystrate2(dv_dt, 1,df_dt_min=0.0,df_dt_max=200.0,dv_dt_scale=dv_dt_scale)
    df_dt_dwn = hystrate2(dv_dt,-1,df_dt_min=0.0,df_dt_max=200.0,dv_dt_scale=dv_dt_scale)
    
    #xlim = range(dv_dt)
    xlim = c(0,100)
    ylim = range(df_dt_up,df_dt_dwn)

    #ylab  = bquote("dT/dt [K / 10"^"6"*" a"*"]")
    ylab  = "dT/dt [K per million yr]"

    myfigure("./","hyster_rate",asp=1.1,pointsize=20,type="png")

    par(plt=c(0.13,0.95,0.13,0.95))
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    axis(1)
    axis(2)
    mtext(side=1,line=1.5,las=0,"dV/dt [Gt/a]")
    mtext(side=2,line=1.8,las=0,ylab)
    abline(h=0,lwd=5,col="grey80")
    abline(v=dv_dt_scale,lwd=2,col=1)

    lines(dv_dt, df_dt_up,lwd=4)
    lines(dv_dt,df_dt_dwn,lwd=4)

    box()
    graphics.off()
}


