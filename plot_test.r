


# Check hyst results 
dat = read.table("tmp")
names(dat) = c("label","time","var","dv_dt","dv_dt_lim","df_dt","f_now")

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
