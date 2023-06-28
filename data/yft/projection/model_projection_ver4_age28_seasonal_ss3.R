library(ss3om)
library(FLasher)
library(patchwork)
library(FLash)
library(FLBRP)
library(ggplot2)
#stock synthesis release  https://github.com/nmfs-stock-synthesis/stock-synthesis/releases
#download in same folders ast ss3 files (four input files ) from Dan FU - click on exe takes a while to run
# --- YFT

#which you can also load from ss3
# LOAD SS_output
# - Using own ss3om function to deal with compressed files

out        = readOutputss3('C:\\alexIRD\\iotc_tuna_fads\\data\\yft\\')

# build stock object from r4ss::SS_output(),
# BUILD FLFisheries + FLBiol
# - Returned as list with elements 'biol' and 'fisheries'
oms        = buildFLBFss330(out)

# EXTRACT biology and fisheries
#The first conatins biology (n,wt, mat, fec, SRR) 
bio        = oms$biol

#FLFishery objects,
#basically fleets with no metiers. Each fleet has effort and capacity
#time series, and then a list of FL:catch objects, one for each stock
#they catch:
fis        = oms$fisheries

#You can create an FLStock for the two objects
stk        = as.FLStock(bio, fis)[,ac(100:296),,,]

#creates a stock object and combines areas into 1, There is also noseason and nounit.simplify(stk, 'area') 
yft_stock  = noarea(stk) 


#Annual rates converted to quarterly
m(yft_stock)[2:29,,,,]       = yft_stock@m[2:29,,,,]*0.25
harvest(yft_stock)[1:29,,,,] = yft_stock@harvest[1:29,,,,]*0.25

#rename fisheries
names(fis) = c('(GI 1a)','(HD 1a)','(LL 1a)','(OT 1a)','(BB 1b)','(FS 1b)','(LL 1b)','(LS 1b)','(TR 1b)','(LL 2)','(LL 3)','(GI 4)','(LL 4)',	'(OT 4)','(TR 4)','(FS 2)','(LS 2)','(TR 2)','(FS 4)','(LS 4)','(FL 4)')

#obtain DFAD catch numbers at age year 2020
DFADcatch.n      = catch.n(fis[['(LS 1b)']])[,ac(296),,,]+catch.n(fis[['(LS 2)']])[,ac(296),,,]+catch.n(fis[['(LS 4)']])[,ac(296),,,]

#calculate revised catch numbers for 100% increase/decrease, 1% FAD sets = 0.29% change in catch weight
##calculate revised catch numbers for 50% increase/decrease, 1% FAD sets = 0.29% change in catch weight
newcatchwt100    = apply(DFADcatch.n*yft_stock@stock.wt[,ac(296),,,,],2,sum)*1.29 #100%
newcatchwt50     = apply(DFADcatch.n*yft_stock@stock.wt[,ac(296),,,,],2,sum)*0.855 #50%

#calculate DFAD proportion by weight
prop_by_weight   = DFADcatch.n*yft_stock@stock.wt[,ac(296),,,,]/c(apply(DFADcatch.n*yft_stock@stock.wt[,ac(296),,,,],2,sum))

#weight at age
newweight_age100 = prop_by_weight*c(newcatchwt100)
newweight_age50  = prop_by_weight*c(newcatchwt50)

#revised catch numbers for calculating F
#100% increase
rev_catch.n100 = unlist((newweight_age100/yft_stock@stock.wt[,ac(296),,,,]) + (yft_stock@catch.n[,ac(296),,,]- DFADcatch.n ))
#50% reduction
rev_catch.n50  = unlist((yft_stock@catch.n[,ac(296),,,]- DFADcatch.n) + unlist(newweight_age50/yft_stock@stock.wt[,ac(296),,,]))
#no FADS
rev_catch.n0   = unlist(yft_stock@catch.n[,ac(296),,,]-DFADcatch.n)

#calculate new F function
newparam       = list()
type           = function(control=yft_stock@catch.n[d,ac(296),,,,])   {
    for (d in (1:29)) {
    x          = as.numeric(yft_stock@harvest[d,ac(296),,,])
    stock      = as.numeric(yft_stock@stock.n[d,ac(296),,,])
    m          = as.numeric(yft_stock@m[d,ac(296),,,])
    catch      = as.numeric(control[d])  #new catch to be put in
    
    fr         = function(x) ((x/(x+m)*stock*(1-exp(-(x+m)))-catch)^2) #^2 minimises the negative number
    opt.ct     = optim(x,fr, method = "L-BFGS-B")#, hessian = TRUE, lower=0, upper=2.0)
 newparam[[d]] = round(as.numeric(opt.ct[1]),5)
    }
    newparam    = unlist(newparam)    
    yft_stock@harvest[,ac(296),,,] = FLQuant(c(newparam),  dimnames=list(0:28, year=ac(296)))  
    return(yft_stock@harvest[,ac(296),,,])
} 

#Here the F's are changed
bau                    = type(control=yft_stock@catch.n[,ac(296),,,,])
inc100                 = type(control=rev_catch.n100 )
dec50                  = type(control=rev_catch.n50)
fob0                   = type(control=rev_catch.n0)

yftsr                  = as.FLSR(yft_stock)

model(yftsr)           = bevholt();
#yftsr                     <-transform(yftsr,rec=rec,ssb=ssb)

yftsr                  = fmle(yftsr);
#yftsr                  = fmle(yftsr,fixed=list(s=0.8));

status_quo             = mean(bau)
inc100x                = mean(inc100)
dec50x                 = mean(dec50)
fob0x                  = mean(fob0)

ctrl_bau               = data.frame(year = 296:336, quantity = "f", val = status_quo)
ctrl_100               = data.frame(year = 296:336, quantity = "f", val = inc100x)
ctrl_50                = data.frame(year = 296:336, quantity = "f", val = dec50x )
ctrl_0                 = data.frame(year = 296:336, quantity = "f", val = fob0x )
ctrl_fbau              = fwdControl(ctrl_bau)
ctrl_f100              = fwdControl(ctrl_100)
ctrl_f50               = fwdControl(ctrl_50)
ctrl_f0                = fwdControl(ctrl_0)

#set the data up for 10 year projection
yft_stf_bau            = stf(yft_stock, nyears=40)
yft_stf_100            = stf(yft_stock, nyears=40)
yft_stf_50             = stf(yft_stock, nyears=40)
yft_stf_0              = stf(yft_stock, nyears=40)

yft_f_sq               = fwd(yft_stf_bau, ctrl = ctrl_fbau, sr = yftsr)
yft_f_100              = fwd(yft_stf_100, ctrl = ctrl_f100, sr = yftsr)
yft_f_50               = fwd(yft_stf_50,  ctrl = ctrl_f50,  sr = yftsr)
yft_f_0                = fwd(yft_stf_0,   ctrl = ctrl_f0,   sr = yftsr)

plot(window(yft_f_sq,  start = 296))
plot(window(yft_f_100, start = 296))
plot(window(yft_f_50,  start = 296))
plot(window(yft_f_0,   start = 296))
plot(FLStocks(yft_f_sq[,ac(296:336),,,],yft_f_100[,ac(296:336),,,],yft_f_50[,ac(296:336),,,],yft_f_0[,ac(296:336),,,]))
test=subset(as.data.frame(window(yft_f_sq,  start = 296) ), slot %in%c("catch","harvest"))
