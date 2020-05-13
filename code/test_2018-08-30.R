## https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html
install.packages('survminer')
source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA.clinical") # data for examples

library(survminer)
library(survival)
require( readxl)
#require( caret )
require( rms)

#### install.packages("readxl")
#setwd("G:/Projects/Kim_YC/Yang, George")
#dd <- read_excel("./data/FLIPI_prognosticscore excel.xlsx")
dd <- read_excel("./data/FLIPI_Lymphopenia_clean_Youngchul.xlsx")
  head(dd)

dd$OS  <- dd$RECALCULATEDSurvivalTimeMonths
dd$death <- dd$VS
  table( dd$FLIPI)
  table( dd$FLIPIALS)
dd$FLIPI <- factor(dd$FLIPI)
dd$FLIPIALS <- factor(dd$FLIPIALS)

colfunc <-colorRampPalette(c("red", "green"))
plot(rep(1,50),col=(colfunc(50)), pch=19,cex=2)
  colfunc(6)


table( dd$FLIPI, dd$FLIPIALS)
table( dd$FLIPI, dd$AbsLymphScore1lt10gt1)
  

### Graphical Calibration of survival probability at 120 months 
par( mfrow=c(1,2))  
#### FLIPI score
sdata <- dd
model.form <- formula("Surv(OS, death) ~  factor(FLIPI)")

#fit <- survfit( model.form, data=sdata )
#fit <- survfit( Surv(OS, death) ~  FLIPI, data=sdata )
time <- 120
fit <- coxph(Surv(OS, death)~FLIPIALS, data = sdata)
bh <- basehaz(fit)
#survdiff( model.form,data=sdata )
## get the culative hazard at the time of 120
H.120 <- bh[ which(bh$time==time), ]$hazard
xb <- predict(fit, type="lp")
coxph.surv.prob.120 <- (exp(-1*H.120))^exp(xb)
#plot( dd$OS, surv.prob.120)


### Define cutpoints to form prognostic groups: 6 groups.
cutpoints <- sort(unique(coxph.surv.prob.120))
sdata$group <- cut(coxph.surv.prob.120, breaks=c(0, cutpoints[1:7]), include.lowest = F, right=T) #, labels=1:length(cutpoints))
  table(sdata$group)
levels(sdata$group) <- (nlevels(sdata$group)-1):0


## Average survival probability at 120 months based on Cox model by group  
p.coxph <- tapply(coxph.surv.prob.120, sdata$group, mean)
  p.coxph

  
## KM survival probability    
##ss <- survfit(Surv(OS, death)~group, data=subset(dd, group==1))  
km.fit <- survfit(Surv(OS, death)~group, data=sdata)
s <- summary(km.fit, times=120, extend=TRUE) #$surv
  s$surv

table(FLIPIALS=sdata$FLIPIALS, PI.group=sdata$group)
  
## Compare Average survival probability of Cox model with KM-based survival probability
plot(p.coxph, s$surv, xlim=c(0,1), ylim=c(0,1), pch=4, type="b", lwd=3,
     xlab="Predicted Survival Probability at 120 months", 
     ylab="Kaplan-Meier Survival Probability at 120 months")  
abline(0,1, lty=2)
 
 
## Discrepancy calculated on the complete dataset.
d <- p.coxph - s$surv
  d
d.observed <-  d
p.coxph.whole <- p.coxph



### u=The time point for which to validate predictions for survival models. 
## For cph fits, you must have specified surv=TRUE, time.inc=u, where u is the constant specifying the time to predict.
par( mfrow=c(2,2))
f <- cph(Surv(OS, death)~FLIPIALS, data=dd, x=TRUE, y=TRUE, time.inc=120, surv=TRUE)
cal <- calibrate(f, u=120, B=200)
plot(cal, xlim=c(0, 1), ylim=c(0,1))
calkm <- calibrate(f, u=120, cmethod='KM', B=200, 
                   cuts=c(0, cutpoints, 1))
plot(calkm, add=T)

validate(f, B=200)
rcorr.cens(as.numeric(dd$FLIPIALS), Surv(dd$OS, dd$death))


## FLIPI model
f <- cph(Surv(OS, death)~FLIPI, data=dd, x=TRUE, y=TRUE, time.inc=120, surv=TRUE)
cal <- calibrate(f, u=120, B=200)
plot(cal, xlim=c(0, 1), ylim=c(0,1))

calkm <- calibrate(f, u=120, cmethod='KM', B=200, 
                   cuts=c(0, cutpoints, 1))
plot(calkm, add=T)

validate(f, B=200)
rcorr.cens(as.numeric(dd$FLIPIALS), Surv(dd$OS, dd$death))





#### Bootstrapping-based Optimism for a model calibration of FLIPI + ALS prognostic index #####
set.seed(90)
d.result <- NULL
km.prob <- NULL

for( B in 1:200){
  sdata <- dd[sample(1:nrow(dd), nrow(dd), replace = T), ]
  fit <- coxph(Surv(OS, death)~FLIPIALS, data = sdata)
  #survdiff( model.form,data=sdata )
  ## get the culative hazard at the time of 120
  H.120 <- bh[ which(bh$time==time), ]$hazard
  xb <- predict(fit, type="lp")
  coxph.surv.prob.120 <- (exp(-1*H.120))^exp(xb)
  
  #plot( dd$OS, surv.prob.120)
  # Use cutpoints defined in the whole data
  sdata$group <- cut(coxph.surv.prob.120, breaks=c(0, cutpoints[2:6], 1), include.lowest=T)
    table(sdata$group)
  levels(sdata$group) <- 1:nlevels(sdata$group)
    
  ## average of Cox-model deriven survival probabilities for each interval
  p.coxph <- tapply(coxph.surv.prob.120, sdata$group, mean)
  
  if( all(!is.na(p.coxph))){
    ##ss <- survfit(Surv(OS, death)~group, data=subset(dd, group==1))  
    km.fit <- survfit(Surv(OS, death)~group, data=sdata)
    s <- summary(km.fit, times=120, extend=TRUE) #$surv
    ## difference in survival probabilites between cox model and KM on bootstrap samples
    d <- p.coxph - s$surv
    d.result <- rbind(d.result, d)
    km.prob <- rbind(km.prob, s$surv)
    #points(p.coxph.whole, s$surv)
  }
}

#### Summary statistics of KM
m <- colMeans(km.prob)
se <- apply(km.prob, 2, sd)
km.max <- apply(km.prob, 2, max)
km.min <- apply(km.prob, 2, min)

### Here p.coxph.whole is a survival probability that is estimated from Cox model on the complete dataset
### m and se are the mean and standard deviation of Kaplan-Meier-based survival probabilities 
points(p.coxph.whole, m, col="red",  type="b", pch=16)
segments(x0=p.coxph.whole, x1=p.coxph.whole, y0=m, y1=km.max, lwd=2)
segments(x0=p.coxph.whole, x1=p.coxph.whole, y0=m, y1=km.min, lwd=2)

#segments(x0=p.coxph.whole, x1=p.coxph.whole, y0=m, y1=m+1.96*se, lwd=2)
#segments(x0=p.coxph.whole, x1=p.coxph.whole, y0=m, y1=m-1.96*se, lwd=2)

d - colMeans(d.result)
d




####################################################
f <- cph(Surv(OS, death)~FLIPIALS, data=dd, x=TRUE, y=TRUE, time.inc=12, surv=TRUE)
  validate(f, B=200)

v <- validate(f, B=5, method="crossvalidation")
  v
  
  
set.seed(200)
rs <- rs2 <- NULL ## result set
for(B in 1:200){
  bdata <- dd[sample(x=1:nrow(dd), size=nrow(dd), replace = T), ]
  s <- rcorr.cens(-1*as.numeric(bdata$FLIPIALS), Surv(bdata$OS, bdata$death))
  rs <- rbind(rs, s)
  
  s <- rcorr.cens(-1*as.numeric(bdata$FLIPI), Surv(bdata$OS, bdata$death))
  rs2 <- rbind(rs2, s)
}
## Dxy
quantile(rs[,2], probs=c(0.025, 0.975))
quantile(rs2[,2], probs=c(0.025, 0.975))

## C-index
quantile(rs[,1], probs=c(0.025, 0.975))
quantile(rs2[,1], probs=c(0.025, 0.975))



f <- cph(Surv(OS, death)~FLIPI, data=dd, x=TRUE, y=TRUE, time.inc=12, surv=TRUE)
validate(f, B=200)
validate(f, B=5, method="crossvalidation")




####### Cross-validated Kaplan-Meier #######
## step 1. Split a complete data into k-folds. 
## step 2. Perform cross-validation and save prediction scores of a test set in each test-set fold.
## step 3. Generate AUC of the saved prediction score against Overall survival of all data

# https://www.rdocumentation.org/packages/bmrm/versions/3.7/topics/balanced.cv.fold
require( bmrm)
require(survivalROC)
## k=4 to have at least two subjects in each FLIPIALS category
set.seed(200)
k = 4
folds <- balanced.cv.fold(y=dd$FLIPIALS, num.cv=k)
table( dd$FLIPIALS, folds)

time <- 60

par(mfrow=c(1,2))  
  tmp <- dd
  tmp$prob <- NA ## survival probability based on the new model 
  tmp$prob2 <- NA ## survival probability based on the old model
  tmp$group <- NA
  tmp$group2 <- NA
  
  for(j in 1:k){
    # j <- 1
    index.test <- which(folds==j)
    train <- tmp[-index.test, ]
    test <- tmp[index.test, ]
    
    ## FLIPI + ALS Model
    fit.train <- coxph(Surv(OS, death)~FLIPIALS, data=train)
    
    ## cutoff points
    xb <- predict(fit.train, type="lp")
    cpoints <- sort(unique(xb))
      cpoints    
      
    # Prediction on test set : predict(fit.train, newdata=test, type="lp")  
    bh <- basehaz(fit.train)
    H.120 <- bh[ which(bh$time==time), ]$hazard
    xb <- predict(fit.train, newdata=test, type="lp")
    surv.prob.120 <- (exp(-1*H.120))^exp(xb)
    tmp$prob[index.test] <- surv.prob.120
    # tmp$group[index.test] <- as.numeric( cut(xb, breaks=c(-Inf, cpoints, Inf), include.lowest=F, right =T) )
    # tmp$group[index.test] <- match( round(xb,5), round(cpoints,5))
    tmp$group[index.test] <- newgroup <- as.numeric( cut(surv.prob.120, breaks=seq(0, 1, 1/6), include.lowest=F, right =T) )
    table(newgroup)    
    
    table(tmp$group)    
     
    ### FLIPI model
    fit.train <- coxph(Surv(OS, death)~FLIPI, data=train)
    cpoints <- sort(unique(predict(fit.train, type="lp")))
    
    # prediction on test set: predict(fit.train, newdata=test, type="lp")  
    bh <- basehaz(fit.train)
    H.120 <- bh[ which(bh$time==time), ]$hazard
    xb <- predict(fit.train, newdata=test, type="lp")
    surv.prob.120 <- (exp(-1*H.120))^exp(xb)
    
    tmp$prob2[index.test] <- surv.prob.120
    #tmp$group2[index.test] <- as.numeric(cut(xb, breaks=c(-Inf, cpoints, Inf), include.lowest=F, right =T))
    #tmp$group2[index.test] <- match(round(xb,5), round(cpoints,5) )
    tmp$group2[index.test] <- newgroup <- as.numeric( cut(surv.prob.120, breaks=seq(0, 1, 1/6), include.lowest=F, right =T) )
    table(newgroup)    
    
      table(tmp$group2)
  }


  colfunc <-colorRampPalette(c("red", "green"))
  par(mfrow=c(1,1))
  ### FLIPI + ALS
  #plot( fit <- survfit( Surv(OS, death) ~ group, data=tmp), lwd=3, col=colfunc(8), xlab="OS", ylab="Survival Probability" )
  
  tmp$gg <- factor(as.numeric( cut(tmp$prob, breaks=seq(0, 1, 1/6), include.lowest=F, right =T) ))
  plot( fit <- survfit( Surv(OS, death) ~ tmp$gg, data=tmp), lwd=3, col=colfunc(8), xlab="Time", ylab="Survival Probability" )

  x <- tapply(tmp$OS, INDEX=tmp$gg, FUN = max)
  s <- summary(fit)
  y <- tapply(s$surv, s$strata, FUN=min)
  text(x+5, y, names(x) )
  title("cross-validated KM: FLIP+ALS")
  sdiff <- survdiff( Surv(OS, death) ~ gg, data=tmp)
  sdiff
  table(tmp$gg)
  cox.1 <- summary(coxph( Surv(OS, death)~factor(tmp$gg), data=tmp))$coef
  ## plot( fit <- survfit( Surv(OS, death) ~ FLIPIALS, data=dd), lwd=3, col=colfunc(8), xlab="OS", ylab="Survival Probability" )
  
  ggsurvplot(fit, title="Cross-validated KM Survival Curve: FLIP+ALS Model", 
             legend.title="Group", 
             legend.labs=seq(2,6,1), risk.table=TRUE)
  
  
  
  ### FLIPI 
  #plot( fit <- survfit( Surv(OS, death) ~ group2, data=tmp), lwd=3, col=colfunc(8), xlab="OS", ylab="Survival Probability" )
  tmp$gg <- as.numeric( cut(tmp$prob2, breaks=seq(0, 1, 1/6), include.lowest=F, right =T) )
  plot( fit <- survfit( Surv(OS, death) ~ tmp$gg, data=tmp), lwd=3, col=colfunc(8), xlab="OS", ylab="Survival Probability" )

  x <- tapply(tmp$OS, INDEX=tmp$gg, FUN = max)
  s <- summary(fit)
  y <- tapply(s$surv, s$strata, FUN=min)
  
  text(x+5, y, names(x) )
  title("cross-validated KM: FLIPI")  
  sdiff <- survdiff( Surv(OS, death) ~ gg, data=tmp)
  sdiff  
  cox.2 <- summary(coxph( Surv(OS, death)~factor(tmp$gg), data=tmp))$coef
  # P< 2e-16
    
  summary(coxph( model.form, data=sdata ))$coef
  s <- summary(fit)
    
  


####### Cross-validated AUC #######
## step 1. Split a complete data into k-folds. 
## step 2. Perform cross-validation and save prediction scores of a test set in each test-set fold.
## step 3. Generate AUC of the saved prediction score against Overall survival of all data

# https://www.rdocumentation.org/packages/bmrm/versions/3.7/topics/balanced.cv.fold
require( bmrm)
require(survivalROC)
## k=4 to have at least two subjects in each FLIPIALS category
set.seed(200)
k = 4
folds <- balanced.cv.fold(y=dd$FLIPIALS, num.cv=k)
  table( dd$FLIPIALS, folds)

par(mfrow=c(1,2))  
for(time in c(60, 120)){
  tmp <- dd
  tmp$prob <- NA ## survival probability based on the new model 
  tmp$prob2 <- NA ## survival probability based on the old model
  
  for(j in 1:k){
    ## FLIPI + ALS Model
    # j <- 1
    index.test <- which(folds==j)
    train <- tmp[-index.test, ]
    test <- tmp[index.test, ]
    fit.train <- coxph(Surv(OS, death)~FLIPIALS, data=train)
    
    # Prediction on test set : predict(fit.train, newdata=test, type="lp")  
    bh <- basehaz(fit.train)
    H.120 <- bh[ which(bh$time==time), ]$hazard
    xb <- predict(fit.train, newdata=test, type="lp")
    surv.prob.120 <- (exp(-1*H.120))^exp(xb)
    tmp$prob[index.test] <- surv.prob.120
    
    ### FLIPI model
    fit.train <- coxph(Surv(OS, death)~FLIPI, data=train)
    # prediction on test set: predict(fit.train, newdata=test, type="lp")  
    bh <- basehaz(fit.train)
    H.120 <- bh[ which(bh$time==time), ]$hazard
    xb <- predict(fit.train, newdata=test, type="lp")
    surv.prob.120 <- (exp(-1*H.120))^exp(xb)
    tmp$prob2[index.test] <- surv.prob.120
  }
  
  nobs <- nrow(tmp)
  xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob, 
              predict.time=time, method="NNE", span=0.25*nobs^(-0.2))
  AUC.old <- xx$AUC
  plot(xx$FP, xx$TP, type="l", xlim=c(0, 1), ylim=c(0,1), 
       xlab=paste( "False Positive"), 
       ylab="True Positive", lwd=3, 
        main=paste("4-Fold Cross-validation,\n Month = ", time))
  abline(0, 1, lty=2)
  text(0.6, 0.2, paste("AUC of FLIPI+ALS=",round(xx$AUC,3)) )
    
  xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob2, 
                    predict.time=time, method="NNE", span=0.25*nobs^(-0.2))
  points(xx$FP, xx$TP, type="l",col="gray", lwd=3)
  AUC.new <- xx$AUC
  text(0.6, 0.1, paste("AUC of FLIPI=",round(xx$AUC,3)) )
}



## Prob: Prognostic index of FLIPI + ALS
plot(survfit(Surv(OS, death) ~ cut(tmp$prob, breaks=2), data=tmp))
survdiff(Surv(OS, death) ~ cut(tmp$prob, breaks=2), data=tmp)

## Prob: Prognostic index of FLIPI only
plot(survfit(Surv(OS, death) ~ cut(tmp$prob2, breaks=2), data=tmp))
survdiff(Surv(OS, death) ~ cut(tmp$prob2, breaks=2), data=tmp)






#### Confidence interval of cross-validated AUC
set.seed(200)

k = 4
B <- 200
boot.AUC <- matrix(NA, nrow=B, ncol=3)

for(b in 1:B){## bootstrap 
  
  dd.auc <- as.data.frame( dd[sample(1:nrow(dd), nrow(dd), replace = T), ] )
  folds <- balanced.cv.fold(y=dd.auc$FLIPIALS, num.cv=k)
  table( dd.auc$FLIPIALS, folds)
  
  par(mfrow=c(1,2))  
  for(time in c(60)){
    tmp <- dd.auc
    tmp$prob <- NA ## survival probability based on the new model 
    tmp$prob2 <- NA ## survival probability based on the old model
    
    for(j in 1:k){
      ## FLIPI + ALS Model
      # j <- 4
      index.test <- which(folds==j)
      train <- tmp[-index.test, ]
      test <- tmp[index.test, ]
      
      
      if( exists("w")) rm(w)  
      w <- tryCatch(
        fit.train <- coxph(Surv(OS, death)~FLIPIALS, data=train, iter.max=100),
        error=function(e) e, warning=function(w) w
      )  
      w <- ifelse( is(w, "warning"), 1, 0 ) 
                    
      #fit.train <- coxphf(Surv(OS, death)~FLIPIALS, data=train, maxit=1000)
      
      # Prediction on test set : predict(fit.train, newdata=test, type="lp")  
      bh <- basehaz(fit.train)
      H.120 <- bh[ which.min(abs(bh$time-time)), ]$hazard
      xb <- predict(fit.train, newdata=test, type="lp")
      surv.prob.120 <- (exp(-1*H.120))^exp(xb)
      tmp$prob[index.test] <- surv.prob.120
      
      ### FLIPI model
      fit.train <- coxph(Surv(OS, death)~FLIPI, data=train)
      #fit.train <- coxphf(Surv(OS, death)~FLIPI, data=train, maxit=1000)
      
      # prediction on test set: predict(fit.train, newdata=test, type="lp")  
      bh <- basehaz(fit.train)
      H.120 <- bh[ which.min(abs(bh$time-time)), ]$hazard
      xb <- predict(fit.train, newdata=test, type="lp")
      surv.prob.120 <- (exp(-1*H.120))^exp(xb)
      tmp$prob2[index.test] <- surv.prob.120
    }
    
    nobs <- nrow(tmp)
    xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob, 
                      predict.time=time, method="NNE", span=0.25*nobs^(-0.2))
    AUC.new <- xx$AUC
    #plot(xx$FP, xx$TP, type="l", xlim=c(0, 1), ylim=c(0,1), 
    #     xlab=paste( "FP", "", "AUC = ",round(xx$AUC,3)), 
    #     ylab="TP",
    #     main=paste("4-Fold Cross-validation,\n Year = ", time))
    #abline(0, 1, lty=2)
    #text(0.6, 0.2, paste("AUC of FLIPI+ALS=",round(xx$AUC,3)) )
    
    xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob2, 
                      predict.time=time, method="NNE", span=0.25*nobs^(-0.2))
    #points(xx$FP, xx$TP, type="l",col="gray")
    AUC.old <- xx$AUC
    #text(0.6, 0.1, paste("AUC of FLIPI=",round(xx$AUC,3)) )
  }
  boot.AUC[b, ] <- c(AUC.new, AUC.old, w)
}

## only 190 results with no warnings for 60 months
table( boot.AUC[,3])
quantile( boot.AUC[which( boot.AUC[,3]==0), 1], probs=c(0.025, 0.975)) ## 0.601 - 0.737
quantile( boot.AUC[which( boot.AUC[,3]==0), 2], probs=c(0.025, 0.975)) ## 0.583 - 0.713

### 60 months: 2018-08-31
#  FLIPIALS: 0.603 to 0.731
#  FLIPI   : 0.583 to 0.713 

## only 190 results with no warnings for 120 months
table( boot.AUC[,3])
quantile( boot.AUC[which( boot.AUC[,3]==0), 1], probs=c(0.025, 0.975)) ## 0.613 - 0.740
quantile( boot.AUC[which( boot.AUC[,3]==0), 2], probs=c(0.025, 0.975)) ## 0.601 - 0.725

### 120 months: 2018-08-31
#  FLIPIALS: 0.6182153 0.7416418 
#  FLIPI: 0.6080163 0.7265946 






##### Test for the hypothesis AUC of 0.5
require( bmrm)
require(survivalROC)
## k=4 to have at least two subjects in each FLIPIALS category
set.seed(200)
k = 4
tmp <- dd
perm.AUC.new <- NULL
perm.AUC.old <- NULL

for(b in 549:1000){
  id.perm <- sample(1:nrow(dd), nrow(dd), replace = F)
  tmp$OS <- dd$OS[id.perm]
  tmp$death <- dd$death[id.perm]
  
  folds <- balanced.cv.fold(y=tmp$FLIPIALS, num.cv=k)
    table( tmp$FLIPIALS, folds)
  
  par(mfrow=c(1,2))  
  for(time in c(60, 120)){
    tmp$prob <- NA ## survival probability based on the new model 
    tmp$prob2 <- NA ## survival probability based on the old model
    
    for(j in 1:k){
      # j <- 1
      index.test <- which(folds==j)
      train <- tmp[-index.test, ]
      test <- tmp[index.test, ]
      fit.train <- coxph(Surv(OS, death)~FLIPIALS, data=train)
      #predict(fit.train, newdata=test, type="lp")  
      bh <- basehaz(fit.train)
      H.120 <- bh[ which.min(abs(bh$time-time)), ]$hazard
      xb <- predict(fit.train, newdata=test, type="lp")
      surv.prob.120 <- (exp(-1*H.120))^exp(xb)
      tmp$prob[index.test] <- surv.prob.120
      
      fit.train <- coxph(Surv(OS, death)~FLIPI, data=train)
      #predict(fit.train, newdata=test, type="lp")  
      bh <- basehaz(fit.train)
      H.120 <- bh[ which.min(abs(bh$time-time)), ]$hazard
      xb <- predict(fit.train, newdata=test, type="lp")
      surv.prob.120 <- (exp(-1*H.120))^exp(xb)
      tmp$prob2[index.test] <- surv.prob.120
    }
    
    nobs <- nrow(tmp)
    xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob, 
                      predict.time=50, method="NNE", span=0.25*nobs^(-0.2))
    perm.AUC.new[b] <- xx$AUC
    
    xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob2, 
                      predict.time=50, method="NNE", span=0.25*nobs^(-0.2))
    perm.AUC.old[b] <- xx$AUC
  }
}


sum(perm.AUC.old>0.72) / 1000
#[1] 0
sum(perm.AUC.old>0.77) / 1000
#[1] 0
##sum(perm.AUC.new > perm.AUC.old)
##[1] 521
##sum(perm.AUC.new < perm.AUC.old)
##[1] 479
save.image("~/Documents/Project/Yang, George/RData/2018-07-31.RData")














####### Cross-validated AUC for Transformation #######
## step 1. Split a complete data into k-folds. 
## step 2. Perform cross-validation and save prediction scores of a test set in each test-set fold.
## step 3. Generate AUC of the saved prediction score against Overall survival of all data

# https://www.rdocumentation.org/packages/bmrm/versions/3.7/topics/balanced.cv.fold
require( bmrm)
require(survivalROC)
## k=4 to have at least two subjects in each FLIPIALS category
set.seed(200)
k = 4
folds <- balanced.cv.fold(y=dd$FLIPIALS, num.cv=k)
table( dd$FLIPIALS, folds)

dd$y <- dd$TransformedDisease1
dd$FLIPALS.num <- as.numeric(dd$FLIPIALS)
dd$FLIPI.num <- as.numeric(dd$FLIPI)


dd$FLIPALS.num <- (dd$FLIPIALS)
dd$FLIPI.num <- (dd$FLIPI)


barplot( t(prop.table( table(dd$FLIPIALS, dd$y), margin=1 )[,2:1]) )

  tmp <- dd
  tmp$prob <- NA ## survival probability based on the new model 
  tmp$prob2 <- NA ## survival probability based on the old model
  
  for(j in 1:k){
    ## FLIPI + ALS Model
    # j <- 1
    index.test <- which(folds==j)
    train <- tmp[-index.test, ]
    test <- tmp[index.test, ]
    #fit.train <- glm(y~FLIPIALS, data=train, family="binomial")
    fit.train <- glm(y~FLIPALS.num, data=train, family="binomial")
      xb.train <- predict(fit.train, newdata=train)
      #plot(fit.train)
      summary(fit.train)
        
    # Prediction on test set : predict(fit.train, newdata=test, type="lp")  
    xb <- predict(fit.train, newdata=test, type="response")
      plot(test$FLIPALS.num, xb)
      mean.y <- tapply(test$y, xb, mean)
      plot( names(mean.y), mean.y)
    tmp$prob[index.test] <- xb    
    
    #fit.train <- glm(y~FLIPIALS, data=train, family="binomial")
    fit.train <- glm(y~FLIPI.num, data=train, family="binomial")
    xb.train <- predict(fit.train, newdata=train)
    #plot(fit.train)
    summary(fit.train)
    
    # Prediction on test set : predict(fit.train, newdata=test, type="lp")  
    xb <- predict(fit.train, newdata=test, type="response")
      plot(test$FLIPI, xb)
    mean.y <- tapply(test$y, xb, mean)
    plot( names(mean.y), mean.y)
    tmp$prob2[index.test] <- xb    
  }

  par(mfrow=c(1,1))  
  require( "ROCR")
  ## New model  
  xx <- ROCR::prediction(predictions=tmp$prob, labels=tmp$y)
  perf.A <- ROCR::performance(prediction.obj = xx, "tpr", "fpr")  
  plot( perf.A, type="l", lwd=3, col="black", cex.lab=1.3, cex.axis=1.3)
  abline(0,1, lty=2)
  s <- ROCR::performance(prediction.obj = xx, measure = "auc" )
  s@y.values ## 0.608
  
  xx <- ROCR::prediction(predictions=tmp$prob2, labels=tmp$y)
  perf.B <- ROCR::performance(prediction.obj = xx, "tpr", "fpr" )  
  points( perf.B@x.values[[1]], perf.B@y.values[[1]], type="l", lwd=3, col="darkgray")
  s <- ROCR::performance(prediction.obj = xx, measure = "auc" )
  s@y.values ## 0.579

  legend(0.4, 0.3, legend=c("FLIPI+ALS", "FLIPI"), col=c("black", "darkgray"), lwd=5, bty="n")
  
  require(cvAUC)
  cvAUC::ci.cvAUC()
  cvAUC::ci.cvAUC(predictions=xx, labels=tmp$y, confidence=0.95)
    
 
  
 
######## Calculate cross-validated AUC ########
  iid_example <- function(data, V=10){
    .cvFolds <- function(Y, V){  #Create CV folds (stratify by outcome)
      Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
      Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
      folds <- vector("list", length=V)
      for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
      return(folds)
    }
    .doFit <- function(v, folds, data){  #Train/test glm for each fold
      fit <- glm(Y~., data=data[-folds[[v]],], family=binomial)
      pred <- predict(fit, newdata=data[folds[[v]],], type="response")
      return(pred)
    }
    folds <- .cvFolds(Y=data$Y, V=V)  #Create folds
    predictions <- unlist(sapply(seq(V), .doFit, folds=folds, data=data))  #CV train/predict
    predictions[unlist(folds)] <- predictions  #Re-order pred values
    # Get CV AUC and confidence interval
    out <- ci.cvAUC(predictions=predictions, labels=data$Y, folds=folds, confidence=0.95)
    return(out)
  }
  
  set.seed(1)
  zz <- data.frame(Y=tmp$y, x=tmp$FLIPALS.num)
  out <- iid_example(data=zz)
  out  
  
  zz <- data.frame(Y=tmp$y, x=tmp$FLIPI.num)
  out.2 <- iid_example(data=zz)
  out.2

    
  
  
    require(AUC)
  auc( roc(predictions = tmp$prob, labels = tmp$y) )
   
    
  
    
  AUC.old <- xx$AUC
  plot(xx$FP, xx$TP, type="l", xlim=c(0, 1), ylim=c(0,1), 
       xlab=paste( "False Positive"), 
       ylab="True Positive", lwd=3, 
       main=paste("4-Fold Cross-validation,\n Year = ", time))
  abline(0, 1, lty=2)
  text(0.6, 0.2, paste("AUC of FLIPI+ALS=",round(xx$AUC,3)) )
  
  xx <- survivalROC(Stime=tmp$OS, status=tmp$death, marker=-1*tmp$prob2, 
                    predict.time=50, method="NNE", span=0.25*nobs^(-0.2))
  points(xx$FP, xx$TP, type="l",col="gray", lwd=3)
  AUC.new <- xx$AUC
  text(0.6, 0.1, paste("AUC of FLIPI=",round(xx$AUC,3)) )
}













colMeans(d.result)

d.observed





bh <- basehaz(fit)
plot( bh[,2], bh[,1], main="Cumulative Hazard Function")
bh


log.Cum.Hazard <- predict(fit, type="expected")
surv.prob <- exp(-1*log.Cum.Hazard)
cutpoints <- quantile(surv.prob, probs=seq(0, 1, 0.1))


f <- cph(Surv(OS, death)~FLIPIALS, data = dd, x=TRUE, y=TRUE, 
         time.inc=12, surv=TRUE)
cal <- calibrate(f, u=120, method="crossvalidation", B=200)
plot(cal)

calkm <- calibrate(f, u=120, cmethod="KM", B=200)
plot( calkm, add=TRUE)

## KM curve
plot(fit, col=colfunc(6), lwd=4)
  summary(coxph( model.form, data=sdata ))$coef
  s <- summary(fit)
  x <- tapply(sdata$OS, INDEX=sdata$FLIPI, FUN = max)
  y <- tapply(s$surv, s$strata, FUN=min)
  text(x+5, y, names(x) )
#ggsurvplot(fit, data = sdata, risk.table = TRUE, pval=TRUE)


#### FLIPIALS model  
model.form <- formula("Surv(OS, death) ~  FLIPIALS")
fit <- survfit( Surv(OS, death) ~ FLIPIALS,data=sdata )
survdiff( model.form,data=sdata )
## KM curve
plot(fit, col=colfunc(7), lwd=5)
summary(coxph( model.form, data=sdata ))$coef
s <- summary(fit)
x <- tapply(sdata$OS, INDEX=sdata$FLIPIALS, FUN = max)
y <- tapply(s$surv, s$strata, FUN=min)
text(x+5, y, names(x) )
#ggsurvplot(fit, data = sdata, risk.table = TRUE, pval=TRUE)



##### https://www.rdocumentation.org/packages/caret/versions/6.0-80/topics/createDataPartition
##### permuted
dd.perm <- dd[sample(nrow(dd)),]
nfold<-5
k<-5

folds <- cut(seq(1,nrow(dd.perm)),breaks=5,labels=FALSE)

surv_prob <- rep(list(rep(NA)),k)
time <- rep(list(rep(NA)),k)
cens <- rep(list(rep(NA)),k)

#pdf( file="KM.pdf" )

calibration.slope <- NULL
for(i in 1:5){
	# i <- 3
	#Segment data by fold using the which() function
	testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- dd.perm[testIndexes, ]
	trainData <- dd.perm[-testIndexes, ]

	#### training model
	#model <- coxph(Surv(OS, death)~  factor(FLIPI) + AbsLymphScore1lt10gt1,data=trainData)# fit the model with train data
	# n=588
	par( mfrow=c(2,3))
	fit <- survfit(Surv(OS, death)~  FLIPIALS,data=trainData)# fit the model with train data
  plot(fit, col=colfunc(7), lwd=5, xlim=c(0, 250))
  title("Train (n=588)")
  	summary(coxph( model.form, data=sdata ))$coef
  	s <- summary(fit)
  	x <- tapply(sdata$OS, INDEX=sdata$FLIPIALS, FUN = max)
  	y <- tapply(s$surv, s$strata, FUN=min)
  	text(x+5, y, names(x) )
  
	model.train <- coxph(Surv(OS, death)~  FLIPIALS,data=trainData)# fit the model with train data
	s <- summary(model.train)
	s
	coef.train <- s$coefficients  
	## linear predictor value
	xb.train <- predict(model.train, newdata=trainData, type="lp")
	
	
	
	fit.test <- survfit(Surv(OS, death)~  FLIPIALS,data=testData)# fit the model with train data
	plot(fit.test, col=colfunc(7), lwd=5, xlim=c(0, 250))
	title("Test (n=147)")
  	#summary(coxph( model.form, data=sdata ))$coef
  	s <- summary(fit.test)
  	x <- tapply(sdata$OS, INDEX=sdata$FLIPIALS, FUN = max)
  	y <- tapply(s$surv, s$strata, FUN=min)
  	text(x+5, y, names(x) )	

  	
  ### The higher the FLIPIALS is, the higher the prognostic index is.
  plot( trainData$FLIPIALS, xb.train)
  #par( mfrow=c(2,1))
  hist(xb.train, xlim=c(-2.5, 2.5), main="trainingset" )
  coxph( Surv(OS, death) ~ xb.train, data=trainData) ## coefficients should be 1.
  
  
	## generate a prediction for test set
	results_prob <- predict(model,newdata=testData,type="lp")	
	surv_prob[[i]]<-(results_prob)
	time[[i]]<- testData$RECALCULATEDSurvivalTimeMonths
	cens[[i]]<-testData$VS

	
  #### test set analysis	
	testData <- data.frame( testData)
	xb <- predict(model.train,newdata=testData,type="lp", xlim=c(-2.5, 2.5),main="validationset")
  	hist( xb )
	## cox regression with the predicted linear predictor
	## If the model is well calibrated, this model should have a coefficient of 1 for xb.
  fit.valid <- coxph( Surv(OS, death) ~ xb, data=testData)
  s <- summary(fit.valid)
  
  ### Test for Calibration Slope
  beta <- s$coefficients[1]
  se.beta <- s$coefficients[3]
  beta0 <- 1
  p <- 2*(1 - pnorm( abs(beta-beta0)/se.beta) )
  p ## 0.6139
  
  calibration.slope <- rbind(calibration.slope,  c(i, beta, p))
  
  ### offset
  #fit.offset <- coxph( Surv(OS, death) ~ offset(xb), data=testData)
  #  fit.offset
  fit.offset <- coxph( Surv(OS, death) ~ FLIPIALS + offset(xb), data=testData)
    fit.offset
    summary(fit.offset) ## see its score (log rank) test
  coxph( Surv(OS, death) ~ FLIPIALS, data=testData)
  coxph.fit( y=Surv(testData$OS, testData$death), x=matrix(FLIPIALS), offset=offset)

}




	par( mfrow=c(2,2))
	for(num.group in 3:6){
		# num.group <- 3
	  #ff <- survfit( Surv(OS, death) ~ xb, data=testData)
	  #plot( ff, lwd=3, col=1:12, mark =letters[1:12] )
	  
	  # xb <- rnorm( nrow(testData), 0, 1)
		# testData$group <- cut(xb, breaks=num.group)
	  ##testData$group <- cut(xb, breaks=c(min(xb), -0.47, 0.2, max(xb)), include.lowest=T )
	  testData$group <- cut(xb, breaks=quantile(xb, probs=seq(0,num.group, 1)/num.group), include.lowest=T )
	  
	  table(testData$group)
		
		## Kaplan-Meier 
		ff <- survfit( Surv(OS, death) ~ group, data=testData)
			ff
		plot( ff, mark=letters[1:3],
			lwd=3,
			col=c("red", "black", "green", "pink", "orange", "darkgray")[1:num.group] )
	}
	
}
dev.off()



ggsurvplot(fit, data = sdata, risk.table = TRUE)
plot(fit, col=colfunc(7), lwd=4)




fit <- survfit( Surv(OS, death)~  factor(FLIPI) + AbsLymphScore1lt10gt1,data=sdata )
ggsurvplot(fit, data = sdata, risk.table = TRUE)


