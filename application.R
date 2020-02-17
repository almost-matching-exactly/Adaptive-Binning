setwd(path.expand("~/Dropbox/Duke/projects/FLAME-binning/Adaptive-Binning/"))
source("AB_MIPs.R")
library(ggplot2)
library(reshape2)
library(dbarts)
library(RColorBrewer)
require(dbarts)
require(MatchIt)
require(cem)
require(tidyverse)
require(foreign)
require(glmnet)
source("sims.R")
set.seed(42069)

standardize <- function(x){
  list(x=(x - mean(x, na.rm=T))/sd(x, na.rm=T), m=mean(x, na.rm=T), sddev=sd(x, na.rm=T))
}

n_train_t <- 100
train_pct <- .8
n_train_c <- 100

# Preprocess data and adapt to our format
nsw <- read.dta("Lalonde/nsw_dw.dta")
nswt <- nsw[nsw$treat==1, ]
nswc <- nsw[nsw$treat==0, ]
cps1 <- read.dta('Lalonde/cps_controls.dta')
cps2 <- read.dta('Lalonde/cps_controls2.dta')
cps3 <- read.dta('Lalonde/cps_controls3.dta')
psid1 <- read.dta('Lalonde/psid_controls.dta')
psid2 <- read.dta('Lalonde/psid_controls2.dta')
psid3 <- read.dta('Lalonde/psid_controls3.dta')

#all_dfs = list(nswc, cps1, cps2, cps3, psid1, psid2, psid3)
all_dfs = list(nswc, rbind(cps1, cps2, cps3), rbind(psid1, psid2, psid3))

# Treated units
nswt <- data.frame(nswt[, 3:(ncol(nswt)-1)], treated=nswt$treat, Y=nswt$re78)
train_t_ids <- sample(nrow(nswt), n_train_t, replace = F)
nswt_tr <- nswt[train_t_ids, ]
nswt_ts <- nswt[-train_t_ids, ] 

all_bart_cates <- list()
all_mip_cates <- list()
all_mip_cates_both <- list()
all_mip_objs <- list()
all_mip_bins <- list()
message(paste('Starting analysis, experimental ATE =', mean(nswt$Y) - mean(nswc$re78)))
for (ds in 1:length(all_dfs)){
  message("###########################################")
  message(paste('Analyzing dataset', ds, 'of', length(all_dfs)))
  message("###########################################")

  cdf <- all_dfs[[ds]]
  cdf <- data.frame(age=cdf$age, education=cdf$education, black=cdf$black, hispanic=cdf$hispanic,
                   married=cdf$married, nodegree=cdf$nodegree, re74=cdf$re74, re75=cdf$re75,
                   treated=cdf$treat, Y=cdf$re78)
  
  n_train_c <- floor(train_pct * nrow(cdf))
  n_train = n_train_c
  cdf_tr_ids <- sample(nrow(cdf), n_train_c, replace = F)
  cdf_tr <- cdf[cdf_tr_ids, ]
  cdf_ts <- cdf[-cdf_tr_ids, ]
  
  #df <- rbind(nswt_tr, cdf_tr, nswt_ts, cdf_ts)
  df <- rbind(cdf_tr, nswt, cdf_ts)
  df$treated <- as.logical(df$treated)
  
  n = nrow(df)
  inputs <- estimator_inputs(df, n_train, n, 'BART', cv = FALSE)
  c(df, f, n, n_train, p,
    train_df, train_covs, train_control, train_treated,
    test_df, test_covs, test_control, test_treated, 
    n_test_control, n_test_treated, 
    bart_fit, counterfactuals) %<-% inputs
  
  # BART
  #fhat1 = colMeans(predict(bart_fit, test = as.matrix(cbind(test_covs, treated=1))))
  #fhat0 = colMeans(predict(bart_fit, test = as.matrix(cbind(test_covs, treated=0))))
  fhat0 = colMeans(predict(bart_fit, test = as.matrix(test_covs)))
  
  ## GLMNET v2
  # m = glmnet(x=as.matrix(cbind(train_covs, train_df$treated)),
  #            y=train_df$Y, family="gaussian", alpha=0.3)
  # l = cv.glmnet(x=as.matrix(cbind(train_covs, train_df$treated)),
  #               y=train_df$Y, family="gaussian", alpha=0.3)$lambda.min
  # fhat0 = predict(m, newx=as.matrix(cbind(test_covs, treated=0)), type='response', s=l)
  # fhat1 = predict(m, newx=as.matrix(cbind(test_covs, treated=1)), type='response', s=l)

  # m = glmnet(x=as.matrix(train_covs),
  #            y=train_df$Y, family="gaussian", alpha=0.5)
  # l = cv.glmnet(x=as.matrix(train_covs),
  #               y=train_df$Y, family="gaussian", alpha=0.5)$lambda.min
  # fhat0 = predict(m, newx=as.matrix(test_covs), type='response', s=l)
  fhat1 = fhat0
  # ## OLS
  # dw_fmla <- as.formula(Y ~ treated + age + age**2 + education +
  #                         black + hispanic + married + nodegree + re74 + re75)
  # ols <- (lm(dw_fmla, data=train_df))
  # fhat1 <- predict(ols, newdata=(cbind(test_covs, treated=TRUE)))
  # fhat0 <- predict(ols, newdata=(cbind(test_covs, treated=FALSE)))

  # bcates = (fhat1 - fhat0)
  # all_bart_cates[[ds]] <- bcates[test_treated]
  # message(paste("Black-box ATE =", mean(bcates[test_treated])))
  
  bcates = (test_df$Y - fhat0)
  all_bart_cates[[ds]] <- bcates[test_treated]
  message(paste("Black-box ATE =", mean(bcates[test_treated])))

  alpha=0
  lambda0=0
  lambda1=0
  gamma0=7.5
  gamma1=0
  beta=2
  m=1
  M=1e5
  
  ## MIP
  mip_cates = vector('numeric', n_test_treated)
  mip_cates_both = vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  mip_objs = vector('numeric', n_test_treated)
  message("Running MIP")
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    if(n_test_treated + n_test_control > 300){
      ord = order(abs(fhat1[i] - fhat1) * gamma1/sd(fhat1) + abs(fhat0[i] - fhat0) * gamma0/sd(fhat0))
      #cands = c(which(test_df$treated==1), ord[test_df$treated[ord]==0][1:50])
      cands = ord[test_df$treated[ord]==0][1:100]
    }else{
      cands = (1:(n - n_train))[-i]
    }
      
    mip_pars =  setup_miqp_fhat(xi = as.numeric(test_covs[i, ]),
                                x_test = as.matrix(test_covs[cands, ]),
                                z_test = test_df$treated[cands],
                                fhati1=fhat1[i],
                                fhati0=fhat0[i],
                                fhat1=fhat1[cands],
                                fhat0=fhat0[cands], 
                                alpha=alpha, lambda0=lambda0, lambda1=lambda1,
                                gamma0=gamma0/sd(fhat0[cands]), gamma1=gamma1/sd(fhat1[cands]), 
                                beta=beta, m=m, M=M)

    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, length(cands), p)
    
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b
    
    mg = make_mg(test_covs, mip_out$a, mip_out$b)
    mg = mg[mg %in% cands]
    #mg = cands[mip_out$w > 1/M]
    mip_cates[l] = test_df$Y[i] - mean(test_df$Y[mg][!test_df$treated[mg]])
    mip_cates_both[l] = mean(test_df$Y[mg][test_df$treated[mg]]) - mean(test_df$Y[mg][!test_df$treated[mg]])
    mip_objs[l] = mip_out$obj
  }
  message(paste('MIP ATE =', mean(mip_cates, na.rm=T), 'MIP ATE Both =', mean(mip_cates_both, na.rm=T)))
  all_mip_bins[[ds]] <- mip_bins
  all_mip_cates[[ds]] <- mip_cates
  all_mip_cates_both[[ds]] <- mip_cates_both
  all_mip_objs[[ds]] <- mip_objs
}

mnsw = mean(all_mip_cates[[1]][all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95)], na.rm=T)
mcps = mean(all_mip_cates[[2]][all_mip_objs[[2]] < quantile(all_mip_objs[[2]], .95)], na.rm=T)
mpsid = mean(all_mip_cates[[3]][all_mip_objs[[3]] < quantile(all_mip_objs[[3]], .95)], na.rm=T)

CATEs <- cbind(test_covs[test_treated, ], Y=nswt$Y,
               data.frame(nsw=all_mip_cates[[1]], cps=all_mip_cates[[1]], psid=all_mip_cates[[2]]))
CATEs %>%
  group_by(education) %>%
  summarize(nsw = mean(nsw[all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95)], na.rm = TRUE),
            cps = mean(cps[all_mip_objs[[2]] < quantile(all_mip_objs[[2]], .95)], na.rm = TRUE),
            psid = mean(psid[all_mip_objs[[3]] < quantile(all_mip_objs[[3]], .95)], na.rm = TRUE))

cov <- 7
bins <- all_mip_bins[[1]][all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95),,]
bindata <- data.frame(unique(bins[,cov,]))
names(bindata) <- c('a', 'b')
simdata = data.frame(x=test_df[, cov], tr=test_df$treated, y=test_df$Y)
ggplot(simdata) + 
  geom_rect(data=bindata, aes(xmin=a, ymin=min(simdata$y[simdata$tr==1], na.rm=T), xmax=b, 
                              ymax=max(simdata$y[simdata$tr==1],na.rm=T)), 
            color="black", size=0.5, alpha=0.1, fill="grey") +
  geom_point(aes(x=x, y=y, color=tr), size=2) + 
  ylim(c(min(simdata$y[simdata$tr==1]), max(simdata$y[simdata$tr==1]))) +
  xlim(c(min(simdata$x[simdata$tr==1]), max(simdata$x[simdata$tr==1]))) +
  scale_color_manual(values=c(brewer.pal(3, 'Set1')[1], brewer.pal(3, 'Set1')[2]), labels=c('control', 'treated')) +
  xlab("Income 1974") + ylab("Income 1978 (Outcome)") + theme_bw() + 
  theme(legend.position = c(0.95,0.9), legend.background = element_rect(color="black", size=0.5),
        legend.title = element_blank())

