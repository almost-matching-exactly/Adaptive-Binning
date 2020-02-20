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
library(latex2exp)
require(foreign)
require(xtable)
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

m_estimators <- c('Full Matching', 'Prognostic', 'CEM',
                  'Mahalanobis', 'Nearest Neighbor')
all_means <- data.frame(matrix(NA, 6, 4))
names(all_means) = c('Method', 'NSW', 'CPS', 'PSID')
all_means$Method = c(m_estimators, 'AHB')
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
  fhat0 = colMeans(predict(bart_fit, test = as.matrix(test_covs)))
  fhat1 = fhat0
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
all_means[6, 2:4] = c(mnsw, mcps, mpsid)

CATEs <- cbind(test_covs[test_treated, ], Y=nswt$Y,
               data.frame(nsw=all_mip_cates[[1]], cps=all_mip_cates[[1]], psid=all_mip_cates[[2]]))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
modal_bins <- sapply(4:16, function(v) apply(round(matrix(all_mip_bins[[2]][CATEs$education==v, 2, ], ncol=2)), 2, getmode))
ed_table <- CATEs %>%
  group_by(education) %>%
  summarize(nsw = mean(nsw[all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95)], na.rm = TRUE),
            cps = mean(cps[all_mip_objs[[2]] < quantile(all_mip_objs[[2]], .95)], na.rm = TRUE),
            psid = mean(psid[all_mip_objs[[3]] < quantile(all_mip_objs[[3]], .95)], na.rm = TRUE))
ed_table = cbind(ed_table, modal_lb=modal_bins[1, ], modal_ub=modal_bins[2, ])
cat(print(xtable(ed_table, digits = 0)), file='ed_table.tex')



cov <- 7
bins <- all_mip_bins[[1]][all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95),,]
bindata <- data.frame(unique(bins[,cov,]))
names(bindata) <- c('a', 'b')
simdata = data.frame(x=test_df[, cov], tr=test_df$treated, y=test_df$Y)
ggplot(simdata) + 
  geom_rect(data=bindata, aes(xmin=a, ymin=min(simdata$y[simdata$tr==1], na.rm=T), xmax=b, 
                              ymax=max(simdata$y[simdata$tr==1],na.rm=T)), 
            color="black", size=0.5, alpha=0.1, fill="grey") +
  geom_point(aes(x=x, y=y, color=tr, size=tr)) + 
  ylim(c(min(simdata$y[simdata$tr==1]), max(simdata$y[simdata$tr==1]))) +
  xlim(c(min(simdata$x[simdata$tr==1]), max(simdata$x[simdata$tr==1]))) +
  scale_size_manual(values=c(1, 2), guide='none') + 
  scale_color_manual(values=c(brewer.pal(3, 'Set1')[1], 'black'), labels=c('PSID control', 'NSW treated')) +
  xlab("Income 1974 in USD") + ylab("Income 1978 in USD (Outcome)") + theme_bw() + 
  theme(legend.position = c(0.87,0.9), legend.background = element_rect(color="black", size=0.5),
        legend.title = element_blank(), text = element_text(size=20))
ggsave('lalonde_income_bins.png', width = 10, height = 6, units = "in")



cov <- 7
bins <- all_mip_bins[[1]][all_mip_objs[[1]] < quantile(all_mip_objs[[1]], .95),,]
bindata <- data.frame(unique(bins[,cov,]))
names(bindata) <- c('a', 'b')
simdata = data.frame(x=CATEs[, cov],  y=CATEs$Y)
ggplot(simdata) + 
  #geom_rect(data=bindata, aes(xmin=a, ymin=min(simdata$y, na.rm=T), xmax=b,ymax=max(simdata$y,na.rm=T)),
  #         color="black", size=0.5, alpha=0, fill="grey", linetype=2) +
  geom_linerange(data=bindata, aes(x=b, ymin=0, ymax=40000), linetype=2) +
  geom_point(aes(x=x, y=y)) + 
  geom_smooth(aes(x=x,y=fhat0[test_treated]), se = FALSE, linetype=1, size=1.5, color="black") +
  annotate("text", x = 31000, y=21000, label=TeX("\\hat{f}_0(\\mathbf{x}) (Smoothed)"), size=6) + 
  #scale_color_manual(values=c(brewer.pal(3, 'Set1')[1], 'black'), labels=c('PSID control', 'NSW treated')) +
  ylim(c(0, 40000)) + 
  xlab("Income 1974 in USD") + ylab("Estimated ITE (PSID) in US$") + theme_minimal() + 
  theme(legend.position = c(0.87,0.9), legend.background = element_rect(color="black", size=0.5),
        legend.title = element_blank(), text = element_text(size=20))
ggsave('lalonde_income_bins.png', width = 10, height = 6, units = "in")


## other methods
naive_means = rep(NA, length(all_dfs))
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
  
  ## Other estimators
  nn_out = matchit(f, df, method='nearest')
  nn_est = (df$Y[df$treated==1] %*% nn_out$weights[df$treated==1])/sum(nn_out$weights[df$treated==1]) - 
    (df$Y[df$treated==0] %*% nn_out$weights[df$treated==0])/sum(nn_out$weights[df$treated==0])
  all_means[1, ds+1] = nn_est
  fm_out = matchit(f, df, method='full', ratio=1)
  fm_est = (df$Y[df$treated==1] %*% fm_out$weights[df$treated==1])/sum(fm_out$weights[df$treated==1]) - 
    (df$Y[df$treated==0] %*% fm_out$weights[df$treated==0])/sum(fm_out$weights[df$treated==0])
  all_means[2, ds+1] = fm_est
  cem_out = matchit(f, df, method='cem', ratio=1)
  cem_est = (df$Y[df$treated==1] %*% cem_out$weights[df$treated==1])/sum(cem_out$weights[df$treated==1]) - 
    (df$Y[df$treated==0] %*% cem_out$weights[df$treated==0])/sum(cem_out$weights[df$treated==0])
  all_means[3, ds+1] = cem_est
  mah_out = matchit(f, df,distance = 'mahalanobis', method = "nearest", ratio=1)
  mah_est = (df$Y[df$treated==1] %*% mah_out$weights[df$treated==1])/sum(mah_out$weights[df$treated==1]) - 
    (df$Y[df$treated==0] %*% mah_out$weights[df$treated==0])/sum(mah_out$weights[df$treated==0])
  all_means[4, ds+1] = mah_est
  prg_est = mean(est_prog(df, counterfactuals))
  all_means[5, ds+1] = prg_est
  
  naive_means[ds] = mean(df$Y[df$treated==1]) - mean(df$Y[df$treated==0])
}

all_means = rbind(all_means, c('Naive', naive_means))
all_means = data.frame('Method'=all_means[, 1], apply(all_means[, 2:4], 2, as.numeric))
all_err = all_means[, 2:4] - all_means[7,2]
err_tab = data.frame(Method=all_means$Method, CPS=all_means$CPS, CPS_err=all_err$CPS, 
                     PSID=all_means$PSID, PSID_err=all_err$PSID)
cat(print(xtable(err_tab, digits=0)), file='lalonde_tab.tex')


#units = c(66,10, 68)
units = which(CATEs$black==1 & CATEs$married==0 & CATEs$age < 25 & 
                (CATEs$education == 11|CATEs$education==12) & CATEs$nodegree==0)
units = which((CATEs$education == 11|CATEs$education==12) )
cov1 = 2
cov2 = 1
bins = all_mip_bins[[2]][units, ,]
bindata = data.frame(a1=bins[, cov1, 1], a2=bins[, cov2, 1], b1=bins[, cov1, 2], b2=bins[, cov2, 2])
#bindata$nodegree = CATEs$psid[units] > 0
bindata$nodegree = CATEs$nodegree[units]
simdata = data.frame(x1=test_df[, cov1], x2=test_df[, cov2], tr=test_df$treated, y=test_df$Y)
simdata$selected = 1:nrow(test_df) %in% test_treated[units]
ggplot(simdata) + 
  geom_rect(data=bindata, aes(xmin=a1, ymin=a2, xmax=b1, ymax=b2, color=as.logical(nodegree)), 
            size=0.5, alpha=0) + 
  geom_point(aes(x=x1, y=x2,size=as.logical(selected))) +
  ylim(c(min(bindata$a2), max(bindata$b2))) + 
  xlim(c(min(bindata$a1), max(bindata$b1))) + 
  scale_size_manual(values=c(1, 4)) + 
  scale_color_manual(values=c('red', 'black')) + 
  xlab("x1") + ylab("x2") +  labs(color="Degree") +  theme_bw() + 
  theme(legend.position = c(0.95,0.85), 
        legend.background = element_rect(color="black", size=0.5))


cov=2
bindata <- data.frame(bins[,cov,])
names(bindata) <- c('a', 'b')
simdata = data.frame(x=CATEs[units, cov], y=CATEs$cps[units], a=bindata$a, b=bindata$b, 
                     nodegree=CATEs$nodegree[units])
ggplot(simdata) + 
  geom_rect(aes(xmin=a, ymin=min(simdata$y, na.rm=T), xmax=b, 
                ymax=max(simdata$y,na.rm=T), color=as.logical(nodegree)), 
            size=0.5, alpha=0, fill="grey") +
  geom_point(aes(x=x, y=y), size=2) + 
  scale_color_manual(values=c('red', 'black'), labels=c('control', 'treated')) +
  xlab("x") + ylab("y") + theme_bw() + 
  theme(legend.position = c(0.95,0.9), legend.background = element_rect(color="black", size=0.5),
        legend.title = element_blank())


mg = make_mg(test_covs, all_mip_bins[[2]][66,,1], all_mip_bins[[2]][66,,2])
mg1 = rbind(test_df[66, ], test_df[mg, ][c(which(!test_df$treated[mg])), ])
  
mg = make_mg(test_covs, all_mip_bins[[2]][181,,1], all_mip_bins[[2]][181,,2])
mg2 = test_df[mg, ]

cat(print(xtable(rbind(mg1, mg2), digits=0)), file='mg_tables.tex')

