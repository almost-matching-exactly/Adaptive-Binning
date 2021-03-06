library(RColorBrewer)
library(ggiraph)

size_vs_error_plot <- function(sim_out, estimator, sim_number = 1) {
  CATEs <- sim_out$CATEs %>%
    filter(estimator == !!enquo(estimator))
  bins <- sim_out$bins[[sim_number]] %>%
    extract2(estimator)
  
  areas <- 
    (bins[, , 2] - bins[, , 1]) %>%
    apply(1, prod)
  
  ggdata <- 
    data.frame(Error = abs(CATEs$actual - CATEs$predicted), 
               Area = areas)
  
  p <- 
    ggplot(data = ggdata, aes(x = Area, y = Error)) + 
    geom_point() + 
    labs(x = 'Box Area', y = 'CATE Error')
  print(p)
}


bin_plot <- function(bins, test_df, cov1, cov2, fhat, background = FALSE, ...) {
  # browser()
  # Input covariates / dimensions for which you would like bins
  boundaries <- c(...)
  test_control <- 
    test_df %>%
    filter(!treated) %>%
    dplyr::select(cov1, cov2) %>%
    `colnames<-`(c('X1', 'X2'))
  test_treated <- 
    test_df %>%
    filter(treated)
  bins <- 
    cbind(bins[, c(cov1, cov2), 1], 
          bins[, c(cov1, cov2), 2],
          test_treated[, cov1],
          test_treated[, cov2]) %>%
    as.data.frame() %>%
    `colnames<-`(c('lower1', 'lower2', 
                   'upper1', 'upper2',
                   'point1', 'point2')) %>%
    mutate(id = 1:nrow(.),
           unmatched = is.na(lower1)) 
  
  min_dim1 <- min(bins$point1, na.rm = TRUE)
  min_dim2 <- min(bins$point2, na.rm = TRUE)
  max_dim1 <- max(bins$point1, na.rm = TRUE)
  max_dim2 <- max(bins$point2, na.rm = TRUE)
  
  if (background) {
    predict_points <- 
      expand.grid(seq(min_dim1, max_dim1, length.out = 50),
                  seq(min_dim2, max_dim2, length.out = 50)) %>%
      as.matrix() %>%
      `colnames<-`(c('X1', 'X2'))
    # browser()
      ###### FIX THIS FOR WHEN DIM > 2 and you need to integrate over other dimensions
    preds <- colMeans(predict(fhat, predict_points))
    
    predict_points %<>%
      as.data.frame() %>%
      `colnames<-`(c('X1', 'X2')) %>%
      mutate(preds = preds)
  }
  if (sum(is.na(bins$lower1)) == 0) {
    unmatched_units <- NULL
  }
  else {
    unmatched_units <- filter(bins, is.na(lower1)) 
  }
  
  
  p <- ggplot(data = bins)
  if (background) {
    p <- p + geom_raster(data = predict_points, aes(x = X1, y = X2, fill = preds), 
                         interpolate = TRUE)
  }
  p <- p + 
    geom_rect_interactive(aes(xmin = lower1, xmax = upper1,
                              ymin = lower2, ymax = upper2,
                              data_id = id),
                          fill = 'grey', color = 'NA',
                          size = 0.5, alpha = 0) +
    geom_point_interactive(data = filter(bins, !is.na(lower1)),
                           aes(x = point1, y = point2, data_id = id, color = 'Matched')) +
    geom_point(data = test_control, aes(x = X1, y = X2, color = 'Control'), size = 0.2) +
    scale_fill_gradient(guide="colorbar") + 
    theme_minimal() + 
    labs(x = 'x1', y = 'x2')
  
  if (is.null(unmatched_units)) {
    p <- p + 
      scale_color_manual(values = c('red', 'black'),
                         name = 'Type',
                         labels = c('Control', 'Matched'))
  }
  else {
    p <- 
      p + 
      geom_point(data = unmatched_units, aes(x = point1, y = point2, color = 'Unmatched')) +
      scale_color_manual(values = c('red', 'black', 'blue'),
                         name = 'Type',
                         labels = c('Control', 'Matched', 'Unmatched'))
  }
  
  if (length(boundaries) != 0) {
    boundaries %<>% 
      matrix(ncol = 4, byrow = TRUE) %>%
      as.data.frame() %>%
      `colnames<-`(c('Start1', 'End1', 'Start2', 'End2'))
    # browser()
    p <- 
      p + 
      geom_rect(data = boundaries, 
                aes(xmin = Start1, xmax = End1,
                ymin = Start2, ymax = End2),
                fill = NA, color = '#E9CCD1')
  }
  
  girafe(ggobj = p,
         options = list(
           opts_hover(
             css = girafe_css(
               css = 'fill:orange',
               area = 'fill:#E8C3D6;fill-opacity:0.8',
               point = 'fill:black;r:5;stroke:black'
             ))))
}

CATE_error_plot <- function(res) {
  ATT <- mean(res$actual)
  
  res$estimator %<>% 
    factor(levels = c('BART', 'Best CF', 
                      'MIP', 'Approximate', 
                      'CEM', 'Full Matching', 'Nearest Neighbor', 'Prognostic', 'Mahal', 'GenMatch'))
  
  res %<>% 
    mutate(error = abs(predicted - actual))
    
  perc_missing <- 
    res %>%
    group_by(estimator) %>%
    summarize(percent_missing = round(100 * mean(is.na(predicted)), 1))
  
  estimators <- perc_missing$estimator
  n_estimators <- length(estimators)
  
  group_means <- 
    res %>%
    group_by(estimator) %>%
    summarize(mean = mean(error, na.rm = TRUE))
  
  baseline_estimators <- intersect(estimators, 
                                   c('Approximate', 'MIP'))
  
  baseline_mean <- 
    group_means %>%
    filter(estimator %in% baseline_estimators) %>%
    pull(mean) %>%
    min()
  
  max_error_plotted <- 
    quantile(res$error, probs = .975, na.rm = T) / ATT
  
  res %<>%
    filter(is.na(predicted) | error / ATT <= max_error_plotted)
    
  lower_than <- group_means$mean <= baseline_mean
  lower_than[is.na(lower_than)] = FALSE
  
  colors = brewer.pal(3, "Set1")[if_else(lower_than, 2, 1)]
  colors[which(estimators %in% c('Greedy', 'MIP', 'MIP_both', 'MIP-Explain', 
                                 'MIP-Predict', 'MIQP-Variance'))] = brewer.pal(3, "Set1")[3]
  
  p <- 
    ggplot(data = res) + 
    geom_violin(aes(y = error / ATT, x = estimator, fill = estimator),
                color = "black", draw_quantiles = c(0.5), size = 0.2) + 
    geom_segment(x = 4.5, xend = 4.5, y = 0, yend = max_error_plotted) + 
    geom_segment(x = 2.5, xend = 2.5, y = 0, yend = max_error_plotted) + 
    geom_text(data = perc_missing, 
              aes(x = estimator, y = max_error_plotted + 0.25, label = percent_missing)) +
    annotate('text', x = 5.5, y = max_error_plotted + 0.5, label = 'Percent Missing') +
    scale_fill_manual(values = colors) + 
    ylim(c(0, max_error_plotted + 0.7)) +
    xlab("") + 
    ylab("Absolute estimation error\nRelative to ATT") + 
    theme_minimal() + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}

CATE_scatter_plot <- function(res){
  
  p <- ggplot(data=res) + 
       geom_point(aes(x=actual, y=predicted)) + 
       geom_abline(intercept = 0, slope=1, color = 'red', size = 1.5) + 
       facet_wrap(estimator~.) + 
       theme_minimal()
  
  print(p)
}

static_2d_bin_plot <- function(bins, test_df, cov1, cov2){
  bindata = data.frame(a1=bins[, cov1, 1], a2=bins[, cov2, 1], b1=bins[, cov1, 2], b2=bins[, cov2, 2])
  simdata = data.frame(x1=test_df[, cov1], x2=test_df[, cov2], y=test_df$treated)
  ggplot(simdata) + 
    geom_rect(data=bindata, aes(xmin=a1, ymin=a2, xmax=b1, ymax=b2), fill="grey", 
              color="black", size=0.5, alpha=0) + 
    geom_point(aes(x=x1, y=x2, color=y), size=2) +
    ylim(c(min(bindata$a2), max(bindata$b2))) + 
    xlim(c(min(bindata$a1), max(bindata$b1))) + 
    scale_color_manual(values=c('black', 'red')) + 
    xlab("x1") + ylab("x2") +  labs(color="Treated") +  theme_bw() + 
    theme(legend.position = c(0.95,0.85), 
          legend.background = element_rect(color="black", size=0.5))
}

static_1d_bin_plot <- function(bins, test_df, cov){
  bindata <- data.frame(unique(bins[,cov,]))
  names(bindata) <- c('a', 'b')
  simdata = data.frame(x=test_df[, cov], tr=test_df$treated, y=test_df$Y)
  ggplot(simdata) + 
    geom_rect(data=bindata, aes(xmin=a, ymin=min(simdata$y[simdata$tr==1], na.rm=T), xmax=b, 
                                ymax=max(simdata$y[simdata$tr==1],na.rm=T)), 
              color="black", size=0.5, alpha=0, fill="grey") +
    geom_point(aes(x=x, y=y, color=tr), size=2) + 
    ylim(c(min(simdata$y[simdata$tr==1]), max(simdata$y[simdata$tr==1]))) +
    xlim(c(min(simdata$x[simdata$tr==1]), max(simdata$x[simdata$tr==1]))) +
    scale_color_manual(values=c('red', 'black'), labels=c('control', 'treated')) +
    xlab("x") + ylab("y") + theme_bw() + 
    theme(legend.position = c(0.95,0.9), legend.background = element_rect(color="black", size=0.5),
          legend.title = element_blank())
}


# CATE_error_plot(res)
# dat <- matrix(c(1, 3, 2, 5, 2, 4, 0, 2, 1, 2.5, 1.5, 1.5), ncol = 6, byrow = TRUE) %>%
#   as.data.frame() %>%
#   `colnames<-`(c('A', 'B', 'C', 'D', 'E', 'G'))  %>%
#   mutate(id = c(1:2))
#   
#   p <- 
#     ggplot(data = dat) + 
#     geom_rect_interactive(aes(xmin = A, xmax = B, ymin = C, ymax = D,
#                               data_id = id),
#               fill = 'grey', color = 'black', 
#               size = 0.5, alpha = 0.3) + 
#     geom_point_interactive(aes(x = E, y = G, data_id = id)) 
#   
#   girafe(ggobj = p,
#          options = list(
#            opts_hover(
#            css = girafe_css(
#              css = 'fill:orange',
#              area = 'fill:red',
#              point = 'fill:green;r:5'
#          ))))
# 
# plot_bins <- function(simdata) {
#   dims <- length(simdata[[1]][['a']]) 
#   if (dims == 1) {
#     simdata %<>% 
#       lapply(function(x) {
#         x[['a']] <- c(x[['a']], min(y))
#         x[['b']] <- c(x[['b']], max(y))
#         x
#       })
#   }
#   
#   simdata <- 
#     sapply(1:n, function(i) {
#       x <- simdata[[i]]
#       y_est <- mean(y_test[-i][which(mip_outputs[[i]]$w >= 0.1)])
#       c(x[['a']][1], x[['a']][2], x[['b']][1], x[['b']][2], y_test[i], y_est, x_test[i])
#     }) %>%
#     as.data.frame() %>%
#     `colnames<-`(c('a1', 'a2', 'b1', 'b2', 'y', 'y_est', 
#                    paste0(rep('x', dims), 1:dims)))
#   
  # p <- 
  #   ggplot(simdata) + 
  #   geom_rect(aes(xmin = a1, xmax = b1, ymin = a2, ymax = b2),
  #             fill = 'grey', color = 'black')
    
  # ##### 2D
  # p <- 
  #   ggplot(simdata) + 
  #   geom_rect(aes(xmin = a1, ymin = a2, xmax = b1, ymax = b2), 
  #             fill = "grey", color = "black", 
  #             size = 0.5, alpha = 0.3) + 
  #   geom_point(aes(x = x1, y = x2, color = abs(yest - y) / mean(abs(y))), 
  #              size=2) + 
  #   scale_color_gradient(low = "blue", high = "red") + 
  #   xlab("x1") + ylab("x2") +  labs(color = "% Abs. error") +  
  #   theme_bw() + 
  #   theme(legend.position = c(0.9, 0.5), 
  #         legend.background = element_rect(color = "black", size = 0.5))
  #   
  ### 1D
#   ggplot(simdata) + 
#     geom_rect(aes(xmin=a1, ymin=min(y), xmax=b1, ymax=max(y)), 
#               color="black", size=0.5, alpha=0.3, fill="grey") +
#     geom_point(aes(x=x1, y=y), color='red', size=2) + 
#     geom_point(aes(x=x1, y=yest), color="blue") + 
#     geom_point(data=data.frame(x_train, y_train),aes(x=x_train, y=y_train), pch=18) +
#     xlab("x") + ylab("y") + theme_bw() + 
#     theme(legend.position = c(0.8,0.2), legend.background = element_rect(color="black", size=0.5),
#           legend.title = element_blank())
#   print(p)
# }







########################################################################
# unique_HTEs <- unique(HTE$actual)
# if (length(unique_HTEs) < nrow(HTE) / n_estimators) { # Constant treatment effect 
#   gg <- 
#     ggplot(HTE, aes(x = as.factor(actual), y = predicted, color = estimator)) + 
#     geom_boxplot()
#   for (i in 1:length(unique_HTEs)) {
#     gg <- gg + 
#       geom_hline(yintercept = unique_HTEs[i])
#   }
#   gg <- gg + 
#     labs(x = 'Actual', 
#          y = 'Predicted',
#          title = '(Piecewise-)Constant Treatment Effect')
#   print(gg)
# } else {
#   ggplot(HTE, aes(x = actual, y = predicted)) + 
#     geom_point(aes(color = estimator)) + 
#     geom_abline(intercept = 0, slope = 1, color = 'black') + 
#     labs(title = 'Heterogeneous Treatment Effect')
# }
# 
# partitions <- 
#   bins %>%
#   c() %>%
#   unique()
# 
# g <- ggplot(data = df, aes(x = X, y = Y)) + 
#   geom_point(aes(color = treated))
# for (partition in partitions) {
#   g <- g + geom_vline(xintercept = partition)
# }
# plot(g)