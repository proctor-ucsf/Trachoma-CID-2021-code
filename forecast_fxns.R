library("EnvStats") #for ppareto, dpareto
library("tidyverse")
#library("Renext")
#library("gridExtra")
#library("searcher")
IGN_STEP <- 0.5

# To add forecast regressions, modify: Init_par, pdf_fxn, setup_pdf, cdf_fxn

clean_rawdata <- function (raw_data, yrs = NULL, zero_thresh,key_cols, bin_width = 0, keep_na_surveys = TRUE, duplicates = FALSE) {
  # Remove duplicates.  (Some TF prevalences were assigned to groups of regions)
  if (keep_na_surveys == FALSE) {
    raw_data <- raw_data %>% filter(survey_type != "")
  }
  if (duplicates == FALSE) {
    if (keep_na_surveys == TRUE)
      raw_data <- distinct(raw_data,admin0_id, admin1_id, survey_year, tf_prev, .keep_all= TRUE)
    else
      raw_data <- distinct(raw_data,admin0_id, admin1_id, survey_year, tf_prev, survey_type, .keep_all= TRUE)
  }
  if (keep_na_surveys == FALSE)
  # Remove NA
  raw_data <- raw_data[is.finite(raw_data$tf_prev),]
  # Remove unncessary columns
  raw_data <- raw_data %>% dplyr::select(key_cols)  
  # Remove distrcits with essentially no TF
  raw_data <- raw_data %>% filter(tf_prev >= zero_thresh) # same as Pinset et al.
  # Select just relevant years and data values
  raw_data <- raw_data %>% filter(survey_year %in% yrs) # same as Pinset et al.
  # Indexing data into bins
  if (bin_width > 0) {
    raw_data <- cbind(raw_data,bin_num = floor((raw_data$tf_prev-zero_thresh)/bin_width+1)) #same as Pinset et al.
  }
  raw_data
}

logit <- function (p) {
  #log is taking during regression
  p/(1-p)
}

inv_logit <- function (a) {
  #exponentiation is taking post-regression
  a/(1+a)
}

det_init_par <- function (idata,fit_dist) {
  mean_data <- mean(idata$tf_prev)
  sd_data <- sd(idata$tf_prev)
  if (fit_dist == "renorm_exp") {
    init_par <- list(rate = 1/mean_data)
  } else if (fit_dist == "renorm_exp_trunc") {
    init_par <- list(rate = 1/mean_data)
  } else if (fit_dist == "renorm_gamma") {
    init_par <- list(shape = (mean_data/sd_data)^2, rate = mean_data/sd_data^2)
  } else if (fit_dist == "renorm_norm") {
    init_par <- list(mean = 0, sd = sd_data)
  } else if (fit_dist == "renorm_norm_lin") {
    init_par <- list(exp_mean = 1e-4, sd = sd_data)
  } else if (fit_dist == "renorm_beta") {
    beta_mean = mean_data/100
    beta_sd = sd_data/100
    shape1 <- ((1-beta_mean)/beta_sd^2 - 1/beta_mean)*beta_mean^2
    init_par <- list(shape1 = shape1, shape2 = shape1 * (1/beta_mean -1))
  } else if (fit_dist == "renorm_pareto") {
    # Approximation to avoid dealing with numeric solution
    init_par <- list(location = 2*mean_data/3, shape = 3)
  } else if (fit_dist == "renorm_lomax") {
    # Approximation to avoid dealing with numeric solution
    init_par <- list(scale = 2*mean_data/3, shape = 3)
  } else if (fit_dist == "moment_gamma") {
    init_par <- list(mean = mean_data, sd = sd_data)
  } else if (fit_dist == "moment_beta") {
    init_par <- list(mean = mean_data, sd = sd_data)
  } else if (fit_dist == "mlogit_beta") {
    init_par <- list(logit_mean = logit(mean_data/100), logit_sd = logit(sd_data/25)) #100 and 25 are the maximum avg and sd of the TF respectively
  } else if (fit_dist == "moment_pareto") {
    # adj_shape + 2 = shape.  This way, min(adj_shape) = 0 corresponds to Std Dev asymptomically approaching Inf and we can regress on logarithm
    init_par <- list(mean = mean_data, adj_shape = 1)
  } else if (fit_dist == "offset_exp") {
    init_par <- list(rate = 1/mean_data)
  } else if (fit_dist == "offset_gamma") {
    init_par <- list(mean = mean_data, sd = sd_data)
  } else if (fit_dist == "offset_beta") {
    init_par <- list(logit_mean = logit((mean_data)/100), logit_sd = logit(sd_data/25)) #100 and 25 are the maximum avg and sd of the TF respectively
  } else if (fit_dist == "offset_pareto") {
    init_par <- list(mean = mean_data, adj_shape = 1)
  } else if (fit_dist == "mlogit_beta_B") {
    init_par <- list(logit_mean = logit(mean_data/100), logit_sd = logit(sd_data/25)) #100 and 25 are the maximum avg and sd of the TF respectively
  } else if (fit_dist == "moment_pareto_B") {
    # adj_shape + 2 = shape.  This way, min(adj_shape) = 0 corresponds to Std Dev asymptomically approaching Inf and we can regress on logarithm
    init_par <- list(mean = mean_data, adj_shape = 1)
  } else if (fit_dist == "exp_norm_mix") {
    init_par <- list(mean_e = mean_data, mean_g = 35, sd_g = 10, prop_g = 0.1)
  }
  init_par
}

# Calculate PDF
pdf_fxn <- function (x,pdf_xval,fit_dist,zero_thresh) {
  if (any(is.na(x))) { 
    pdf_val <- 0
  } else if (fit_dist == "renorm_exp") {
    pdf_val <- dexp(pdf_xval,rate = x) / (pexp(100,rate = x) - pexp(zero_thresh, rate = x))
  } else if (fit_dist == "renorm_exp_trunc") {
    pdf_val <- dexp(pdf_xval,rate = x) / (pexp(EXP_THRESH,rate = x) - pexp(zero_thresh, rate = x))
  } else if (fit_dist == "renorm_gamma") {
    pdf_val <- dgamma(pdf_xval,shape = x[[1]], rate = x[[2]]) / 
      (pgamma(100,shape = x[[1]], rate = x[[2]]) - pgamma(zero_thresh,shape = x[[1]], rate = x[[2]]))
  } else if (fit_dist == "renorm_norm") {
    pdf_val <- dnorm(pdf_xval,mean = x[[1]], sd = x[[2]]) /
      (pnorm(100,mean = x[[1]], sd = x[[2]]) - pnorm(zero_thresh,mean = x[[1]], sd = x[[2]]))
  } else if (fit_dist == "renorm_norm_lin") {
    pdf_val <- dnorm(pdf_xval,mean = log(x[[1]]), sd = x[[2]]) /
      (pnorm(100,mean = log(x[[1]]), sd = x[[2]]) - pnorm(zero_thresh,mean = log(x[[1]]), sd = x[[2]]))
  } else if (fit_dist == "renorm_beta") {
    pdf_val <- dbeta(pdf_xval/100,shape1 = x[[1]], shape2 = x[[2]]) /
      (pbeta(100/100,shape1 = x[[1]], shape2 = x[[2]]) - pbeta(zero_thresh/100,shape1 = x[[1]], shape2 = x[[2]])) / 100
  } else if (fit_dist == "renorm_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      pdf_val <- dpareto(pdf_xval + x[[1]],location = x[[1]], shape = x[[2]]) /
        (ppareto(100 + x[[1]],location = x[[1]], shape = x[[2]]) - ppareto(zero_thresh + x[[1]],location = x[[1]], shape = x[[2]]))
    }
  } else if (fit_dist == "renorm_lomax") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      pdf_val <- dlomax(pdf_xval,scale = x[[1]], shape = x[[2]]) /
        (plomax(100,scale = x[[1]], shape = x[[2]]) - plomax(zero_thresh,scale = x[[1]], shape = x[[2]]))
    }
  } else if (fit_dist == "moment_gamma") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      sd <- x[[2]]
      shape <- (avg/sd)^2
      rate <- avg/sd^2
      pdf_val <- dgamma(pdf_xval,shape = shape, rate = rate) / 
        (pgamma(100,shape = shape, rate = rate) - pgamma(zero_thresh,shape = shape, rate = rate))
    }
  } else if (fit_dist == "moment_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]/100
      sd <- x[[2]]/100
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      pdf_val <- dbeta(pdf_xval/100,shape1 = shape1, shape2 = shape2) /
        (pbeta(100/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2)) / 100
    }
  } else if (fit_dist == "mlogit_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- inv_logit(x[[1]])
      sd <- inv_logit(x[[2]])*.25 # 0.25 is the max sd of beta
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      pdf_val <- dbeta(pdf_xval/100,shape1 = shape1, shape2 = shape2) /
        (pbeta(100/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2)) / 100
    }
  } else if (fit_dist == "moment_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      # Lowmax: Mean_l = scale / (shape-1)
      # Lomax (shape, scale) = Pareto(shape, scale) - scale
      # Pareto: Mean_p = shape*scale / (shape -1) [scale is same as location]
      #
      # mean_l + scale = mean_p = shape*scale/(shape-1)  
      # mean_l = scale/(shape-1) [aveg corresponds to lomax average since we want min value to be 0]
      avg <- x[[1]]
      shape <- x[[2]] + 2
      location <- (shape -1)*avg
      pdf_val <- dpareto(pdf_xval + location,location = location, shape = shape) /
        (ppareto(100 + location,location = location, shape = shape) - ppareto(zero_thresh + location,location = location, shape = shape))
    }
  } else if (fit_dist == "offset_exp") {
    pdf_val <- dexp(pdf_xval-zero_thresh,rate = x) / (pexp(100-zero_thresh,rate = x))
  } else if (fit_dist == "offset_gamma") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      sd <- x[[2]]
      shape <- (avg/sd)^2
      rate <- avg/sd^2
      pdf_xval<-pdf_xval[pdf_xval>zero_thresh]
      pdf_val <- dgamma(pdf_xval-zero_thresh,shape = shape, rate = rate) / 
        (pgamma(100-zero_thresh,shape = shape, rate = rate))
    }
  } else if (fit_dist == "offset_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- inv_logit(x[[1]])
      sd <- inv_logit(x[[2]])*.25 # 0.25 is the max sd of beta
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      pdf_xval<-pdf_xval[pdf_xval>zero_thresh]
      pdf_val <- dbeta((pdf_xval-zero_thresh)/(100-zero_thresh),shape1 = shape1, shape2 = shape2) /
        (pbeta((100-zero_thresh)/(100-zero_thresh),shape1 = shape1, shape2 = shape2)) / (100-zero_thresh)
    }
  } else if (fit_dist == "offset_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      shape <- x[[2]] + 2
      location <- (shape -1)*avg
      pdf_xval<-pdf_xval[pdf_xval>zero_thresh]
      pdf_val <- dpareto(pdf_xval-zero_thresh + location,location = location, shape = shape) /
        (ppareto(100-zero_thresh + location,location = location, shape = shape))
    }
  } else if (fit_dist == "exp_norm_mix") {
    if(x[[2]] <= 0) {pdf_val <- 0}
    else {
      p_g <- x[[4]]
      pdf_val <- (dexp(pdf_xval, rate = 1/x[[1]])*(1-p_g) + dnorm(pdf_xval,mean = x[[2]],sd = x[[3]])*p_g) /
        (pexp(100, rate = 1/x[[1]])*(1-p_g) - pexp(zero_thresh, rate = 1/x[[1]])*(1-p_g) + 
           pnorm(100,mean = x[[2]],sd = x[[3]])*p_g - pnorm(zero_thresh,mean = x[[2]],sd = x[[3]])*p_g) 
    }
  }
}

# Dummy function to just help pass parameter functions around
setup_for_pdf <- function(par_data) {
  x<-NA
  if (par_data$fit[1] == "renorm_exp") {
    x <- as.numeric(par_data[par_data$par_name == "rate","par_value"])
  } else if (par_data$fit[1] == "renorm_exp_trunc") {
    x <- as.numeric(par_data[par_data$par_name == "rate","par_value"])
  } else if (par_data$fit[1] == "renorm_gamma") {
    x <- as.numeric(c(par_data[par_data$par_name == "shape","par_value"],par_data[par_data$par_name == "rate","par_value"]))
  } else if (par_data$fit[1] == "renorm_norm") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "sd","par_value"]))
  } else if (par_data$fit[1] == "renorm_norm_lin") {
    x <- as.numeric(c(par_data[par_data$par_name == "exp_mean","par_value"],par_data[par_data$par_name == "sd","par_value"]))
  } else if (par_data$fit[1] == "renorm_beta") {
    x <- as.numeric(c(par_data[par_data$par_name == "shape1","par_value"],par_data[par_data$par_name == "shape2","par_value"]))
  } else if (par_data$fit[1] == "renorm_pareto") {
    x <- as.numeric(c(par_data[par_data$par_name == "location","par_value"],par_data[par_data$par_name == "shape","par_value"]))
  } else if (par_data$fit[1] == "renorm_lomax") {
    x <- as.numeric(c(par_data[par_data$par_name == "scale","par_value"],par_data[par_data$par_name == "shape","par_value"]))
  } else if (par_data$fit[1] == "moment_gamma") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "sd","par_value"]))
  } else if (par_data$fit[1] == "moment_beta") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "sd","par_value"]))
  } else if (par_data$fit[1] == "mlogit_beta") {
    x <- as.numeric(c(par_data[par_data$par_name == "logit_mean","par_value"],par_data[par_data$par_name == "logit_sd","par_value"]))
  } else if (par_data$fit[1] == "moment_pareto") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "adj_shape","par_value"]))
  } else if (par_data$fit[1] == "offset_exp") {
    x <- as.numeric(par_data[par_data$par_name == "rate","par_value"])
  } else if (par_data$fit[1] == "offset_gamma") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "sd","par_value"]))
  } else if (par_data$fit[1] == "offset_beta") {
    x <- as.numeric(c(par_data[par_data$par_name == "logit_mean","par_value"],par_data[par_data$par_name == "logit_sd","par_value"]))
  } else if (par_data$fit[1] == "offset_pareto") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean","par_value"],par_data[par_data$par_name == "adj_shape","par_value"]))
  } else if (par_data$fit[1] == "exp_norm_mix") {
    x <- as.numeric(c(par_data[par_data$par_name == "mean_e","par_value"],par_data[par_data$par_name == "mean_g","par_value"],
                      par_data[par_data$par_name == "sd_g","par_value"],par_data[par_data$par_name == "prop_g","par_value"]))
    }
  x
}

# Calculate CDF
cdf_fxn <- function (x,cdf_xval,fit_dist, zero_thresh) {
  if (fit_dist == "renorm_exp") {
    cdf_val <- (pexp(cdf_xval,rate = x) - pexp(zero_thresh,rate = x)) / (pexp(100,rate = x) - pexp(zero_thresh, rate = x))
  } else if (fit_dist == "renorm_gamma") {
    cdf_val <- (pgamma(cdf_xval,shape = x[[1]], rate = x[[2]]) - pgamma(zero_thresh,shape = x[[1]], rate = x[[2]])) / 
      (pgamma(100,shape = x[[1]], rate = x[[2]]) - pgamma(zero_thresh,shape = x[[1]], rate = x[[2]]))
  } else if (fit_dist == "renorm_norm") {
    cdf_val <- (pnorm(cdf_xval,mean = x[[1]], sd = x[[2]]) - pnorm(zero_thresh,mean = x[[1]], sd = x[[2]])) /
      (pnorm(100,mean = x[[1]], sd = x[[2]]) - pnorm(zero_thresh,mean = x[[1]], sd = x[[2]]))
  } else if (fit_dist == "renorm_norm_lin") {
    cdf_val <- (pnorm(cdf_xval,mean = log(x[[1]]), sd = x[[2]]) - pnorm(zero_thresh,mean = log(x[[1]]), sd = x[[2]])) /
      (pnorm(100,mean = log(x[[1]]), sd = x[[2]]) - pnorm(zero_thresh,mean = log(x[[1]]), sd = x[[2]]))
  } else if (fit_dist == "renorm_beta") {
    cdf_val <- (pbeta(cdf_xval/100,shape1 = x[[1]], shape2 = x[[2]]) - pbeta(zero_thresh/100,shape1 = x[[1]], shape2 = x[[2]])) /
      (pbeta(100/100,shape1 = x[[1]], shape2 = x[[2]]) - pbeta(zero_thresh/100,shape1 = x[[1]], shape2 = x[[2]]))
  } else if (fit_dist == "renorm_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {cdf_val <- 0}
    else {
      cdf_val <- (ppareto(cdf_xval + x[[1]],location = x[[1]], shape = x[[2]]) - ppareto(zero_thresh + x[[1]],location = x[[1]], shape = x[[2]])) /
        (ppareto(100 + x[[1]],location = x[[1]], shape = x[[2]]) - ppareto(zero_thresh + x[[1]],location = x[[1]], shape = x[[2]]))
    }
  } else if (fit_dist == "renorm_lomax") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {cdf_val <- 0}
    else {
      cdf_val <- (plomax(cdf_xval,scale = x[[1]], shape = x[[2]]) - plomax(zero_thresh,scale = x[[1]], shape = x[[2]])) /
        (plomax(100,scale = x[[1]], shape = x[[2]]) - plomax(zero_thresh,scale = x[[1]], shape = x[[2]]))
    }
  } else if (fit_dist == "moment_gamma") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      sd <- x[[2]]
      shape <- (avg/sd)^2
      rate <- avg/sd^2
      cdf_val <- (pgamma(cdf_xval,shape = shape, rate = rate) - pgamma(zero_thresh,shape = shape, rate = rate)) / 
        (pgamma(100,shape = shape, rate = rate) - pgamma(zero_thresh,shape = shape, rate = rate))
    }
  } else if (fit_dist == "moment_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {cdf_val <- 0}
    else {
      avg <- x[[1]]/100
      sd <- x[[2]]/100
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      cdf_val <- (pbeta(cdf_xval/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2)) /
        (pbeta(100/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2))
    }
  } else if (fit_dist == "mlogit_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {cdf_val <- 0}
    else {
      avg <- inv_logit(x[[1]])
      sd <- inv_logit(x[[2]]) *.25 # 0.25 is max sd of n.b distribution
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      cdf_val <- (pbeta(cdf_xval/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2)) /
        (pbeta(100/100,shape1 = shape1, shape2 = shape2) - pbeta(zero_thresh/100,shape1 = shape1, shape2 = shape2))
    }
  } else if (fit_dist == "moment_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {cdf_val <- 0}
    else {
      avg <- x[[1]]
      shape <- x[[2]] + 2
      location <- (shape - 1)*avg
      cdf_val <- (ppareto(cdf_xval + location,location = location, shape = shape) - ppareto(zero_thresh + location,location = location, shape = shape)) /
        (ppareto(100 + location,location = location, shape = shape) - ppareto(zero_thresh + location,location = location, shape = shape))
    }
  } else if (fit_dist == "offset_exp") {
    cdf_val <- pexp(cdf_xval-zero_thresh,rate = x) / (pexp(100-zero_thresh,rate = x))
  } else if (fit_dist == "offset_gamma") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      sd <- x[[2]]
      shape <- (avg/sd)^2
      rate <- avg/sd^2
      cdf_val <- pgamma(cdf_xval-zero_thresh,shape = shape, rate = rate) / 
        (pgamma(100-zero_thresh,shape = shape, rate = rate))
    }
  } else if (fit_dist == "offset_beta") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- inv_logit(x[[1]])
      sd <- inv_logit(x[[2]])*.25 # 0.25 is the max sd of beta
      shape1 <- ((1-avg)/sd^2 - 1/avg)*avg^2
      shape2 <-  shape1 * (1/avg -1)
      cdf_val <- pbeta((cdf_xval-zero_thresh)/(100-zero_thresh),shape1 = shape1, shape2 = shape2) /
        (pbeta((100-zero_thresh)/(100-zero_thresh),shape1 = shape1, shape2 = shape2)) / (100-zero_thresh)
    }
  } else if (fit_dist == "offset_pareto") {
    if(x[[1]] <= 0 | x[[2]] <= 0) {pdf_val <- 0}
    else {
      avg <- x[[1]]
      shape <- x[[2]] + 2
      location <- (shape -1)*avg
      cdf_val <- ppareto(cdf_xval-zero_thresh + location,location = location, shape = shape) /
        (ppareto(100-zero_thresh + location,location = location, shape = shape))
    }
  } else if (fit_dist == "exp_norm_mix") {
    p_g <- x[[4]]
    cdf_val <- (pexp(cdf_xval, rate = 1/x[[1]])*(1-p_g) + pnorm(cdf_xval,mean = x[[2]],sd = x[[3]])*p_g) /
      (pexp(100, rate = 1/x[[1]])*(1-p_g) - pexp(zero_thresh, rate = 1/x[[1]])*(1-p_g) + 
         pnorm(100,mean = x[[2]],sd = x[[3]])*p_g - pnorm(zero_thresh,mean = x[[2]],sd = x[[3]])*p_g)
  }
  #  cdf_val[!is.finite(cdf_val) | is.na(cdf_val)] <- 1e-300
  cdf_val
}


# Generate binned pdfs
pdf_gen_bindata <- function (model_data, zero_thresh = ZERO_THRESH, bin_width = BIN_WIDTH) {
  xvar = seq(zero_thresh + 0.5*bin_width, 100 -0.5*bin_width, bin_width)
  pdf_data<-model_data %>%
    group_by(fit,year) %>%
    do({
      x <- setup_for_pdf(.)
      if (any(is.na(x))) {
        data_frame()
      } else {
        delta = (max(xvar)-min(xvar))/(length(xvar)-1)/2
        yvar_pdf <- cdf_fxn(x, cdf_xval =  xvar + delta, fit_dist = .$fit[1], zero_thresh = zero_thresh) -
          cdf_fxn(x, cdf_xval =  xvar - delta, fit_dist = .$fit[1], zero_thresh = zero_thresh)
        if(any(!is.finite(yvar_pdf))) {
          data_frame()
        } else {
          data_frame(fit = .$fit[1], year = .$year[1], xvar = xvar, yvar_pdf = yvar_pdf)
        }
      }
    })
}

calc_logL_pdf <- function(x,focal_data,model,weighted=FALSE,zero_thresh=ZERO_THRESH){
  num_var <- length(x)/2
  num_years <- length(unique(focal_data$survey_year))
  x<-as.double(x)
  L_scores <- focal_data %>% group_by(survey_year) %>% do({
    val_year <- exp(x[1:num_var] + (.$survey_year[1]-2020)*x[(num_var+1):(2*num_var)])
    if (weighted == TRUE) {
      logL <- log(pdf_fxn(val_year,.$tf_prev,model,zero_thresh))/nrow(.)
    } else {
      logL <- log(pdf_fxn(val_year,.$tf_prev,model,zero_thresh))
    }
    data_frame(logL = logL)
  })
  logL <- sum(L_scores$logL)
}

mle_fit <- function (data,fits) {
  data_frame (fit = fits) %>% group_by(fit) %>% do({
    # Initialize optimizer
    print(.$fit[1])
    init_par <- det_init_par(data,.$fit[1])
    init_val <- c(as.double(init_par), rep(1,length(init_par)))
    # Optimize - note the log of the coefficients are passed along
    optim_res <- optim(init_val, function (x) - calc_logL_pdf(log(x),data,.$fit[1]))
    opt_par <- log(optim_res$par)
    num_val <- length(opt_par)/2
    data_frame(par_name = rep(names(init_par), 2), coef = rep(c("intercept","slope"),each = num_val), par_value = opt_par)
  })
}

det_coef <- function (mle_input, yrs = CLEAN_YRS) {
  mle_input %>% group_by(fit,par_name) %>% do({
    intercept <- .[.$coef == "intercept","par_value"]
    slope <- .[.$coef == "slope","par_value"]
    data_frame(year = yrs, par_value = exp(intercept$par_value + (yrs-2020)*slope$par_value))
  })
}

# Calculate log likelihood of a model (including time variation)

det_mle_ignorance <- function (idata,fit_coef) {
  fit_coef %>% filter(year == idata$survey_year[1]) %>% group_by(fit) %>% do({
    # Initialize optimizer
    par_val <- setup_for_pdf(.)
    # Optimize - note the log of the coefficients are passed along
    optim_res <- optimize(function(x) - calc_logL_bin_1yr(par_val,idata,.$fit[1], ignorance = x),interval = c(-1,10))
    mle_ignorance <- optim_res$minimum
    data_frame(mle_ignorance = mle_ignorance)
  })
}

calc_pdf_bin <- function(par_val, fit, ignorance = 0, zero_thresh = ZERO_THRESH, pdf_width = PDF_WIDTH) {
#  browser()
  xvar = seq(zero_thresh + 0.5*pdf_width, 100 -0.5*pdf_width, pdf_width)
  delta = (max(xvar)-min(xvar))/(length(xvar)-1)/2
  yvar_pdf <- cdf_fxn(par_val, cdf_xval =  xvar + delta, fit_dist = fit, zero_thresh = zero_thresh) -
    cdf_fxn(par_val, cdf_xval =  xvar - delta, fit_dist = fit, zero_thresh = zero_thresh)
  
  # Adjust for ignorance
  if (ignorance != 0) {
    yvar_pdf <- yvar_pdf^(1+ignorance)
    yvar_pdf <- yvar_pdf / sum (yvar_pdf)
  }
  data_frame(xvar = xvar, yvar_pdf = yvar_pdf)
}

calc_logL_bin_1yr <- function (par_val, idata, fit, ignorance) {
  pdf_tbl <- calc_pdf_bin(par_val = par_val, fit = fit, ignorance = ignorance)
  # calculate pdf
  pdf_table <- data.frame(bin_num = 1:nrow(pdf_tbl), xvar = pdf_tbl$xvar, logL = log(pdf_tbl$yvar_pdf))
  idata <- left_join(idata,pdf_table, by = "bin_num")
  
  logL <- sum(idata$logL)
}

# ################## TEMP - bin evrything in unit percentiles for likelihood calculation.  Much slower
# # determine mle fits
# mle_fit_bin <- function (data,fits) {
#   data_frame (fits = fits) %>% group_by(fits) %>% do({
#     # Initialize optimizer
#     print(.$fits[1])
#     init_par <- det_init_par(data,.$fits[1])
#     init_val <- c(as.double(init_par), rep(1,length(init_par)))
#     # Optimize - note the log of the coefficients are passed along
#     optim_res <- optim(init_val, function (x) - calc_logL_bin(log(x),data,.$fits[1]))
#     opt_par <- log(optim_res$par)
#     num_val <- length(opt_par)/2
#     data_frame(fit = .$fits[1],par_name = rep(names(init_par), num_val), coef = rep(c("intercept","slope"),each = num_val), par_value = opt_par)
#   })
# }
# 
# ################## TEMP - bin evrything in unit percentiles for likelihood calculation.  Much slower
# calc_logL_bin <- function(x,focal_data,model,zero_thresh=ZERO_THRESH,bin_width=BIN_WIDTH){
#   num_var <- length(x)/2
#   num_years <- length(unique(focal_data$survey_year))
#   x<-as.double(x)
#   if(any(is.nan(x))) {
#     logL <- Inf
#   } else {
#     L_scores <- focal_data %>% group_by(survey_year) %>% do({
#       val_year <- exp(x[1:num_var] + (.$survey_year[1]-2020)*x[(num_var+1):(2*num_var)])
# 
#       xvar = seq(zero_thresh + 0.5*bin_width, 100 -0.5*bin_width, bin_width)
#       delta = (max(xvar)-min(xvar))/(length(xvar)-1)/2
#       yvar_pdf <- cdf_fxn(val_year, cdf_xval =  xvar + delta, fit_dist = model, zero_thresh = zero_thresh) -
#         cdf_fxn(val_year, cdf_xval =  xvar - delta, fit_dist = model, zero_thresh = zero_thresh)
#       pdf_model <- data_frame(fit = model, year = .$survey_year[1], xvar = xvar, yvar_pdf = yvar_pdf)
# 
#       # Score fits
#       if (nrow(pdf_model) > 0 ){
#         score <- score_by_bin(.,pdf_model)
#       }
#       data_frame(logL = score$score)
#     })
#     logL <- sum(L_scores$logL)
#   }
# }
# 
# # Score data according to bins
# ################## TEMP - bin evrything in unit percentiles for likelihood calculation.  Much slower
# score_by_bin <- function(data,model_pdf,ignorance = 0, zero_thresh = ZERO_THRESH, bin_width = BIN_WIDTH) {
#   data$bin_num <- floor((data$tf_prev-zero_thresh)/bin_width) * bin_width + zero_thresh + bin_width/2
#   model_pdf$xvar <- round(model_pdf$xvar*1e6)/1e6 # Weirdly there can be some rounding errors
#   data$bin_num <- round(data$bin_num*1e6)/1e6 # Weirdly there can be some rounding errors
#   if(length(unique(model_pdf$year)) > 1) {
#     print("POTENTIAL BUG:  score_by_bin asked to score data containing more than one year")
#   }
#   model_pdf %>% group_by(fit) %>% do ({
#     counts <-count(data[data$survey_year == .$year[1],], vars = bin_num)
#     if (ignorance != 0) {
#       .$yvar_pdf <- .$yvar_pdf^(1+ignorance)
#       .$yvar_pdf <- .$yvar_pdf / sum (.$yvar_pdf)
#     }
#     scorecard <- left_join(counts, ., by = c("vars" = "xvar")) -> tt
#     score <- sum(scorecard$n*log(scorecard$yvar_pdf/bin_width))
#     data_frame(fit = .$fit[1], score = score)
#   })
# }
# 
# # Determine point-pdf values for data
# pdf_gen_pointdata <- function (data, model_pars, zero_thresh) {
#   print("XXX pdf_gen_pointdata")
#   data<-data[data$zero_bin == FALSE,]
#   pdf_frame<-model_pars %>%
#     group_by(fit,year) %>%
#     do({
#       x <- setup_for_pdf(.)
#       xvar <- data[data$survey_year == .$year[1],"tf_prev"]
#       if (any(is.na(x)) | length(xvar) == 0) {
#         data_frame()
#       } else {
#         yvar_pdf <- pdf_fxn(x, pdf_xval =  xvar, fit_dist = .$fit[1], zero_thresh = zero_thresh)
#         if(any(!is.finite(yvar_pdf))){
#           data_frame()
#         } else {
#           data_frame(fit = .$fit[1], year = .$year[1], xvar = xvar, yvar_pdf = yvar_pdf)
#         }
#       }
#     })
# }
# 
# # Score data using point-value of pdf
# score_by_pointpdf <- function(data, model_pars, zero_thresh) {
#   print("XXX score_by_pointpdf")
#   pdf_pointdata <- pdf_gen_pointdata(data = data, model_pars = model_pars, zero_thresh = zero_thresh)
#   pdf_pointdata %>% group_by(fit, year) %>% 
#     summarize(score = sum(log(yvar_pdf)))
# }
# 
# 

# # Best fit based on direct optimization of likelihood
# optimize_fit <- function (data, distributions) {
#   focal_data <- data[data$zero_bin == FALSE,]
#   data_frame(fit = distributions) %>% group_by(fit) %>% do({ 
#     # Determine initial
#     init_par <- det_init_par(focal_data,.$fit[1])
#     init_val <- c(as.double(init_par), rep(1,length(init_par)))
#     # Optimize
#     optim_res <- optim(init_val, function (x) - calc_logL_bin(log(x),focal_data,.$fit[1]))
#     opt_par <- log(optim_res$par)
#     num_val <- length(opt_par)/2
#     data_frame(par_name = rep(names(init_par), num_val), coef = rep(c("intercept","slope"),each = num_val), par_value = opt_par, logL = -optim_res$value)
#   })
# }

# # Model creation based on direct optimization of temporal likelihood
# optimize_likelihood <- function (data,ipars,weighted=FALSE) {
#   ipars %>% group_by(forecast_len, num_prior_years, forecast_year, fit) %>% do({
#     train_yrs <- (.$forecast_year[1]-.$forecast_len[1]-.$num_prior_years[1]+1):(.$forecast_year[1]-.$forecast_len[1])
#     # Determine initial
#     focal_data <- data[data$survey_year %in% train_yrs,]
#     focal_data <- focal_data[focal_data$zero_bin == FALSE,]
#     init_par <- det_init_par(focal_data,.$fit)
#     ####################################
#     # Use next blox if optimizing log of coeff is wanted
#     ####################################
#     # init_val <- c(log(as.double(init_par)), rep(0,length(init_par))) # values are parameters of log_regression
#     # # Optimize
#     # optim_res <- optim(init_val, function (x) - calc_logL(x,focal_data,.$fit,weighted))
#     # opt_par <- optim_res$par
#     ####################################
#     init_val <- c(as.double(init_par), rep(1,length(init_par)))
#     # Optimize
#     optim_res <- optim(init_val, function (x) - calc_logL_bin(log(x),focal_data,.$fit[1]))
#     opt_par <- log(optim_res$par)
#     num_val <- length(opt_par)/2
#     if(.$num_prior_years == 1) {
#       # to prevent runaway regression when only one year is used for training
#       opt_par[1:num_val] <- opt_par[1:num_val] + (.$forecast_year[1]-.$forecast_len[1]-2020)*opt_par[num_val+(1:num_val)]
#       opt_par[num_val+(1:num_val)] <- 0
#     }
#     data_frame(par_name = rep(names(init_par), num_val), coef = rep(c("intercept","slope"),each = num_val), par_value = opt_par)
#   })
# }

# # Using direct likelihood calculation, compute forecasts, looping through: each model, each forecast interval, each length of prior data interval and each possible forecast year
# compute_forecast_byL <- function (data, weighted, forecast_years = NULL, models = DISTRIBUTIONS,forecast_len_arr = FORECAST_LEN, num_prior_years_arr = NUM_PRIOR_YEARS, clean_yrs = CLEAN_YRS) {
#   forecast_setup <- 
#     data_frame(forecast_len = rep(forecast_len_arr, each = length (num_prior_years_arr)), num_prior_years = rep(num_prior_years_arr, time = length(forecast_len_arr))) %>% 
#     group_by(forecast_len, num_prior_years) %>% do ({
#       # To allow for consistent comparison of scores, keeping the range of forecast years the same throughout
#       if(is.null(forecast_years)) {
#         forecast_years = (min(clean_yrs) + max(forecast_len_arr) + max(num_prior_years_arr) -1):max(clean_yrs)
#       }
#       data_frame(forecast_year = rep(forecast_years, each = length(models)), fit = rep(models, time = length(forecast_years)))
#     })
#   # Perform regression on models
#   regress_pars <- optimize_likelihood(data, forecast_setup,  weighted)
#   
# }

# # Calculate the pdf for a forecast determined from the method of direct likelihood optimization
# calc_forecast_pdf_L <- function(input, zero_thresh = ZERO_THRESH, bin_width = BIN_WIDTH) {
#   intercepts <- input[input$coef == "intercept",]
#   slopes <- input[input$coef == "slope",]
#   forecast_par <- exp(intercepts$par_value + (input$forecast_year[1]-2020)*slopes$par_value)
#   
#   if (any(is.na(forecast_par))) {
#     data_frame()
#   } else {
#     xvar = seq(zero_thresh + 0.5*bin_width, 100 -0.5*bin_width, bin_width)
#     delta = (max(xvar)-min(xvar))/(length(xvar)-1)/2
#     yvar_pdf <- cdf_fxn(forecast_par, cdf_xval =  xvar + delta, fit_dist = input$fit[1], zero_thresh = zero_thresh) -
#       cdf_fxn(forecast_par, cdf_xval =  xvar - delta, fit_dist = input$fit[1], zero_thresh = zero_thresh)
#     if(any(!is.finite(yvar_pdf))) {
#       data_frame()
#     } else {
#       #          yvar_pdf[yvar_pdf == 0] <- 1e-300
#       data_frame(fit = input$fit[1], year = input$forecast_year[1], xvar = xvar, yvar_pdf = yvar_pdf)
#     }
#   }
# }

# # For direct likelihood optimzation, compute forecasts, looping through: each model, each forecast interval, each length of prior data interval and each possible forecast year
# score_forecast_L <- function (forecast, data, ignorance = 0) {
#   forecast_res <- forecast %>% group_by (forecast_len, num_prior_years, forecast_year, fit) %>% do ({
#     pdf_model<-calc_forecast_pdf_L(.)
#     # Score fits
#     if (nrow(pdf_model) > 0 ){
#       score <- score_by_bin(data,pdf_model,ignorance = ignorance)
#       score_addition <- pdf_stats(pdf_model)
#       right_join(score,score_addition)
#     } else {
#       data_frame()
#     }
#   })
#   forecast_res
# }

# score_sum_by_coef <- function (forecast, data, ignorance, zero_thresh, bin_width) {
#   pdf_data_model <- pdf_gen_bindata(model_data = forecast, zero_thresh = zero_thresh, bin_width = bin_width)
#   # Score fits
#   scores <- score_by_bin(data,pdf_data_model,bin_width = bin_width,ignorance = ignorance) # insert years
#   sum(scores$score)
# }
# 
# forecast_ignorance_adj_by_coef <- function (par_arr,forecast, data, bin_width = BIN_WIDTH, zero_thresh = ZERO_THRESH) {
#   forecast %>% group_by(forecast_len, num_prior_years, fit) %>% do({
#     opt_par <- optimize(function (x) - score_sum_by_coef(.,data,ignorance = x, zero_thresh = zero_thresh, bin_width = bin_width), interval = c(-1,10))
#     alpha <- opt_par$minimum
#     pdf_data_model <- pdf_gen_bindata(model_data = ., zero_thresh = zero_thresh, bin_width = bin_width)
#     # Score fits
#     scores <- score_by_bin(data,pdf_data_model,bin_width = bin_width,ignorance = alpha) # insert years
#     scores$alpha <- alpha
#     scores$forecast_year <- scores$year
#     scores
#   })
# }
# 
# score_sum <- function (forecast, data, ignorance) {
#   scores <- score_forecast_L(forecast = forecast ,data = data,ignorance = ignorance)
#   sum(scores$score)
# }
# 
# forecast_ignorance_adj <- function (forecast, data) {
#   forecast %>% group_by(forecast_len, num_prior_years, fit) %>% do({
#     opt_par <- optimize(function (x) - score_sum(.,data,ignorance = x), interval = c(-1,10))
#     alpha <- opt_par$minimum
#     scores <- score_forecast_L(forecast = .,data = data,ignorance = alpha)
#     scores$alpha <- alpha
#     scores
#   })
# }
