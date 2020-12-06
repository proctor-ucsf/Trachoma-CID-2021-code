rm (list = ls())
source("0 - Config.R")
raw_data <- read.csv("../UCSF_DataUpdate_21Jan2020_get2020_prevalence_tab.csv")

####################################################
# Clean data
####################################################
clean_data <- clean_rawdata(raw_data, yrs = CLEAN_YRS, zero_thresh = ZERO_THRESH, key_cols = KEY_COL_NAMES)
write.csv(clean_data,"Output/clean_data.csv",row.names = FALSE)
clean_data <- clean_data %>% mutate(year = survey_year)

####################################################
# Determine best fit pdfs for overall data
####################################################
if (file.exists("Output/global_fit.csv")) {
  global_fit <- read.csv("Output/global_fit.csv")
} else {
  global_fit <- mle_fit(data = clean_data, fits = DISTRIBUTIONS)
  write.csv(global_fit, "Output/global_fit.csv")
}

fit_coef <- det_coef(global_fit)
fit_pdfs <- pdf_gen_bindata(fit_coef)

#normalize to number in year
year_count<-clean_data %>% group_by(year) %>% summarize(count = n())
fit_pdfs<- right_join(fit_pdfs,year_count, by = "year")

####################################################
# Determine best fit pdfs for yearly data
####################################################
if (file.exists("Output/fit_coef_yr.csv")) {
  fit_coef_yr <- read.csv("Output/fit_coef_yr.csv")
} else {
  fit_coef_yr <- clean_data %>% group_by(year) %>% do({
    yr_fit <- mle_fit(data = ., fits = DISTRIBUTIONS)
    det_coef(yr_fit,yrs = .$year[1])
  })
  write.csv(fit_coef_yr, "Output/fit_coef_yr.csv")
}

fit_pdfs_yr <- pdf_gen_bindata(fit_coef_yr)
fit_pdfs_yr<- right_join(fit_pdfs_yr,year_count, by = "year")

####################################################
# Determine best fit pdfs for year cohorts
####################################################
clean_data <- clean_data %>% mutate(cohort = floor((year - min(CLEAN_YRS))/COHORT_SIZE))
if (file.exists("Output/fit_coef_cohort.csv")) {
  fit_coef_cohort <- read.csv("Output/fit_coef_cohort.csv")
  fit_coef_cohort_sh <- read.csv("Output/fit_coef_cohort_sh.csv")
} else {
  fit_coef_cohort_sh <- clean_data %>% group_by(cohort) %>% do({
    yr_fit <- mle_fit(data = ., fits = DISTRIBUTIONS)
    det_coef(yr_fit,yrs = .$year[1]) %>% select(-year)
  })
  fit_coef_cohort <- tibble(year = CLEAN_YRS) %>% group_by(year) %>% do({
    c <- floor((.$year[1] - min(CLEAN_YRS))/COHORT_SIZE)
    fit_coef_cohort_sh %>% filter(cohort == c)
  })
  fit_coef_cohort_sh <- fit_coef_cohort_sh %>% mutate(year = cohort)
  write.csv(fit_coef_cohort_sh, "Output/fit_coef_cohort_sh.csv")
  write.csv(fit_coef_cohort, "Output/fit_coef_cohort.csv")
}
cohort_count<-clean_data %>% group_by(cohort) %>% summarize(year = min(cohort), count = n())
fit_pdfs_cohort_sh <- pdf_gen_bindata(fit_coef_cohort_sh)
fit_pdfs_cohort_sh <- right_join(fit_pdfs_cohort_sh,cohort_count, by = "year")


####################################################
# Plot data and best fit
####################################################

# Plot pdf + raw data
plot_fits <- function(pdfs, coefs, string) {
  ggplot () +
    geom_histogram(data = clean_data, aes(x = tf_prev),binwidth = DISPLAY_WIDTH) +
    #geom_annotate() +
    #geom_freqpoly(data = clean_data, aes(x = tf_prev),binwidth = 5) +
    #geom_line(data = pdfs, aes(x = xvar, y = yvar_pdf*count*DISPLAY_WIDTH, col = fit)) +
    #ggtitle(paste("Data fits: ",string, sep = "")) +
    facet_wrap(~ year, nrow = 4) +
    xlab("TF prevalence (%)") + ylab("Count")
  
  ggsave(paste("Figs/Fig_data_",string,".jpg",sep=''))
  
#  ggplot(coefs %>% filter(fit == 'exp_norm_mix', par_name == 'prop_g')) +
#    geom_line(aes(x=year, y = par_value))
  
  ggplot(coefs %>% filter(fit == 'exp_norm_mix')) +
    geom_line(aes(x=year, y = par_value)) +
    scale_y_continuous(trans='log10') +
    ggtitle(paste("Mixture coef: ",string, sep = "")) +
    facet_wrap(~par_name)
  ggsave(paste("Figs/Fig_mixturecoef_",string,".jpg",sep=''))
  
  ggplot(coefs %>% filter(par_name == 'prop_g')) +
    geom_line(aes(x=year, y = par_value)) +
    ggtitle(paste("Gaussian proportion: ",string, sep = "")) 
  ggsave(paste("Figs/Fig_propg_",string,".jpg",sep=''))
}
plot_fits(fit_pdfs,fit_coef,'all')
plot_fits(fit_pdfs_yr,fit_coef_yr,'yearly')
clean_data <- clean_data %>% mutate (year = cohort)
plot_fits(fit_pdfs_cohort_sh,fit_coef_cohort_sh,'cohort')
clean_data <- clean_data %>% mutate (year = survey_year)

####################################################
# Q-value
####################################################
q_eval <- function(pdfs,string) {
  delta <- pdfs$xvar[2] - pdfs$xvar[1]
  qval_res <- pdfs %>% filter(fit == "renorm_exp") %>% group_by(year,xvar) %>% do({
    x <- .$xvar[1]
    yr <- .$year[1]
    data_year <- clean_data %>% filter(year == yr)
    data_cdf <- data_year %>% filter(tf_prev >= x - delta/2)
    pdf_vals <- pdfs %>% filter(fit == 'renorm_exp', year == yr,xvar >= x)
    tibble(theor_cdf = sum(pdf_vals$yvar_pdf), obs_cdf = nrow(data_cdf)/nrow(data_year))
  })
  
  ggplot(qval_res) +
    geom_line(aes(x = xvar, y = (obs_cdf-theor_cdf)/theor_cdf)) +
    ggtitle(paste("Q value: ",string, sep = "")) +
    facet_wrap(~year)
  ggsave(paste("Figs/Fig_qval_",string,".jpg",sep=''))
  
  ggplot(qval_res) +
    geom_line(aes(x = xvar, y = obs_cdf-theor_cdf)) +
    ggtitle(paste("Excess above geometric: ",string, sep = "")) +
    facet_wrap(~year)
  ggsave(paste("Figs/Fig_excess_",string,".jpg",sep=''))
}
q_eval(fit_pdfs,'all')
q_eval(fit_pdfs_yr,'yearly')
clean_data <- clean_data %>% mutate (year = cohort)
q_eval(fit_pdfs_cohort_sh,'cohort')
clean_data <- clean_data %>% mutate (year = survey_year)

ggplot(fit_coef_cohort_sh %>% filter(par_name == 'rate')) + geom_line(aes(x = cohort, y = 1/par_value))
ggsave("Figs/Fig_mean_geom_by4yr.jpg")