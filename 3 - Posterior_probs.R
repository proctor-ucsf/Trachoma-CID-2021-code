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
# Determine best fit pdfs for year cohorts
####################################################
# Divide data into cohorts
# Then fit gamma to full range and exponential to data with TF < EXP_THRESH
clean_data <- clean_data %>% mutate(cohort = floor((year - min(CLEAN_YRS))/COHORT_SIZE))
if (file.exists("Output/fit_coef_post_sh.csv")) {
  fit_coef_post_sh <- read.csv("Output/fit_coef_post_sh.csv")
} else {
  fit_coef_post_sh <- clean_data %>% group_by(cohort) %>% do({
    .$survey_year = .$cohort[1]
    gamma_fit <- mle_fit(data = ., fits = 'renorm_gamma')
    gamma_coef <- det_coef(gamma_fit,yrs = .$year[1])
    subdata <- as_tibble(.) %>% filter(tf_prev < EXP_THRESH)
    exp_fit <- mle_fit(data = subdata, fits = 'renorm_exp_trunc')
    exp_coef <- det_coef(exp_fit,yrs = .$year[1])
    bind_rows(gamma_coef, exp_coef)
  })
  write.csv(fit_coef_post_sh, "Output/fit_coef_post_sh.csv")
}
fit_coef_post_sh$fit[fit_coef_post_sh$fit == 'renorm_exp_trunc'] = 'renorm_exp'
fit_pdfs_post_sh <- pdf_gen_bindata(fit_coef_post_sh)
fit_pdfs_post_sh <- fit_pdfs_post_sh %>% mutate(cohort = floor((year - min(CLEAN_YRS))/COHORT_SIZE))

####################################################
# Plot fits
####################################################
cohort_count<-clean_data %>% group_by(cohort) %>% summarize(count = n())
fit_pdfs_post_sh <- left_join(fit_pdfs_post_sh,cohort_count, by = 'cohort')
ggplot () +
  geom_histogram(data = clean_data, aes(x = tf_prev),binwidth = DISPLAY_WIDTH) +
  #geom_freqpoly(data = clean_data, aes(x = tf_prev),binwidth = 5) +
  geom_line(data = fit_pdfs_post_sh, aes(x = xvar, y = yvar_pdf*count*DISPLAY_WIDTH, col = fit)) +
#  ggtitle(paste("Data fits: ",string, sep = "")) +
  facet_wrap(~ cohort)

####################################################
# Calculate posterior probability
####################################################
prob_wider <- fit_pdfs_post_sh %>% pivot_wider(names_from = fit, values_from = yvar_pdf) %>%
  mutate (posterior = (renorm_gamma-renorm_exp)/renorm_gamma)

####################################################
# Plot posterior
####################################################
cohort_name <- tibble(cohort = 0:3, Years = c('2004-2007','2008-2011','2012-2015','2016-2019'))
prob_wider <- left_join(prob_wider, cohort_name, by = 'cohort')
ggplot(prob_wider %>% filter(posterior> .1) %>% filter(xvar > 5)) +
  xlab("TF prevalence (%)") + ylab("Probability of supercritical transmission") +
  geom_line(aes(x= xvar, y = posterior, col = Years)) + xlim(0,100) + theme(text = element_text(size = 20)) +
  theme(legend.position = c(.75,.25))
ggsave("Figs/Fig_posterior_probability.jpg")

####################################################
# Toy plot
####################################################
x <- 1:100
y1 <- dexp(1:100,1/50)
y2 <- dgamma(1:100,rate = 1/100, shape = 0.5)
aa<- tibble(Prevalence = x, Density = y1, y2 = y2)

ggplot(aa)  + theme(text = element_text(size = 40)) +
  geom_line(aes(x =Prevalence, y = Density), col = 'blue') +
  geom_line(aes(x =Prevalence, y = y2), col = 'black')