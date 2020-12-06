rm (list = ls())
source("0 - Config.R")
clean_data<- read_csv("Output/clean_data.csv")

####################################################
# 2004-2019
####################################################
data_stats <- clean_data %>% group_by (survey_year) %>% do ({
  yr_data <- .
  bs_stats <- data_frame(iter = 1:NUM_RESAMP) %>% group_by(iter) %>% do({
    ind <- sample(x = 1:nrow(yr_data),size = nrow(yr_data),replace = TRUE)
    resamp_data <- yr_data[ind,]
    data_frame (TF_mean = mean(resamp_data$tf_prev), TF_5 = sum(resamp_data$tf_prev < 5)/nrow(yr_data))
  })
  data_frame(TF_5_2p5 = quantile(bs_stats$TF_5,0.025), TF_5_5 = quantile(bs_stats$TF_5,0.05), TF_5_50 = quantile(bs_stats$TF_5,0.5), TF_5_95 = quantile(bs_stats$TF_5,0.95), TF_5_97p5 = quantile(bs_stats$TF_5,0.975),
             TF_mean_2p5 = quantile(bs_stats$TF_mean,0.025), TF_mean_5 = quantile(bs_stats$TF_mean,0.05), TF_mean_50 = quantile(bs_stats$TF_mean,0.5), TF_mean_95 = quantile(bs_stats$TF_mean,0.95), TF_mean_97p5 = quantile(bs_stats$TF_mean,0.975))
})

####################################################
# 2020-2030
####################################################
baseline_stats <- read.csv("from_rigel_51320/baseline_stats.csv")
bs_stats <- read_csv("from_rigel_51320/admin0_stats_1000.csv")

bs_stat_sum <- bs_stats %>% filter(year >= 2020) %>% group_by(year,fit) %>% do({
  data_frame(TF_5_2p5 = quantile(.$TF_5,0.025), TF_5_5 = quantile(.$TF_5,0.05), TF_5_50 = quantile(.$TF_5,0.5), TF_5_95 = quantile(.$TF_5,0.95), TF_5_97p5 = quantile(.$TF_5,0.975),
             TF_mean_2p5 = quantile(.$TF_mean,0.025), TF_mean_5 = quantile(.$TF_mean,0.05), TF_mean_50 = quantile(.$TF_mean,0.5), TF_mean_95 = quantile(.$TF_mean,0.95), TF_mean_97p5 = quantile(.$TF_mean,0.975))
})

####################################################
# Plot TF_mean and TF_5 statistics - exponential
####################################################
baseline_plot_stats <- baseline_stats %>% filter(fit == "renorm_exp", year > 2019)
bs_plot_sum <- bs_stat_sum %>% filter(fit == "renorm_exp")

p1 <- ggplot () + xlab("Year") + ylab("Average TF") +
  geom_errorbar(data = data_stats, aes(x = survey_year, ymin = TF_mean_2p5, ymax = TF_mean_97p5)) +
  geom_point(data = data_stats, aes(x = survey_year, y = TF_mean_50)) +
  geom_errorbar(data = bs_plot_sum, aes(x = year, ymin = TF_mean_2p5, ymax = TF_mean_97p5), col = "blue") +
#  geom_point(data = baseline_plot_stats, aes(x = year, y = TF_mean), col = "grey") +
  geom_point(data = bs_plot_sum, aes(x = year, y = TF_mean_50), shape = 1, col = "blue")
#(p1)

p2 <- ggplot ()  + xlab("Year") + ylab("Proportion TF< 5%") +
  geom_errorbar(data = data_stats, aes(x = survey_year, ymin = TF_5_2p5, ymax = TF_5_97p5)) +
  geom_point(data = data_stats, aes(x = survey_year, y = TF_5_50)) +
  geom_errorbar(data = bs_plot_sum, aes(x = year, ymin = TF_5_2p5, ymax = TF_5_97p5), col = "blue") +
#  geom_point(data = baseline_plot_stats, aes(x = year, y = TF_5), col = "grey") +
  geom_point(data = bs_plot_sum, aes(x = year, y = TF_5_50), shape = 1, col = "blue")
#(p2)

(p3 <- grid.arrange(p1,p2))
ggsave("Figs/Fig_trends_exp.jpg", plot = p3)

####################################################
# Data for paper
####################################################
exp2025 <- bs_stats %>% filter(year == 2025, fit == "renorm_exp")
exp2030 <- bs_stats %>% filter(year == 2030, fit == "renorm_exp")

# ####################################################
# # Plot R estimate
# ####################################################
# bs_stats_2020 <- bs_stats %>% filter(year == 2020, fit == "renorm_exp")
# bs_stats_2021 <- bs_stats %>% filter(year == 2021, fit == "renorm_exp")
# R_2020 <- tibble(R = sqrt(bs_stats_2021$TF_mean / bs_stats_2020$TF_mean))
# ggplot(R_2020) + geom_histogram(aes(x=R),breaks = c(seq(0.9,1.0,.0025)))
# (paste("Mean R: ",mean(R_2020$R)))
# (paste("Std Dev R: ",sd(R_2020$R)))
# 
# ####################################################
# # Plot other trends
# ####################################################
# p1_all <- ggplot () +
#   geom_errorbar(data = data_stats, aes(x = year, ymin = TF_mean_2p5, ymax = TF_mean_97p5)) +
#   geom_point(data = data_stats, aes(x = year, y = TF_mean_50)) +
#   geom_errorbar(data = bs_stat_sum, aes(x = year, ymin = TF_mean_2p5, ymax = TF_mean_97p5), col = "blue") +
#   geom_point(data = bs_stat_sum, aes(x = year, y = TF_mean_50), col = "blue") +
# #  geom_point(data = baseline_stats, aes(x = year, y = TF_mean), col = "grey") +
#   facet_wrap (~fit, nrow = 1)
# (p1_all)
# 
# p2_all <- ggplot () +
#   geom_errorbar(data = data_stats, aes(x = year, ymin = TF_5_2p5, ymax = TF_5_97p5)) +
#   geom_point(data = data_stats, aes(x = year, y = TF_5_50)) +
#   geom_errorbar(data = bs_stat_sum, aes(x = year, ymin = TF_5_2p5, ymax = TF_5_97p5), col = "blue") +
#   geom_point(data = bs_stat_sum, aes(x = year, y = TF_5_50), col = "blue") +
# #  geom_point(data = baseline_plot_stats, aes(x = year, y = TF_5), col = "grey") +
#   facet_wrap (~fit, nrow = 1)
# (p2_all)
# 
# (p3_all <- grid.arrange(p1_all,p2_all))
# #ggsave("Figs/Fig_trends_all.jpg", plot = p3_all)

