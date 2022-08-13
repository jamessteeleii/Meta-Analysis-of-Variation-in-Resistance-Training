### This script runs all the analyses for the paper

# title: "Meta-Analysis of Variation in Sport and Exercise Science"
# subtitle: "Examples of Application Within Resistance Training Research"
# author: 
#   - James Steele
# - James Fisher
# - Dave Smith
# - Andrew Vigotsky
# - Brad Schoenfeld
# - Yefeng Yang
# - Shinichi Nakagawa

# Open packages
library(metafor)
library(tidyverse)
library(patchwork)
library(europepmc)

##### Trend of Meta-Analysis in Resistance Training

SES_meta_trend <- epmc_hits_trend(query = 'TITLE: (sport OR exercise) AND "meta analysis"',
                                  period = 1976:2022) %>%
  mutate(`Search string` = "(sport OR exercise) AND meta analysis")

RT_meta_trend <- epmc_hits_trend(query = 'TITLE: ("resistance training" OR "strength training") AND "meta analysis"',
                                 period = 1976:2022) %>%
  mutate(`Search string` = "(resistance training OR strength training) AND meta analysis")

trends <- rbind(SES_meta_trend, RT_meta_trend) %>%
  mutate(`Search string` = factor(`Search string`, levels = c('(sport OR exercise) AND meta analysis',
                                                              '(resistance training OR strength training) AND meta analysis')))

trends_plot <- trends %>%
  ggplot(aes(x = year, y = (query_hits / all_hits)*100, linetype = `Search string`)) +
  geom_line(size = 1) +
  labs(x = "Year", y = "Proportion of Published Articles (%)",
       title = "Articles in Europe Pub Med Central Database") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

save(trends_plot, file = "plots/trends_plot")

##### Read csv as data frame into environment - Note: change source address
Data <- read.csv(here::here("data","Polito et al. RT Extracted Data.csv"), na.strings=c(""," ","NA"))

# Calculate pre-post SDs from SEs
Data$RT_pre_sd <- ifelse(is.na(Data$RT_pre_se), Data$RT_pre_sd, Data$RT_pre_se * sqrt(Data$RT_n))
Data$CON_pre_sd <- ifelse(is.na(Data$CON_pre_se), Data$CON_pre_sd, Data$CON_pre_se * sqrt(Data$CON_n))
Data$RT_post_sd <- ifelse(is.na(Data$RT_post_se), Data$RT_post_sd, Data$RT_post_se * sqrt(Data$RT_n))
Data$CON_post_sd <- ifelse(is.na(Data$CON_post_se), Data$CON_post_sd, Data$CON_post_se * sqrt(Data$CON_n))

# Convert p to t (Change scores)
Data$RT_delta_t_value <- replmiss(Data$RT_delta_t_value, with(Data, qt(RT_delta_p_value/2, df=RT_n-1, lower.tail=FALSE)))
Data$CON_delta_t_value <- replmiss(Data$CON_delta_t_value, with(Data, qt(CON_delta_p_value/2, df=CON_n-1, lower.tail=FALSE)))

# Convert t to SE (Change scores)
Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, ifelse(is.na(RT_delta_m), 
                                                                 (RT_post_m - RT_pre_m)/RT_delta_t_value, RT_delta_m/RT_delta_t_value)))
Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, ifelse(is.na(CON_delta_m), 
                                                                 (CON_post_m - CON_pre_m)/CON_delta_t_value, CON_delta_m/CON_delta_t_value)))
# Make positive
Data$RT_delta_se <- ifelse(Data$RT_delta_se < 0, Data$RT_delta_se * -1, Data$RT_delta_se)
Data$CON_delta_se <- ifelse(Data$CON_delta_se < 0, Data$CON_delta_se * -1, Data$CON_delta_se)

# Convert CI to SE (Change scores)
Data$RT_delta_se <- replmiss(Data$RT_delta_se, with(Data, (RT_delta_CI_upper - RT_delta_CI_lower)/3.92))
Data$CON_delta_se <- replmiss(Data$CON_delta_se, with(Data, (CON_delta_CI_upper - CON_delta_CI_lower)/3.92))

# Convert SE to SD (Change scores)
Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, RT_delta_se * sqrt(RT_n)))
Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, CON_delta_se * sqrt(CON_n)))

# Calculate pre-post correlation coefficient for those with pre, post, and delta SDs
Data$RT_ri <- (Data$RT_pre_sd^2 + Data$RT_post_sd^2 - Data$RT_delta_sd^2)/(2 * Data$RT_pre_sd * Data$RT_post_sd)
Data$CON_ri <- (Data$CON_pre_sd^2 + Data$CON_post_sd^2 - Data$CON_delta_sd^2)/(2 * Data$CON_pre_sd * Data$CON_post_sd)

# Remove values outside the range of -1 to +1 as they are likely due to misreporting or miscalculations in original studies
Data$RT_ri <- ifelse(between(Data$RT_ri,-1,1) == FALSE, NA, Data$RT_ri)
Data$CON_ri <- ifelse(between(Data$CON_ri,-1,1) == FALSE, NA, Data$CON_ri)

# Then we'll convert using Fishers r to z, calculate a meta-analytic point estimate, and impute that across the studies with missing correlations
Data <- escalc(measure = "ZCOR", ri = RT_ri, ni = RT_n, data = Data)

Meta_RT_ri <- rma.mv(yi, V=vi, data=Data,
                     slab=paste(label),
                     random = list(~ 1 | study, ~ 1 | group), method="REML",
                     control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEstMeta_RT_ri <- robust(Meta_RT_ri, Data$study)

z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b[1])

Data$RT_ri <- ifelse(is.na(Data$RT_ri), z2r_RT, Data$RT_ri)

Data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = Data)

### Note, data is coded with study and group as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | group)
Meta_CON_ri <- rma.mv(yi, V=vi, data=Data,
                     slab=paste(label),
                     random = list(~ 1 | study, ~ 1 | group), method="REML",
                     control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEstMeta_CON_ri <- robust(Meta_CON_ri, Data$study)

z2r_CON <- psych::fisherz2r(RobuEstMeta_CON_ri$b[1])

Data$CON_ri <- ifelse(is.na(Data$CON_ri), z2r_CON, Data$CON_ri)

# Estimate change score difference SD where only pre-post data available
Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, sqrt(RT_pre_sd^2 + RT_post_sd^2 - (2*RT_ri*RT_pre_sd*RT_post_sd))))
Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, sqrt(CON_pre_sd^2 + CON_post_sd^2 - (2*CON_ri*CON_pre_sd*CON_post_sd))))

###### Section - Detecting the presence of interindividual response variation to resistance training intervention

###### Comparison of standardised mean changes between RT and CON ######

### Standardised mean difference effect size calculations 
Data$SD_pool <- sqrt(((Data$RT_n - 1)*Data$RT_pre_sd^2 + (Data$CON_n - 1)*Data$CON_pre_sd^2) / (Data$RT_n + Data$CON_n - 2))

Data_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                   m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=RT_ri, data = Data)
Data_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                          m2i=CON_pre_m, sd1i=SD_pool, ni=CON_n, ri=CON_ri, data = Data)
Data_SMD <- Data

Data_SMD$yi <- (Data_RT$yi - Data_CON$yi)
Data_SMD$vi <- (Data_RT$vi + Data_CON$vi)

### Strength
Data_SMD_strength <- Data_SMD %>% 
    filter(!is.na(yi) &  outcome == "strength")

MultiLevelModel_strength <- rma.mv(yi, V=vi, data=Data_SMD_strength,
                                         slab=paste(label),
                                         random = list(~ 1 | study, ~ 1 | group), method="REML",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength, file = "models/MultiLevelModel_strength")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength$vi)
X <- model.matrix(MultiLevelModel_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength <- 100 * sum(MultiLevelModel_strength$sigma2) / (sum(MultiLevelModel_strength$sigma2) + (MultiLevelModel_strength$k-MultiLevelModel_strength$p)/sum(diag(P)))
I2bw_strength <- 100 * MultiLevelModel_strength$sigma2 / (sum(MultiLevelModel_strength$sigma2) + (MultiLevelModel_strength$k-MultiLevelModel_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength <- robust(MultiLevelModel_strength, Data_SMD_strength$study)

save(RobuEstMultiLevelModel_strength, file = "models/RobuEstMultiLevelModel_strength")

### Caterpillar plot 

# Overall estimate
diamond_strength <- data.frame(x = c(RobuEstMultiLevelModel_strength$b[1] + (RobuEstMultiLevelModel_strength$se*1.96),
                            RobuEstMultiLevelModel_strength$b[1],
                            RobuEstMultiLevelModel_strength$b[1] - (RobuEstMultiLevelModel_strength$se*1.96),
                            RobuEstMultiLevelModel_strength$b[1]),
                      y = c(-15,-25,-15,-5))

# Prediction interval
PI_strength <- as.data.frame(predict(RobuEstMultiLevelModel_strength))

# I2 labels
I2_strength_lab <- data.frame(level = c("study", "group"),
                              I2 = I2bw_strength) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_strength <- Data_SMD_strength %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=es)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-2,7.5), breaks = c(-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  scale_y_discrete(limits=rev, expand = expansion(mult = c(0.075,0))) +
  geom_text(data = mutate_if(PI_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 5, y = 90), hjust = "centre", size = 3) +
  geom_text(data = I2_strength_lab,
            aes(label = glue::glue("I Squared [Study = {round(study,1)}%; Group = {round(group,1)}%]"),
                x = 5, y = 70), hjust = "centre", size = 3) +
  geom_segment(data = PI_strength, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_strength, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standardised Mean Difference (Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

### Hypertrophy
Data_SMD_hypertrophy <- Data_SMD %>% 
    filter(!is.na(yi) &  outcome == "hypertrophy")

MultiLevelModel_hypertrophy <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy,
                                   slab=paste(label),
                                   random = list(~ 1 | study, ~ 1 | group), method="REML",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy, file = "models/MultiLevelModel_hypertrophy")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy$vi)
X <- model.matrix(MultiLevelModel_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy <- 100 * sum(MultiLevelModel_hypertrophy$sigma2) / (sum(MultiLevelModel_hypertrophy$sigma2) + (MultiLevelModel_hypertrophy$k-MultiLevelModel_hypertrophy$p)/sum(diag(P)))
I2bw_hypertrophy <- 100 * MultiLevelModel_hypertrophy$sigma2 / (sum(MultiLevelModel_hypertrophy$sigma2) + (MultiLevelModel_hypertrophy$k-MultiLevelModel_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy <- robust(MultiLevelModel_hypertrophy, Data_SMD_hypertrophy$study)

save(RobuEstMultiLevelModel_hypertrophy, file = "models/RobuEstMultiLevelModel_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_hypertrophy$b[1] + (RobuEstMultiLevelModel_hypertrophy$se*1.96),
                            RobuEstMultiLevelModel_hypertrophy$b[1],
                            RobuEstMultiLevelModel_hypertrophy$b[1] - (RobuEstMultiLevelModel_hypertrophy$se*1.96),
                            RobuEstMultiLevelModel_hypertrophy$b[1]),
                      y = c(-15,-25,-15,-5))

# Prediction interval
PI_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_hypertrophy))

# I2 labels
I2_hypertrophy_lab <- data.frame(level = c("study", "group"),
                              I2 = I2bw_hypertrophy) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_hypertrophy <- Data_SMD_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=es)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-2,7.5), breaks = c(-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  scale_y_discrete(limits=rev, expand = expansion(mult = c(0.075,0))) +
  geom_text(data = mutate_if(PI_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 5, y = 90), hjust = "centre", size = 3) +
  geom_text(data = I2_hypertrophy_lab,
            aes(label = glue::glue("I Squared [Study = {round(study,1)}%; Group = {round(group,1)}%]"),
                x = 5, y = 70), hjust = "centre", size = 3) +
  geom_segment(data = PI_hypertrophy, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_hypertrophy, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standardised Mean Difference (Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())


### Combine Plots
forest_plots <- (forest_strength / forest_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_plots, file = "plots/forest_plots")

forest_plots

ggsave("plots/SMD_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)


###### Comparison of variance of change between RT and CON ######

### Using the SDir - standard deviation for individual response
Data_SDir <- Data %>%
  filter(!is.na(RT_delta_sd) | !is.na(CON_delta_sd) | CON_delta_sd == 0) %>%
  mutate(SDir = sqrt(pmax(0, RT_delta_sd^2 - CON_delta_sd^2)),
         SDir_se = sqrt(2*((RT_delta_sd^4)/(RT_n-1)+(CON_delta_sd^4)/(CON_n-1)))) %>%
  filter(!is.na(SDir))

### Strength
Data_SDir_strength <- Data_SDir %>% 
  filter(outcome == "strength")

MultiLevelModel_SDir_strength <- rma.mv(SDir, V=SDir_se^2, data=Data_SDir_strength,
                                         slab=paste(label),
                                         random = list(~ 1 | study, ~ 1 | group), method="REML",
                                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SDir_strength, file = "models/MultiLevelModel_SDir_strength")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SDir_strength$SDir_se^2)
X <- model.matrix(MultiLevelModel_SDir_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_SDir_strength <- 100 * sum(MultiLevelModel_SDir_strength$sigma2) / (sum(MultiLevelModel_SDir_strength$sigma2) + (MultiLevelModel_SDir_strength$k-MultiLevelModel_SDir_strength$p)/sum(diag(P)))
I2bw_SDir_strength <- 100 * MultiLevelModel_SDir_strength$sigma2 / (sum(MultiLevelModel_SDir_strength$sigma2) + (MultiLevelModel_SDir_strength$k-MultiLevelModel_SDir_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SDir_strength <- robust(MultiLevelModel_SDir_strength, Data_SDir_strength$study)

save(RobuEstMultiLevelModel_SDir_strength, file = "models/RobuEstMultiLevelModel_SDir_strength")

### Hypertrophy
Data_SDir_hypertrophy <- Data_SDir %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_SDir_hypertrophy <- rma.mv(SDir, V=SDir_se^2, data=Data_SDir_hypertrophy,
                                            slab=paste(label),
                                            random = list(~ 1 | study, ~ 1 | group), method="REML",
                                            control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SDir_hypertrophy, file = "models/MultiLevelModel_SDir_hypertrophy")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SDir_hypertrophy$SDir_se^2)
X <- model.matrix(MultiLevelModel_SDir_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_SDir_hypertrophy <- 100 * sum(MultiLevelModel_SDir_hypertrophy$sigma2) / (sum(MultiLevelModel_SDir_hypertrophy$sigma2) + (MultiLevelModel_SDir_hypertrophy$k-MultiLevelModel_SDir_hypertrophy$p)/sum(diag(P)))
I2bw_SDir_hypertrophy <- 100 * MultiLevelModel_SDir_hypertrophy$sigma2 / (sum(MultiLevelModel_SDir_hypertrophy$sigma2) + (MultiLevelModel_SDir_hypertrophy$k-MultiLevelModel_SDir_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SDir_hypertrophy <- robust(MultiLevelModel_SDir_hypertrophy, Data_SDir_hypertrophy$study)

save(RobuEstMultiLevelModel_SDir_hypertrophy, file = "models/RobuEstMultiLevelModel_SDir_hypertrophy")

### Using the Log Variability Ratio
Data_logVR <- escalc(measure = "VR", sd1i = RT_delta_sd, sd2i = CON_delta_sd, n1i = RT_n, n2i = CON_n, data = Data)
Data_logVR <- Data_logVR %>% 
    filter(!is.na(yi))

### Strength
Data_logVR_strength <- Data_logVR %>% 
  filter(outcome == "strength")

MultiLevelModel_logVR_strength <- rma.mv(yi, V=vi, data=Data_logVR_strength,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), method="REML",
                                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logVR_strength, file = "models/MultiLevelModel_logVR_strength")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_logVR_strength$vi)
X <- model.matrix(MultiLevelModel_logVR_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logVR_strength <- 100 * sum(MultiLevelModel_logVR_strength$sigma2) / (sum(MultiLevelModel_logVR_strength$sigma2) + (MultiLevelModel_logVR_strength$k-MultiLevelModel_logVR_strength$p)/sum(diag(P)))
I2bw_logVR_strength <- 100 * MultiLevelModel_logVR_strength$sigma2 / (sum(MultiLevelModel_logVR_strength$sigma2) + (MultiLevelModel_logVR_strength$k-MultiLevelModel_logVR_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_strength <- robust(MultiLevelModel_logVR_strength, Data_logVR_strength$study)

save(RobuEstMultiLevelModel_logVR_strength, file = "models/RobuEstMultiLevelModel_logVR_strength")

### Hypertrophy
Data_logVR_hypertrophy <- Data_logVR %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_logVR_hypertrophy <- rma.mv(yi, V=vi, data=Data_logVR_hypertrophy,
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), method="REML",
                                            control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logVR_hypertrophy, file = "models/MultiLevelModel_logVR_hypertrophy")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_logVR_hypertrophy$vi)
X <- model.matrix(MultiLevelModel_logVR_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logVR_hypertrophy <- 100 * sum(MultiLevelModel_logVR_hypertrophy$sigma2) / (sum(MultiLevelModel_logVR_hypertrophy$sigma2) + (MultiLevelModel_logVR_hypertrophy$k-MultiLevelModel_logVR_hypertrophy$p)/sum(diag(P)))
I2bw_logVR_hypertrophy <- 100 * MultiLevelModel_logVR_hypertrophy$sigma2 / (sum(MultiLevelModel_logVR_hypertrophy$sigma2) + (MultiLevelModel_logVR_hypertrophy$k-MultiLevelModel_logVR_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_hypertrophy <- robust(MultiLevelModel_logVR_hypertrophy, Data_logVR_hypertrophy$study)

save(RobuEstMultiLevelModel_logVR_hypertrophy, file = "models/RobuEstMultiLevelModel_logVR_hypertrophy")

### Combine SMD, SDir, and logVR estimates

SMD_SDir_logVR <- rbind(data.frame(Outcome = "strength",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_strength$ci.ub, 2),
                                    `I2 study` = round(I2bw_strength, 2)[1],
                                    `I2 group` = round(I2bw_strength, 2)[2]),
                        data.frame(Outcome = "strength",
                                   Model = "SDir",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_SDir_strength$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_SDir_strength$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_SDir_strength$ci.ub, 2),
                                   `I2 study` = round(I2bw_SDir_strength, 2)[1],
                                   `I2 group` = round(I2bw_SDir_strength, 2)[2]),
                        data.frame(Outcome = "strength",
                                   Model = "logVR",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_logVR_strength$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_logVR_strength$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_logVR_strength$ci.ub, 2),
                                   `I2 study` = round(I2bw_logVR_strength, 2)[1],
                                   `I2 group` = round(I2bw_logVR_strength, 2)[2]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "SMD",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_hypertrophy$ci.ub, 2),
                                   `I2 study` = round(I2bw_hypertrophy, 2)[1],
                                   `I2 group` = round(I2bw_hypertrophy, 2)[2]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "SDir",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_SDir_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_SDir_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_SDir_hypertrophy$ci.ub, 2),
                                   `I2 study` = round(I2bw_SDir_hypertrophy, 2)[1],
                                   `I2 group` = round(I2bw_SDir_hypertrophy, 2)[2]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "logVR",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_logVR_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_logVR_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_logVR_hypertrophy$ci.ub, 2),
                                   `I2 study` = round(I2bw_logVR_hypertrophy, 2)[1],
                                   `I2 group` = round(I2bw_logVR_hypertrophy, 2)[2])
)

save(SMD_SDir_logVR, file = "models/SMD_SDir_logVR")

###### Section - Mean-variance relationships in muscular strength, endurance, and hypertrophy
### We'll look at just the baseline strength data from the Polito et al. data
Data_long_pre_m <- Data %>%
  select(study, group, es, RT_n, CON_n, RT_pre_m, CON_pre_m, outcome) %>%
  pivot_longer(c(RT_pre_m, CON_pre_m),
                   names_to = "key",
                   values_to = "mean")

Data_long_pre_m$key <- recode(Data_long_pre_m$key, RT_pre_m = "RT", CON_pre_m = "CON")

Data_long_pre_sd <- Data %>%
  select(study, group, es, RT_n, CON_n, RT_pre_sd, CON_pre_sd, outcome) %>%
  pivot_longer(c(RT_pre_sd, CON_pre_sd),
                   names_to = "key",
                   values_to = "sd")

Data_long_pre_sd$key <- recode(Data_long_pre_sd$key, RT_pre_sd = "RT", CON_pre_sd = "CON")

Data_long_pre <- cbind(Data_long_pre_m, sd = Data_long_pre_sd$sd) 

Data_long_pre_RT <- Data_long_pre %>%
  filter(key == "RT") %>%
  mutate(n = RT_n,
         group = as.factor(unclass(factor(unlist(group))))) %>%
  select(study, group, es, outcome, key, mean, sd, n)

Data_long_pre_CON <- Data_long_pre %>%
  filter(key == "CON") %>%
  distinct() %>%
  mutate(n = CON_n,
         group = as.factor(unclass(factor(unlist(group)))+length(unique(Data_long_pre_RT$group)))) %>%
  select(study, group, es, outcome, key, mean, sd, n)

Data_long_pre <- rbind(Data_long_pre_RT, Data_long_pre_CON)

# Calculate log SD and variance of log SD
Data_long_pre$SD_log <- log(Data_long_pre$sd) + (1/(2*(Data_long_pre$n-1)))
Data_long_pre$SD_log_vi <- (1/(2*(Data_long_pre$n-1)))


# Plot raw mean and SD
m_sd_strength <- ggplot(subset(Data_long_pre, outcome == "strength"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Strength Outcomes")

m_sd_hypertrophy <- ggplot(subset(Data_long_pre, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd), alpha = 0.2) +
  labs(x = "Mean", y = "Standard Deviation") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Hypertrophy Outcomes")

m_sd_strength_log <- ggplot(subset(Data_long_pre, outcome == "strength"), aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic() 

m_sd_hypertrophy_log <- ggplot(subset(Data_long_pre, outcome == "hypertrophy"), aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log), alpha = 0.2) +
  labs(x = "Log Mean", y = "Log Standard Deviation") +
  guides(fill = "none") +
  theme_classic() 

# Plot together
mean_variance_pre_plots <- ((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 

save(mean_variance_pre_plots, file = "plots/mean_variance_pre_plots")

mean_variance_pre_plots

ggsave("plots/mean_variance_pre_plots.tiff", width = 10, height = 7.5, device = "tiff", dpi = 300)


### Fitting a random slope mixed effects model for mean-variance
# We will include the outcome type as a moderator with random slopes 
# Then the model is sort of similar to the one we'll use for change scores below to introduce it

MultiLevelModel_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                random = list(~ outcome | study, ~ outcome | group), struct = "HCS",
                                                mods = ~ log(mean) + outcome,
                                                method="REML",
                                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
                                            )

save(MultiLevelModel_log_mean_variance, file = "models/MultiLevelModel_log_mean_variance")


### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_log_mean_variance <- robust(MultiLevelModel_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_log_mean_variance, file = "models/RobuEstMultiLevelModel_log_mean_variance")

### Meta-analytic scatter plot

# get the predicted log values
Data_long_pre <- Data_long_pre %>%
  filter(!is.na(SD_log))

Data_long_pre <- cbind(Data_long_pre, pred = predict(RobuEstMultiLevelModel_log_mean_variance)$pred,
                      ci.lb =  predict(RobuEstMultiLevelModel_log_mean_variance)$ci.lb,
                      ci.ub =  predict(RobuEstMultiLevelModel_log_mean_variance)$ci.ub) %>%
  mutate(wi = 1/sqrt(SD_log_vi),
         size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

model_mean_variance_pre_plots <- ggplot(Data_long_pre, aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log, color = outcome, size = size), alpha = 0.05) +
  geom_ribbon(aes(x=log(mean), ymax=ci.ub, ymin=ci.lb, fill = outcome), alpha = 0.2) +
  geom_line(aes(x=log(mean), y=pred, color = outcome)) +
  scale_fill_manual("Outcome", values = alpha(c("#009E73", "#D55E00"),0.5)) +
  scale_color_manual("Outcome", values = alpha(c("#009E73", "#D55E00"),0.5)) +
  labs(x = "Log Mean", y = "Log Standard Deviation", color = "Outcome", shape = "", fill = "") +
  theme_classic() +
  guides(size = "none", fill = "none")

save(model_mean_variance_pre_plots, file = "plots/model_mean_variance_pre_plots")

model_mean_variance_pre_plots

ggsave("plots/model_mean_variance_pre_plots.tiff", width = 7.5, height = 5, device = "tiff", dpi = 300)

###### Section - Reanalysis of interindividual response variation using $lnCVR$ and the random slope mixed effects model

### Change score mean-variances

# Let's take a look at the raw SD and mean, and also the log SD and log mean relatioships

# Add missing deltas
Data$RT_delta_m <- replmiss(Data$RT_delta_m, with(Data, RT_post_m - RT_pre_m))
Data$CON_delta_m <- replmiss(Data$CON_delta_m, with(Data, CON_post_m - CON_pre_m))

# First, we check to see what the relationship between mean and variance for change scores looks like
Data_long_m <- Data %>%
  select(study, group, es, RT_n, CON_n, RT_delta_m, CON_delta_m, outcome) %>%
  pivot_longer(c(RT_delta_m, CON_delta_m),
               names_to = "key",
               values_to = "mean")

Data_long_m$key <- recode(Data_long_m$key, RT_delta_m = "RT", CON_delta_m = "CON")

Data_long_sd <- Data %>%
  select(study, group, es, RT_n, CON_n, RT_delta_sd, CON_delta_sd, outcome) %>%
  pivot_longer(c(RT_delta_sd, CON_delta_sd),
               names_to = "key",
               values_to = "sd")

Data_long_sd$key <- recode(Data_long_sd$key, RT_delta_sd = "RT", CON_delta_sd = "CON")

Data_long <- cbind(Data_long_m, sd = Data_long_sd$sd) 

Data_long_RT <- Data_long %>%
  filter(key == "RT") %>%
  mutate(n = RT_n,
         group = as.factor(unclass(factor(unlist(group))))) %>%
  select(study, group, es, outcome, key, mean, sd, n)

Data_long_CON <- Data_long %>%
  filter(key == "CON") %>%
  distinct() %>%
  mutate(n = CON_n,
         group = as.factor(unclass(factor(unlist(group)))+length(unique(Data_long_RT$group)))) %>%
  select(study, group, es, outcome, key, mean, sd, n)

Data_long <- rbind(Data_long_RT, Data_long_CON)

# Calculate log SD and variance of log SD
Data_long$SD_log <- log(Data_long$sd) + (1/(2*(Data_long$n-1)))
Data_long$SD_log_vi <- (1/(2*(Data_long$n-1)))

# Plot raw mean and SD
m_sd_strength <- ggplot(subset(Data_long, outcome == "strength"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd, color = key), alpha = 0.2) +
  scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_y_continuous(limits = c(0,1000)) +
  labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Condition") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Strength Outcomes")

m_sd_hypertrophy <- ggplot(subset(Data_long, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd, color = key), alpha = 0.2) +
  scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Condition") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Hypertrophy Outcomes")

# Use the absolute change scores as the relationship between mean and variance seems to roughly hold in both directions about zero
# And because the logarithms are undefined for <=0
Data_long$mean <- ifelse(Data_long$mean < 0, Data_long$mean * -1, Data_long$mean)

# Now let's examine raw unweighted log mean - log variance relationships

m_sd_strength_log <- ggplot(subset(Data_long, outcome == "strength"), aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log, color = key), alpha = 0.2) +
  scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition") +
  guides(fill = "none") +
  theme_classic() 

m_sd_hypertrophy_log <- ggplot(subset(Data_long, outcome == "hypertrophy"), aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log, color = key), alpha = 0.2) +
  scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition") +
  guides(fill = "none") +
  theme_classic() 

# Plot together
mean_variance_delta_plots <- ((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 

save(mean_variance_delta_plots, file = "plots/mean_variance_delta_plots")

mean_variance_delta_plots

ggsave("plots/Mean_variance_delta_plots.tiff", width = 10, height = 7.5, device = "tiff", dpi = 300)

### The assumption of linearity seems to be far better supported for log SD and log mean
### There is far more heteroskedasticity in the raw SD and mean data

### Using the Log Coefficient of Variation Ratios

# First make all the delta scores positive - log CVR only works for ratio data so we need to consider absolute change
Data_logCVR <- Data %>%
  mutate(RT_delta_m = if_else(RT_delta_m < 0, RT_delta_m *-1, RT_delta_m),
         CON_delta_m = if_else(CON_delta_m < 0, CON_delta_m *-1, CON_delta_m))

# Log Coefficient of Variation
Data_logCVR <- escalc(measure = "CVR", m1i = RT_delta_m, m2i = CON_delta_m,
                      sd1i = RT_delta_sd, sd2i = CON_delta_sd, n1i = RT_n, n2i = CON_n, data = Data_logCVR)
Data_logCVR <- Data_logCVR %>% 
    filter(!is.na(yi))

### Strength
Data_logCVR_strength <- Data_logCVR %>% 
  filter(outcome == "strength")

MultiLevelModel_logCVR_strength <- rma.mv(yi, V=vi, data=Data_logCVR_strength,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength, file = "models/MultiLevelModel_logCVR_strength")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_logCVR_strength$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength <- 100 * sum(MultiLevelModel_logCVR_strength$sigma2) / (sum(MultiLevelModel_logCVR_strength$sigma2) + (MultiLevelModel_logCVR_strength$k-MultiLevelModel_logCVR_strength$p)/sum(diag(P)))
I2bw_logCVR_strength <- 100 * MultiLevelModel_logCVR_strength$sigma2 / (sum(MultiLevelModel_logCVR_strength$sigma2) + (MultiLevelModel_logCVR_strength$k-MultiLevelModel_logCVR_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength <- robust(MultiLevelModel_logCVR_strength, Data_logCVR_strength$study)

save(RobuEstMultiLevelModel_logCVR_strength, file = "models/RobuEstMultiLevelModel_logCVR_strength")

### Hypertrophy
Data_logCVR_hypertrophy <- Data_logCVR %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_logCVR_hypertrophy <- rma.mv(yi, V=vi, data=Data_logCVR_hypertrophy,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy, file = "models/MultiLevelModel_logCVR_hypertrophy")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_logCVR_hypertrophy$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy <- 100 * sum(MultiLevelModel_logCVR_hypertrophy$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy$sigma2) + (MultiLevelModel_logCVR_hypertrophy$k-MultiLevelModel_logCVR_hypertrophy$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy <- 100 * MultiLevelModel_logCVR_hypertrophy$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy$sigma2) + (MultiLevelModel_logCVR_hypertrophy$k-MultiLevelModel_logCVR_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy <- robust(MultiLevelModel_logCVR_hypertrophy, Data_logCVR_hypertrophy$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy")


### Let's fit some additional models though addressing some of the possible limitations of log CVR
# See "Limitations of lnCVR and an alternative approach" in https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12309 
# The log CVR assumes SD is proportional to the mean with homoskedasticity

### Using the random slope meta regression model of log SD regressed on log mean with RT/CON as a categorical factor

### Strength
Data_long_strength <- Data_long %>%
    filter(!is.na(sd) &
               !is.na(mean)) %>%
    filter(!is.infinite(SD_log) &
               !is.infinite(SD_log_vi)) %>%
    mutate(mean_log = log(mean)) %>%
    filter(!is.na(mean_log) &
               !is.infinite(mean_log)) %>%
    filter(outcome == "strength")

MultiLevelModel_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                          random = list(~ key | study, ~ key | group), struct = "HCS",
                                          mods = ~ mean_log + key,
                                          method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_log_mean_mod_strength, file = "models/MultiLevelModel_log_mean_mod_strength")


### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_log_mean_mod_strength <- robust(MultiLevelModel_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_log_mean_mod_strength")

### Meta-analytic scatter plot

# get the predicted log values
Data_long_strength_log <- cbind(Data_long_strength, predict(RobuEstMultiLevelModel_log_mean_mod_strength)) %>%
    mutate(wi = 1/sqrt(SD_log_vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

model_m_sd_strength_log <- ggplot(Data_long_strength_log, aes(x=mean_log, y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean_log, y=SD_log, color = key, size = size), alpha = 0.2) +
    geom_ribbon(aes(x=mean_log, ymax=ci.ub, ymin=ci.lb, fill = key), alpha = 0.2) +
    geom_line(aes(x=mean_log, y=pred, color = key)) +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition", shape = "", fill = "") +
    theme_classic() +
    ggtitle("Strength Outcomes") +
    guides(size = "none", fill = "none")

### Hypertrophy
Data_long_hypertrophy <- Data_long %>%
    filter(!is.na(sd) &
               !is.na(mean)) %>%
    filter(!is.infinite(SD_log) &
               !is.infinite(SD_log_vi)) %>%
    mutate(mean_log = log(mean)) %>%
    filter(!is.na(mean_log) &
               !is.infinite(mean_log)) %>%
    filter(outcome == "hypertrophy")

MultiLevelModel_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                             random = list(~ key | study, ~ key | group), struct = "HCS",
                                             mods = ~ mean_log + key,
                                          method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_log_mean_mod_hypertrophy")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_log_mean_mod_hypertrophy <- robust(MultiLevelModel_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_log_mean_mod_hypertrophy")

### Meta-analytic scatter plot

# get the predicted log values
Data_long_hypertrophy_log <- cbind(Data_long_hypertrophy, predict(RobuEstMultiLevelModel_log_mean_mod_hypertrophy)) %>%
    mutate(wi = 1/sqrt(SD_log_vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

model_m_sd_hypertrophy_log <- ggplot(Data_long_hypertrophy_log, aes(x=mean_log, y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean_log, y=SD_log, color = key, size = size), alpha = 0.2) +
    geom_ribbon(aes(x=mean_log, ymax=ci.ub, ymin=ci.lb, fill = key), alpha = 0.2) +
    geom_line(aes(x=mean_log, y=pred, color = key)) +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition", shape = "", fill = "") +
    ggtitle("Hypertrophy Outcomes") +
    theme_classic() +
    guides(size = "none", fill = "none")

model_mean_variance_delta_plots <- (model_m_sd_strength_log | model_m_sd_hypertrophy_log) + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 

save(model_mean_variance_delta_plots, file = "plots/model_mean_variance_delta_plots")

model_mean_variance_delta_plots

ggsave("plots/model_mean_variance_delta_plots.tiff", width = 10, height = 5, device = "tiff", dpi = 300)


###### Stuff for supplementary materials?

####### Summary tables for study/participant characteristics 

### Sample sizes
sample_sizes <- Data %>%
  select(group, RT_n, CON_n) %>%
  group_by(group) %>%
  summarise(RT_n = max(RT_n), 
            CON_n = max(CON_n)) %>%
  summarise(`All` = sum(RT_n),
            `Minumum RT` = min(RT_n),
            `Median RT` = median(RT_n),
            `Maximum RT` = max(RT_n),
            `All CON` = sum(CON_n, na.rm =TRUE), 
            `Minumum CON` = min(CON_n, na.rm =TRUE),
            `Median CON` = median(CON_n, na.rm =TRUE),
            `Maximum CON` = max(CON_n, na.rm =TRUE)) %>%
  pivot_longer(1:8, names_to = "Group", values_to = "Sample Size")

save(sample_sizes, file = "models/sample_sizes")

### Descriptive Tables
Data_characteristics <- Data %>% 
  distinct(group, .keep_all = TRUE) %>%
  select(TESTEX,
         age,
         sex_._male,
         weight,
         bmi,
         train_status,
         healthy_clinical,
         RT_only,
         weeks,
         freq,
         exercises,
         sets_exercise,
         reps,
         load,
         task_failure_y_n
  )

summary <- Data_characteristics %>% gtsummary::tbl_summary(type = c(TESTEX, freq) ~ "continuous")

summary_table <- as.data.frame(summary$table_body) %>%
  mutate(Characteristic = label,
         Summary = stat_0) %>%
  select(Characteristic, Summary) %>%
  filter(Characteristic != "Unknown")

summary_table$Characteristic <- recode(summary_table$Characteristic,
                                       age = "Age",
                                       sex_._male = "Proportion Male",
                                       weight = "Weight",
                                       bmi = "BMI",
                                       train_status = "Training Status",
                                       trained = "Trained",
                                       untrained = "Untrained",
                                       healthy_clinical = "Sample Type",
                                       healthy = "Healthy",
                                       clinical = "Clinical",
                                       RT_only = "RT + Adjuvant Intervention?",
                                       n = "N", y = "Y",
                                       weeks = "Duration (weeks)",
                                       freq = "Weekly Frequency",
                                       exercises = "Number of Exercises",
                                       sets_exercise = "Sets per Exercise",
                                       reps = "Number of Repetitions",
                                       load = "Load (%1RM)",
                                       task_failure_y_n = "Task Failure?"
                                       )
summary_table <- summary_table[-c(6,9,12,21),]

rownames(summary_table) <- NULL


save(summary_table, file = "models/summary_table")


###### Small study bias - for supplementary

# Most studies in our field care about detecting treatment effects, not variance
# So we'll explore small study bias first for these across all effects
# We'll create contour enhanced funnel plots for examination of publication bias 

# run that model
Data_SMD <- Data_SMD %>%
    filter(vi > 0) # filter NAs

MultiLevelModel_all <- rma.mv(yi, V=vi, data=Data_SMD,
                                            slab=paste(label),
                                            random = list(~ 1 | study, ~ 1 | group), method="REML")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_all <- robust(MultiLevelModel_all, Data_SMD$study)

### Contour enhanced funnel plot for examination of publication bias
tiff(here("plots","funnel_plot.tiff"), units="in", width=10, height=7.5, res=300)

### produce funnel plot
my_colors <- c("#009E73", "#D55E00")[(factor(Data$outcome))]
funnelplotdata <- funnel(RobuEstMultiLevelModel_all, yaxis = "sei", 
                         level=c(90, 95, 99), 
                         shade=c("white", "lightgray", "darkgray"), col = alpha(0.01), back = "white",
                         refline=0, legend=TRUE,
                         xlab = "Standardised Mean Difference (Positive Values Favour Resistance Training)") # Contour-enhanced funnel plot

# Plot over points
with(funnelplotdata, points(x, y, col = alpha(my_colors,0.5), pch = 19))

# Add a legend
legend(4.5, 0.3, legend=c("Hypertrophy", "Strength"),
       col=c("#009E73", "#D55E00"), pch = 19, cex=1)

dev.off()


###### Exploring moderators of treatment effects and variances ###### 

### We'll just use logCVR for simplicity here as results were similar to the random slope meta-regression anyway
# This also means that we can calculate I2 for each more simply to compare to the main models

# We don't have complete data on moderators for all studies, so we'll look at each separately

### TESTEX score
# SMD
# Strength
Data_SMD_strength_TESTEX <- Data_SMD_strength %>%
  filter(!is.na(TESTEX))

MultiLevelModel_strength_TESTEX <- rma.mv(yi, V=vi, data=Data_SMD_strength_TESTEX,
                                   slab=paste(label),
                                   random = list(~ 1 | study, ~ 1 | group), mods = ~ TESTEX, method="REML",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_TESTEX, file = "models/MultiLevelModel_strength_TESTEX")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_TESTEX$vi)
X <- model.matrix(MultiLevelModel_strength_TESTEX)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_TESTEX <- 100 * sum(MultiLevelModel_strength_TESTEX$sigma2) / (sum(MultiLevelModel_strength_TESTEX$sigma2) + (MultiLevelModel_strength_TESTEX$k-MultiLevelModel_strength_TESTEX$p)/sum(diag(P)))
I2bw_strength_TESTEX <- 100 * MultiLevelModel_strength_TESTEX$sigma2 / (sum(MultiLevelModel_strength_TESTEX$sigma2) + (MultiLevelModel_strength_TESTEX$k-MultiLevelModel_strength_TESTEX$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_TESTEX <- robust(MultiLevelModel_strength_TESTEX, Data_SMD_strength_TESTEX$study)

save(RobuEstMultiLevelModel_strength_TESTEX, file = "models/RobuEstMultiLevelModel_strength_TESTEX")

# Hypertrophy
Data_SMD_hypertrophy_TESTEX <- Data_SMD_hypertrophy %>%
  filter(!is.na(TESTEX))

MultiLevelModel_hypertrophy_TESTEX <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_TESTEX,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ TESTEX, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_TESTEX, file = "models/MultiLevelModel_hypertrophy_TESTEX")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_TESTEX$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_TESTEX)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_TESTEX <- 100 * sum(MultiLevelModel_hypertrophy_TESTEX$sigma2) / (sum(MultiLevelModel_hypertrophy_TESTEX$sigma2) + (MultiLevelModel_hypertrophy_TESTEX$k-MultiLevelModel_hypertrophy_TESTEX$p)/sum(diag(P)))
I2bw_hypertrophy_TESTEX <- 100 * MultiLevelModel_hypertrophy_TESTEX$sigma2 / (sum(MultiLevelModel_hypertrophy_TESTEX$sigma2) + (MultiLevelModel_hypertrophy_TESTEX$k-MultiLevelModel_hypertrophy_TESTEX$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_TESTEX <- robust(MultiLevelModel_hypertrophy_TESTEX, Data_SMD_hypertrophy_TESTEX$study)

save(RobuEstMultiLevelModel_hypertrophy_TESTEX, file = "models/RobuEstMultiLevelModel_hypertrophy_TESTEX")

# logCVR
Data_logCVR_TESTEX <- Data_logCVR %>%
  filter(!is.na(TESTEX))

# Strength
MultiLevelModel_logCVR_strength_TESTEX <- rma.mv(yi, V=vi, data=subset(Data_logCVR_TESTEX, outcome == "strength"),
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ TESTEX, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_TESTEX, file = "models/MultiLevelModel_logCVR_strength_TESTEX")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_TESTEX, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_TESTEX)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_TESTEX <- 100 * sum(MultiLevelModel_logCVR_strength_TESTEX$sigma2) / (sum(MultiLevelModel_logCVR_strength_TESTEX$sigma2) + (MultiLevelModel_logCVR_strength_TESTEX$k-MultiLevelModel_logCVR_strength_TESTEX$p)/sum(diag(P)))
I2bw_logCVR_strength_TESTEX <- 100 * MultiLevelModel_logCVR_strength_TESTEX$sigma2 / (sum(MultiLevelModel_logCVR_strength_TESTEX$sigma2) + (MultiLevelModel_logCVR_strength_TESTEX$k-MultiLevelModel_logCVR_strength_TESTEX$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_TESTEX <- robust(MultiLevelModel_logCVR_strength_TESTEX, subset(Data_logCVR_TESTEX, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_TESTEX, file = "models/RobuEstMultiLevelModel_logCVR_strength_TESTEX")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_TESTEX <- rma.mv(yi, V=vi, data=subset(Data_logCVR_TESTEX, outcome == "hypertrophy"),
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ TESTEX, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_TESTEX, file = "models/MultiLevelModel_logCVR_hypertrophy_TESTEX")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_TESTEX, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_TESTEX)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_TESTEX <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_TESTEX$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_TESTEX$sigma2) + (MultiLevelModel_logCVR_hypertrophy_TESTEX$k-MultiLevelModel_logCVR_hypertrophy_TESTEX$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_TESTEX <- 100 * MultiLevelModel_logCVR_hypertrophy_TESTEX$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_TESTEX$sigma2) + (MultiLevelModel_logCVR_hypertrophy_TESTEX$k-MultiLevelModel_logCVR_hypertrophy_TESTEX$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX <- robust(MultiLevelModel_logCVR_hypertrophy_TESTEX, subset(Data_logCVR_TESTEX, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX")

### Mean age
# SMD
# Strength
Data_SMD_strength_age <- Data_SMD_strength %>%
  filter(!is.na(age))

MultiLevelModel_strength_age <- rma.mv(yi, V=vi, data=Data_SMD_strength_age,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ age, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_age, file = "models/MultiLevelModel_strength_age")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_age$vi)
X <- model.matrix(MultiLevelModel_strength_age)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_age <- 100 * sum(MultiLevelModel_strength_age$sigma2) / (sum(MultiLevelModel_strength_age$sigma2) + (MultiLevelModel_strength_age$k-MultiLevelModel_strength_age$p)/sum(diag(P)))
I2bw_strength_age <- 100 * MultiLevelModel_strength_age$sigma2 / (sum(MultiLevelModel_strength_age$sigma2) + (MultiLevelModel_strength_age$k-MultiLevelModel_strength_age$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_age <- robust(MultiLevelModel_strength_age, Data_SMD_strength_age$study)

save(RobuEstMultiLevelModel_strength_age, file = "models/RobuEstMultiLevelModel_strength_age")

# Hypertrophy
Data_SMD_hypertrophy_age <- Data_SMD_hypertrophy %>%
  filter(!is.na(age))

MultiLevelModel_hypertrophy_age <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_age,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ age, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_age, file = "models/MultiLevelModel_hypertrophy_age")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_age$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_age)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_age <- 100 * sum(MultiLevelModel_hypertrophy_age$sigma2) / (sum(MultiLevelModel_hypertrophy_age$sigma2) + (MultiLevelModel_hypertrophy_age$k-MultiLevelModel_hypertrophy_age$p)/sum(diag(P)))
I2bw_hypertrophy_age <- 100 * MultiLevelModel_hypertrophy_age$sigma2 / (sum(MultiLevelModel_hypertrophy_age$sigma2) + (MultiLevelModel_hypertrophy_age$k-MultiLevelModel_hypertrophy_age$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_age <- robust(MultiLevelModel_hypertrophy_age, Data_SMD_hypertrophy_age$study)

save(RobuEstMultiLevelModel_hypertrophy_age, file = "models/RobuEstMultiLevelModel_hypertrophy_age")

# logCVR
Data_logCVR_age <- Data_logCVR %>%
  filter(!is.na(age))

# Strength
MultiLevelModel_logCVR_strength_age <- rma.mv(yi, V=vi, data=subset(Data_logCVR_age, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ age, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_age, file = "models/MultiLevelModel_logCVR_strength_age")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_age, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_age)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_age <- 100 * sum(MultiLevelModel_logCVR_strength_age$sigma2) / (sum(MultiLevelModel_logCVR_strength_age$sigma2) + (MultiLevelModel_logCVR_strength_age$k-MultiLevelModel_logCVR_strength_age$p)/sum(diag(P)))
I2bw_logCVR_strength_age <- 100 * MultiLevelModel_logCVR_strength_age$sigma2 / (sum(MultiLevelModel_logCVR_strength_age$sigma2) + (MultiLevelModel_logCVR_strength_age$k-MultiLevelModel_logCVR_strength_age$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_age <- robust(MultiLevelModel_logCVR_strength_age, subset(Data_logCVR_age, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_age, file = "models/RobuEstMultiLevelModel_logCVR_strength_age")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_age <- rma.mv(yi, V=vi, data=subset(Data_logCVR_age, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ age, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_age, file = "models/MultiLevelModel_logCVR_hypertrophy_age")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_age, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_age)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_age <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_age$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_age$sigma2) + (MultiLevelModel_logCVR_hypertrophy_age$k-MultiLevelModel_logCVR_hypertrophy_age$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_age <- 100 * MultiLevelModel_logCVR_hypertrophy_age$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_age$sigma2) + (MultiLevelModel_logCVR_hypertrophy_age$k-MultiLevelModel_logCVR_hypertrophy_age$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_age <- robust(MultiLevelModel_logCVR_hypertrophy_age, subset(Data_logCVR_age, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_age, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_age")

### Proportion of males
# SMD
# Strength
Data_SMD_strength_sex_._male <- Data_SMD_strength %>%
  filter(!is.na(sex_._male))

MultiLevelModel_strength_sex_._male <- rma.mv(yi, V=vi, data=Data_SMD_strength_sex_._male,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ sex_._male, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_sex_._male, file = "models/MultiLevelModel_strength_sex_._male")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_sex_._male$vi)
X <- model.matrix(MultiLevelModel_strength_sex_._male)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_sex_._male <- 100 * sum(MultiLevelModel_strength_sex_._male$sigma2) / (sum(MultiLevelModel_strength_sex_._male$sigma2) + (MultiLevelModel_strength_sex_._male$k-MultiLevelModel_strength_sex_._male$p)/sum(diag(P)))
I2bw_strength_sex_._male <- 100 * MultiLevelModel_strength_sex_._male$sigma2 / (sum(MultiLevelModel_strength_sex_._male$sigma2) + (MultiLevelModel_strength_sex_._male$k-MultiLevelModel_strength_sex_._male$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_sex_._male <- robust(MultiLevelModel_strength_sex_._male, Data_SMD_strength_sex_._male$study)

save(RobuEstMultiLevelModel_strength_sex_._male, file = "models/RobuEstMultiLevelModel_strength_sex_._male")

# Hypertrophy
Data_SMD_hypertrophy_sex_._male <- Data_SMD_hypertrophy %>%
  filter(!is.na(sex_._male))

MultiLevelModel_hypertrophy_sex_._male <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_sex_._male,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ sex_._male, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_sex_._male, file = "models/MultiLevelModel_hypertrophy_sex_._male")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_sex_._male$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_sex_._male)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_sex_._male <- 100 * sum(MultiLevelModel_hypertrophy_sex_._male$sigma2) / (sum(MultiLevelModel_hypertrophy_sex_._male$sigma2) + (MultiLevelModel_hypertrophy_sex_._male$k-MultiLevelModel_hypertrophy_sex_._male$p)/sum(diag(P)))
I2bw_hypertrophy_sex_._male <- 100 * MultiLevelModel_hypertrophy_sex_._male$sigma2 / (sum(MultiLevelModel_hypertrophy_sex_._male$sigma2) + (MultiLevelModel_hypertrophy_sex_._male$k-MultiLevelModel_hypertrophy_sex_._male$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_sex_._male <- robust(MultiLevelModel_hypertrophy_sex_._male, Data_SMD_hypertrophy_sex_._male$study)

save(RobuEstMultiLevelModel_hypertrophy_sex_._male, file = "models/RobuEstMultiLevelModel_hypertrophy_sex_._male")

# logCVR
Data_logCVR_sex_._male <- Data_logCVR %>%
  filter(!is.na(sex_._male))

# Strength
MultiLevelModel_logCVR_strength_sex_._male <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sex_._male, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ sex_._male, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_sex_._male, file = "models/MultiLevelModel_logCVR_strength_sex_._male")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_sex_._male, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_sex_._male)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_sex_._male <- 100 * sum(MultiLevelModel_logCVR_strength_sex_._male$sigma2) / (sum(MultiLevelModel_logCVR_strength_sex_._male$sigma2) + (MultiLevelModel_logCVR_strength_sex_._male$k-MultiLevelModel_logCVR_strength_sex_._male$p)/sum(diag(P)))
I2bw_logCVR_strength_sex_._male <- 100 * MultiLevelModel_logCVR_strength_sex_._male$sigma2 / (sum(MultiLevelModel_logCVR_strength_sex_._male$sigma2) + (MultiLevelModel_logCVR_strength_sex_._male$k-MultiLevelModel_logCVR_strength_sex_._male$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_sex_._male <- robust(MultiLevelModel_logCVR_strength_sex_._male, subset(Data_logCVR_sex_._male, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_sex_._male, file = "models/RobuEstMultiLevelModel_logCVR_strength_sex_._male")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_sex_._male <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sex_._male, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ sex_._male, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_sex_._male, file = "models/MultiLevelModel_logCVR_hypertrophy_sex_._male")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_sex_._male, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_sex_._male)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_sex_._male <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_sex_._male$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_sex_._male$sigma2) + (MultiLevelModel_logCVR_hypertrophy_sex_._male$k-MultiLevelModel_logCVR_hypertrophy_sex_._male$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_sex_._male <- 100 * MultiLevelModel_logCVR_hypertrophy_sex_._male$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_sex_._male$sigma2) + (MultiLevelModel_logCVR_hypertrophy_sex_._male$k-MultiLevelModel_logCVR_hypertrophy_sex_._male$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male <- robust(MultiLevelModel_logCVR_hypertrophy_sex_._male, subset(Data_logCVR_sex_._male, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male")

### Mean weight
# SMD
# Strength
Data_SMD_strength_weight <- Data_SMD_strength %>%
  filter(!is.na(weight))

MultiLevelModel_strength_weight <- rma.mv(yi, V=vi, data=Data_SMD_strength_weight,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ weight, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_weight, file = "models/MultiLevelModel_strength_weight")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_weight$vi)
X <- model.matrix(MultiLevelModel_strength_weight)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_weight <- 100 * sum(MultiLevelModel_strength_weight$sigma2) / (sum(MultiLevelModel_strength_weight$sigma2) + (MultiLevelModel_strength_weight$k-MultiLevelModel_strength_weight$p)/sum(diag(P)))
I2bw_strength_weight <- 100 * MultiLevelModel_strength_weight$sigma2 / (sum(MultiLevelModel_strength_weight$sigma2) + (MultiLevelModel_strength_weight$k-MultiLevelModel_strength_weight$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_weight <- robust(MultiLevelModel_strength_weight, Data_SMD_strength_weight$study)

save(RobuEstMultiLevelModel_strength_weight, file = "models/RobuEstMultiLevelModel_strength_weight")

# Hypertrophy
Data_SMD_hypertrophy_weight <- Data_SMD_hypertrophy %>%
  filter(!is.na(weight))

MultiLevelModel_hypertrophy_weight <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_weight,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ weight, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_weight, file = "models/MultiLevelModel_hypertrophy_weight")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_weight$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_weight)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_weight <- 100 * sum(MultiLevelModel_hypertrophy_weight$sigma2) / (sum(MultiLevelModel_hypertrophy_weight$sigma2) + (MultiLevelModel_hypertrophy_weight$k-MultiLevelModel_hypertrophy_weight$p)/sum(diag(P)))
I2bw_hypertrophy_weight <- 100 * MultiLevelModel_hypertrophy_weight$sigma2 / (sum(MultiLevelModel_hypertrophy_weight$sigma2) + (MultiLevelModel_hypertrophy_weight$k-MultiLevelModel_hypertrophy_weight$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_weight <- robust(MultiLevelModel_hypertrophy_weight, Data_SMD_hypertrophy_weight$study)

save(RobuEstMultiLevelModel_hypertrophy_weight, file = "models/RobuEstMultiLevelModel_hypertrophy_weight")

# logCVR
Data_logCVR_weight <- Data_logCVR %>%
  filter(!is.na(weight))

# Strength
MultiLevelModel_logCVR_strength_weight <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weight, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ weight, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_weight, file = "models/MultiLevelModel_logCVR_strength_weight")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_weight, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_weight)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_weight <- 100 * sum(MultiLevelModel_logCVR_strength_weight$sigma2) / (sum(MultiLevelModel_logCVR_strength_weight$sigma2) + (MultiLevelModel_logCVR_strength_weight$k-MultiLevelModel_logCVR_strength_weight$p)/sum(diag(P)))
I2bw_logCVR_strength_weight <- 100 * MultiLevelModel_logCVR_strength_weight$sigma2 / (sum(MultiLevelModel_logCVR_strength_weight$sigma2) + (MultiLevelModel_logCVR_strength_weight$k-MultiLevelModel_logCVR_strength_weight$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_weight <- robust(MultiLevelModel_logCVR_strength_weight, subset(Data_logCVR_weight, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_weight, file = "models/RobuEstMultiLevelModel_logCVR_strength_weight")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_weight <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weight, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ weight, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_weight, file = "models/MultiLevelModel_logCVR_hypertrophy_weight")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_weight, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_weight)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_weight <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_weight$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_weight$sigma2) + (MultiLevelModel_logCVR_hypertrophy_weight$k-MultiLevelModel_logCVR_hypertrophy_weight$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_weight <- 100 * MultiLevelModel_logCVR_hypertrophy_weight$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_weight$sigma2) + (MultiLevelModel_logCVR_hypertrophy_weight$k-MultiLevelModel_logCVR_hypertrophy_weight$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_weight <- robust(MultiLevelModel_logCVR_hypertrophy_weight, subset(Data_logCVR_weight, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_weight, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_weight")

### Mean BMI
# SMD
# Strength
Data_SMD_strength_bmi <- Data_SMD_strength %>%
  filter(!is.na(bmi))

MultiLevelModel_strength_bmi <- rma.mv(yi, V=vi, data=Data_SMD_strength_bmi,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ bmi, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_bmi, file = "models/MultiLevelModel_strength_bmi")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_bmi$vi)
X <- model.matrix(MultiLevelModel_strength_bmi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_bmi <- 100 * sum(MultiLevelModel_strength_bmi$sigma2) / (sum(MultiLevelModel_strength_bmi$sigma2) + (MultiLevelModel_strength_bmi$k-MultiLevelModel_strength_bmi$p)/sum(diag(P)))
I2bw_strength_bmi <- 100 * MultiLevelModel_strength_bmi$sigma2 / (sum(MultiLevelModel_strength_bmi$sigma2) + (MultiLevelModel_strength_bmi$k-MultiLevelModel_strength_bmi$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_bmi <- robust(MultiLevelModel_strength_bmi, Data_SMD_strength_bmi$study)

save(RobuEstMultiLevelModel_strength_bmi, file = "models/RobuEstMultiLevelModel_strength_bmi")

# Hypertrophy
Data_SMD_hypertrophy_bmi <- Data_SMD_hypertrophy %>%
  filter(!is.na(bmi))

MultiLevelModel_hypertrophy_bmi <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_bmi,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ bmi, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_bmi, file = "models/MultiLevelModel_hypertrophy_bmi")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_bmi$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_bmi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_bmi <- 100 * sum(MultiLevelModel_hypertrophy_bmi$sigma2) / (sum(MultiLevelModel_hypertrophy_bmi$sigma2) + (MultiLevelModel_hypertrophy_bmi$k-MultiLevelModel_hypertrophy_bmi$p)/sum(diag(P)))
I2bw_hypertrophy_bmi <- 100 * MultiLevelModel_hypertrophy_bmi$sigma2 / (sum(MultiLevelModel_hypertrophy_bmi$sigma2) + (MultiLevelModel_hypertrophy_bmi$k-MultiLevelModel_hypertrophy_bmi$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_bmi <- robust(MultiLevelModel_hypertrophy_bmi, Data_SMD_hypertrophy_bmi$study)

save(RobuEstMultiLevelModel_hypertrophy_bmi, file = "models/RobuEstMultiLevelModel_hypertrophy_bmi")

# logCVR
Data_logCVR_bmi <- Data_logCVR %>%
  filter(!is.na(bmi))

# Strength
MultiLevelModel_logCVR_strength_bmi <- rma.mv(yi, V=vi, data=subset(Data_logCVR_bmi, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ bmi, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_bmi, file = "models/MultiLevelModel_logCVR_strength_bmi")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_bmi, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_bmi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_bmi <- 100 * sum(MultiLevelModel_logCVR_strength_bmi$sigma2) / (sum(MultiLevelModel_logCVR_strength_bmi$sigma2) + (MultiLevelModel_logCVR_strength_bmi$k-MultiLevelModel_logCVR_strength_bmi$p)/sum(diag(P)))
I2bw_logCVR_strength_bmi <- 100 * MultiLevelModel_logCVR_strength_bmi$sigma2 / (sum(MultiLevelModel_logCVR_strength_bmi$sigma2) + (MultiLevelModel_logCVR_strength_bmi$k-MultiLevelModel_logCVR_strength_bmi$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_bmi <- robust(MultiLevelModel_logCVR_strength_bmi, subset(Data_logCVR_bmi, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_bmi, file = "models/RobuEstMultiLevelModel_logCVR_strength_bmi")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_bmi <- rma.mv(yi, V=vi, data=subset(Data_logCVR_bmi, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ bmi, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_bmi, file = "models/MultiLevelModel_logCVR_hypertrophy_bmi")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_bmi, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_bmi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_bmi <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_bmi$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_bmi$sigma2) + (MultiLevelModel_logCVR_hypertrophy_bmi$k-MultiLevelModel_logCVR_hypertrophy_bmi$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_bmi <- 100 * MultiLevelModel_logCVR_hypertrophy_bmi$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_bmi$sigma2) + (MultiLevelModel_logCVR_hypertrophy_bmi$k-MultiLevelModel_logCVR_hypertrophy_bmi$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_bmi <- robust(MultiLevelModel_logCVR_hypertrophy_bmi, subset(Data_logCVR_bmi, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_bmi")

### Trained vs untrained
# SMD
# Strength
Data_SMD_strength_train_status <- Data_SMD_strength %>%
  filter(!is.na(train_status))

MultiLevelModel_strength_train_status <- rma.mv(yi, V=vi, data=Data_SMD_strength_train_status,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ train_status - 1, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_train_status, file = "models/MultiLevelModel_strength_train_status")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_train_status$vi)
X <- model.matrix(MultiLevelModel_strength_train_status)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_train_status <- 100 * sum(MultiLevelModel_strength_train_status$sigma2) / (sum(MultiLevelModel_strength_train_status$sigma2) + (MultiLevelModel_strength_train_status$k-MultiLevelModel_strength_train_status$p)/sum(diag(P)))
I2bw_strength_train_status <- 100 * MultiLevelModel_strength_train_status$sigma2 / (sum(MultiLevelModel_strength_train_status$sigma2) + (MultiLevelModel_strength_train_status$k-MultiLevelModel_strength_train_status$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_train_status <- robust(MultiLevelModel_strength_train_status, Data_SMD_strength_train_status$study)

save(RobuEstMultiLevelModel_strength_train_status, file = "models/RobuEstMultiLevelModel_strength_train_status")

# Hypertrophy
Data_SMD_hypertrophy_train_status <- Data_SMD_hypertrophy %>%
  filter(!is.na(train_status))

MultiLevelModel_hypertrophy_train_status <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_train_status,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ train_status - 1, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_train_status, file = "models/MultiLevelModel_hypertrophy_train_status")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_train_status$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_train_status)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_train_status <- 100 * sum(MultiLevelModel_hypertrophy_train_status$sigma2) / (sum(MultiLevelModel_hypertrophy_train_status$sigma2) + (MultiLevelModel_hypertrophy_train_status$k-MultiLevelModel_hypertrophy_train_status$p)/sum(diag(P)))
I2bw_hypertrophy_train_status <- 100 * MultiLevelModel_hypertrophy_train_status$sigma2 / (sum(MultiLevelModel_hypertrophy_train_status$sigma2) + (MultiLevelModel_hypertrophy_train_status$k-MultiLevelModel_hypertrophy_train_status$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_train_status <- robust(MultiLevelModel_hypertrophy_train_status, Data_SMD_hypertrophy_train_status$study)

save(RobuEstMultiLevelModel_hypertrophy_train_status, file = "models/RobuEstMultiLevelModel_hypertrophy_train_status")

# logCVR
Data_logCVR_train_status <- Data_logCVR %>%
  filter(!is.na(train_status))

# Strength
MultiLevelModel_logCVR_strength_train_status <- rma.mv(yi, V=vi, data=subset(Data_logCVR_train_status, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ train_status - 1, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_train_status, file = "models/MultiLevelModel_logCVR_strength_train_status")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_train_status, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_train_status)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_train_status <- 100 * sum(MultiLevelModel_logCVR_strength_train_status$sigma2) / (sum(MultiLevelModel_logCVR_strength_train_status$sigma2) + (MultiLevelModel_logCVR_strength_train_status$k-MultiLevelModel_logCVR_strength_train_status$p)/sum(diag(P)))
I2bw_logCVR_strength_train_status <- 100 * MultiLevelModel_logCVR_strength_train_status$sigma2 / (sum(MultiLevelModel_logCVR_strength_train_status$sigma2) + (MultiLevelModel_logCVR_strength_train_status$k-MultiLevelModel_logCVR_strength_train_status$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_train_status <- robust(MultiLevelModel_logCVR_strength_train_status, subset(Data_logCVR_train_status, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_train_status, file = "models/RobuEstMultiLevelModel_logCVR_strength_train_status")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_train_status <- rma.mv(yi, V=vi, data=subset(Data_logCVR_train_status, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ train_status - 1, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_train_status, file = "models/MultiLevelModel_logCVR_hypertrophy_train_status")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_train_status, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_train_status)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_train_status <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_train_status$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_train_status$sigma2) + (MultiLevelModel_logCVR_hypertrophy_train_status$k-MultiLevelModel_logCVR_hypertrophy_train_status$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_train_status <- 100 * MultiLevelModel_logCVR_hypertrophy_train_status$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_train_status$sigma2) + (MultiLevelModel_logCVR_hypertrophy_train_status$k-MultiLevelModel_logCVR_hypertrophy_train_status$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_train_status <- robust(MultiLevelModel_logCVR_hypertrophy_train_status, subset(Data_logCVR_train_status, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_train_status")

### Healthy vs clinical sample
# SMD
# Strength
Data_SMD_strength_healthy_clinical <- Data_SMD_strength %>%
  filter(!is.na(healthy_clinical))

MultiLevelModel_strength_healthy_clinical <- rma.mv(yi, V=vi, data=Data_SMD_strength_healthy_clinical,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ healthy_clinical - 1, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_healthy_clinical, file = "models/MultiLevelModel_strength_healthy_clinical")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_healthy_clinical$vi)
X <- model.matrix(MultiLevelModel_strength_healthy_clinical)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_healthy_clinical <- 100 * sum(MultiLevelModel_strength_healthy_clinical$sigma2) / (sum(MultiLevelModel_strength_healthy_clinical$sigma2) + (MultiLevelModel_strength_healthy_clinical$k-MultiLevelModel_strength_healthy_clinical$p)/sum(diag(P)))
I2bw_strength_healthy_clinical <- 100 * MultiLevelModel_strength_healthy_clinical$sigma2 / (sum(MultiLevelModel_strength_healthy_clinical$sigma2) + (MultiLevelModel_strength_healthy_clinical$k-MultiLevelModel_strength_healthy_clinical$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_healthy_clinical <- robust(MultiLevelModel_strength_healthy_clinical, Data_SMD_strength_healthy_clinical$study)

save(RobuEstMultiLevelModel_strength_healthy_clinical, file = "models/RobuEstMultiLevelModel_strength_healthy_clinical")

# Hypertrophy
Data_SMD_hypertrophy_healthy_clinical <- Data_SMD_hypertrophy %>%
  filter(!is.na(healthy_clinical))

MultiLevelModel_hypertrophy_healthy_clinical <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_healthy_clinical,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ healthy_clinical - 1, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_healthy_clinical, file = "models/MultiLevelModel_hypertrophy_healthy_clinical")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_healthy_clinical$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_healthy_clinical)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_healthy_clinical <- 100 * sum(MultiLevelModel_hypertrophy_healthy_clinical$sigma2) / (sum(MultiLevelModel_hypertrophy_healthy_clinical$sigma2) + (MultiLevelModel_hypertrophy_healthy_clinical$k-MultiLevelModel_hypertrophy_healthy_clinical$p)/sum(diag(P)))
I2bw_hypertrophy_healthy_clinical <- 100 * MultiLevelModel_hypertrophy_healthy_clinical$sigma2 / (sum(MultiLevelModel_hypertrophy_healthy_clinical$sigma2) + (MultiLevelModel_hypertrophy_healthy_clinical$k-MultiLevelModel_hypertrophy_healthy_clinical$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_healthy_clinical <- robust(MultiLevelModel_hypertrophy_healthy_clinical, Data_SMD_hypertrophy_healthy_clinical$study)

save(RobuEstMultiLevelModel_hypertrophy_healthy_clinical, file = "models/RobuEstMultiLevelModel_hypertrophy_healthy_clinical")

# logCVR
Data_logCVR_healthy_clinical <- Data_logCVR %>%
  filter(!is.na(healthy_clinical))

# Strength
MultiLevelModel_logCVR_strength_healthy_clinical <- rma.mv(yi, V=vi, data=subset(Data_logCVR_healthy_clinical, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ healthy_clinical - 1, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_healthy_clinical, file = "models/MultiLevelModel_logCVR_strength_healthy_clinical")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_healthy_clinical, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_healthy_clinical)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_healthy_clinical <- 100 * sum(MultiLevelModel_logCVR_strength_healthy_clinical$sigma2) / (sum(MultiLevelModel_logCVR_strength_healthy_clinical$sigma2) + (MultiLevelModel_logCVR_strength_healthy_clinical$k-MultiLevelModel_logCVR_strength_healthy_clinical$p)/sum(diag(P)))
I2bw_logCVR_strength_healthy_clinical <- 100 * MultiLevelModel_logCVR_strength_healthy_clinical$sigma2 / (sum(MultiLevelModel_logCVR_strength_healthy_clinical$sigma2) + (MultiLevelModel_logCVR_strength_healthy_clinical$k-MultiLevelModel_logCVR_strength_healthy_clinical$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_healthy_clinical <- robust(MultiLevelModel_logCVR_strength_healthy_clinical, subset(Data_logCVR_healthy_clinical, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical, file = "models/RobuEstMultiLevelModel_logCVR_strength_healthy_clinical")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_healthy_clinical <- rma.mv(yi, V=vi, data=subset(Data_logCVR_healthy_clinical, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ healthy_clinical - 1, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_healthy_clinical, file = "models/MultiLevelModel_logCVR_hypertrophy_healthy_clinical")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_healthy_clinical, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_healthy_clinical)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_healthy_clinical <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_healthy_clinical$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_healthy_clinical$sigma2) + (MultiLevelModel_logCVR_hypertrophy_healthy_clinical$k-MultiLevelModel_logCVR_hypertrophy_healthy_clinical$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_healthy_clinical <- 100 * MultiLevelModel_logCVR_hypertrophy_healthy_clinical$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_healthy_clinical$sigma2) + (MultiLevelModel_logCVR_hypertrophy_healthy_clinical$k-MultiLevelModel_logCVR_hypertrophy_healthy_clinical$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical <- robust(MultiLevelModel_logCVR_hypertrophy_healthy_clinical, subset(Data_logCVR_healthy_clinical, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical")

### Only RT intervention (y), or additional e.g., supplements, aerobic (n)
# SMD
# Strength
Data_SMD_strength_RT_only <- Data_SMD_strength %>%
  filter(!is.na(RT_only))

MultiLevelModel_strength_RT_only <- rma.mv(yi, V=vi, data=Data_SMD_strength_RT_only,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ RT_only - 1, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_RT_only, file = "models/MultiLevelModel_strength_RT_only")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_RT_only$vi)
X <- model.matrix(MultiLevelModel_strength_RT_only)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_RT_only <- 100 * sum(MultiLevelModel_strength_RT_only$sigma2) / (sum(MultiLevelModel_strength_RT_only$sigma2) + (MultiLevelModel_strength_RT_only$k-MultiLevelModel_strength_RT_only$p)/sum(diag(P)))
I2bw_strength_RT_only <- 100 * MultiLevelModel_strength_RT_only$sigma2 / (sum(MultiLevelModel_strength_RT_only$sigma2) + (MultiLevelModel_strength_RT_only$k-MultiLevelModel_strength_RT_only$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_RT_only <- robust(MultiLevelModel_strength_RT_only, Data_SMD_strength_RT_only$study)

save(RobuEstMultiLevelModel_strength_RT_only, file = "models/RobuEstMultiLevelModel_strength_RT_only")

# Hypertrophy
Data_SMD_hypertrophy_RT_only <- Data_SMD_hypertrophy %>%
  filter(!is.na(RT_only))

MultiLevelModel_hypertrophy_RT_only <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_RT_only,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ RT_only - 1, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_RT_only, file = "models/MultiLevelModel_hypertrophy_RT_only")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_RT_only$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_RT_only)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_RT_only <- 100 * sum(MultiLevelModel_hypertrophy_RT_only$sigma2) / (sum(MultiLevelModel_hypertrophy_RT_only$sigma2) + (MultiLevelModel_hypertrophy_RT_only$k-MultiLevelModel_hypertrophy_RT_only$p)/sum(diag(P)))
I2bw_hypertrophy_RT_only <- 100 * MultiLevelModel_hypertrophy_RT_only$sigma2 / (sum(MultiLevelModel_hypertrophy_RT_only$sigma2) + (MultiLevelModel_hypertrophy_RT_only$k-MultiLevelModel_hypertrophy_RT_only$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_RT_only <- robust(MultiLevelModel_hypertrophy_RT_only, Data_SMD_hypertrophy_RT_only$study)

save(RobuEstMultiLevelModel_hypertrophy_RT_only, file = "models/RobuEstMultiLevelModel_hypertrophy_RT_only")

# logCVR
Data_logCVR_RT_only <- Data_logCVR %>%
  filter(!is.na(RT_only))

# Strength
MultiLevelModel_logCVR_strength_RT_only <- rma.mv(yi, V=vi, data=subset(Data_logCVR_RT_only, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ RT_only - 1, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_RT_only, file = "models/MultiLevelModel_logCVR_strength_RT_only")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_RT_only, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_RT_only)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_RT_only <- 100 * sum(MultiLevelModel_logCVR_strength_RT_only$sigma2) / (sum(MultiLevelModel_logCVR_strength_RT_only$sigma2) + (MultiLevelModel_logCVR_strength_RT_only$k-MultiLevelModel_logCVR_strength_RT_only$p)/sum(diag(P)))
I2bw_logCVR_strength_RT_only <- 100 * MultiLevelModel_logCVR_strength_RT_only$sigma2 / (sum(MultiLevelModel_logCVR_strength_RT_only$sigma2) + (MultiLevelModel_logCVR_strength_RT_only$k-MultiLevelModel_logCVR_strength_RT_only$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_RT_only <- robust(MultiLevelModel_logCVR_strength_RT_only, subset(Data_logCVR_RT_only, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_RT_only, file = "models/RobuEstMultiLevelModel_logCVR_strength_RT_only")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_RT_only <- rma.mv(yi, V=vi, data=subset(Data_logCVR_RT_only, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ RT_only - 1, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_RT_only, file = "models/MultiLevelModel_logCVR_hypertrophy_RT_only")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_RT_only, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_RT_only)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_RT_only <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_RT_only$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_RT_only$sigma2) + (MultiLevelModel_logCVR_hypertrophy_RT_only$k-MultiLevelModel_logCVR_hypertrophy_RT_only$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_RT_only <- 100 * MultiLevelModel_logCVR_hypertrophy_RT_only$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_RT_only$sigma2) + (MultiLevelModel_logCVR_hypertrophy_RT_only$k-MultiLevelModel_logCVR_hypertrophy_RT_only$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only <- robust(MultiLevelModel_logCVR_hypertrophy_RT_only, subset(Data_logCVR_RT_only, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only")

### Duration of intervention in weeks
# SMD
# Strength
Data_SMD_strength_weeks <- Data_SMD_strength %>%
  filter(!is.na(weeks))

MultiLevelModel_strength_weeks <- rma.mv(yi, V=vi, data=Data_SMD_strength_weeks,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ weeks, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_weeks, file = "models/MultiLevelModel_strength_weeks")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_weeks$vi)
X <- model.matrix(MultiLevelModel_strength_weeks)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_weeks <- 100 * sum(MultiLevelModel_strength_weeks$sigma2) / (sum(MultiLevelModel_strength_weeks$sigma2) + (MultiLevelModel_strength_weeks$k-MultiLevelModel_strength_weeks$p)/sum(diag(P)))
I2bw_strength_weeks <- 100 * MultiLevelModel_strength_weeks$sigma2 / (sum(MultiLevelModel_strength_weeks$sigma2) + (MultiLevelModel_strength_weeks$k-MultiLevelModel_strength_weeks$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_weeks <- robust(MultiLevelModel_strength_weeks, Data_SMD_strength_weeks$study)

save(RobuEstMultiLevelModel_strength_weeks, file = "models/RobuEstMultiLevelModel_strength_weeks")

# Hypertrophy
Data_SMD_hypertrophy_weeks <- Data_SMD_hypertrophy %>%
  filter(!is.na(weeks))

MultiLevelModel_hypertrophy_weeks <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_weeks,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ weeks, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_weeks, file = "models/MultiLevelModel_hypertrophy_weeks")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_weeks$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_weeks)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_weeks <- 100 * sum(MultiLevelModel_hypertrophy_weeks$sigma2) / (sum(MultiLevelModel_hypertrophy_weeks$sigma2) + (MultiLevelModel_hypertrophy_weeks$k-MultiLevelModel_hypertrophy_weeks$p)/sum(diag(P)))
I2bw_hypertrophy_weeks <- 100 * MultiLevelModel_hypertrophy_weeks$sigma2 / (sum(MultiLevelModel_hypertrophy_weeks$sigma2) + (MultiLevelModel_hypertrophy_weeks$k-MultiLevelModel_hypertrophy_weeks$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_weeks <- robust(MultiLevelModel_hypertrophy_weeks, Data_SMD_hypertrophy_weeks$study)

save(RobuEstMultiLevelModel_hypertrophy_weeks, file = "models/RobuEstMultiLevelModel_hypertrophy_weeks")

# logCVR
Data_logCVR_weeks <- Data_logCVR %>%
  filter(!is.na(weeks))

# Strength
MultiLevelModel_logCVR_strength_weeks <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weeks, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ weeks, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_weeks, file = "models/MultiLevelModel_logCVR_strength_weeks")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_weeks, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_weeks)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_weeks <- 100 * sum(MultiLevelModel_logCVR_strength_weeks$sigma2) / (sum(MultiLevelModel_logCVR_strength_weeks$sigma2) + (MultiLevelModel_logCVR_strength_weeks$k-MultiLevelModel_logCVR_strength_weeks$p)/sum(diag(P)))
I2bw_logCVR_strength_weeks <- 100 * MultiLevelModel_logCVR_strength_weeks$sigma2 / (sum(MultiLevelModel_logCVR_strength_weeks$sigma2) + (MultiLevelModel_logCVR_strength_weeks$k-MultiLevelModel_logCVR_strength_weeks$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_weeks <- robust(MultiLevelModel_logCVR_strength_weeks, subset(Data_logCVR_weeks, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_weeks, file = "models/RobuEstMultiLevelModel_logCVR_strength_weeks")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_weeks <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weeks, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ weeks, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_weeks, file = "models/MultiLevelModel_logCVR_hypertrophy_weeks")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_weeks, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_weeks)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_weeks <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_weeks$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_weeks$sigma2) + (MultiLevelModel_logCVR_hypertrophy_weeks$k-MultiLevelModel_logCVR_hypertrophy_weeks$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_weeks <- 100 * MultiLevelModel_logCVR_hypertrophy_weeks$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_weeks$sigma2) + (MultiLevelModel_logCVR_hypertrophy_weeks$k-MultiLevelModel_logCVR_hypertrophy_weeks$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_weeks <- robust(MultiLevelModel_logCVR_hypertrophy_weeks, subset(Data_logCVR_weeks, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_weeks")

### Weekly frequency in days
# SMD
# Strength
Data_SMD_strength_freq <- Data_SMD_strength %>%
  filter(!is.na(freq))

MultiLevelModel_strength_freq <- rma.mv(yi, V=vi, data=Data_SMD_strength_freq,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ freq, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_freq, file = "models/MultiLevelModel_strength_freq")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_freq$vi)
X <- model.matrix(MultiLevelModel_strength_freq)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_freq <- 100 * sum(MultiLevelModel_strength_freq$sigma2) / (sum(MultiLevelModel_strength_freq$sigma2) + (MultiLevelModel_strength_freq$k-MultiLevelModel_strength_freq$p)/sum(diag(P)))
I2bw_strength_freq <- 100 * MultiLevelModel_strength_freq$sigma2 / (sum(MultiLevelModel_strength_freq$sigma2) + (MultiLevelModel_strength_freq$k-MultiLevelModel_strength_freq$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_freq <- robust(MultiLevelModel_strength_freq, Data_SMD_strength_freq$study)

save(RobuEstMultiLevelModel_strength_freq, file = "models/RobuEstMultiLevelModel_strength_freq")

# Hypertrophy
Data_SMD_hypertrophy_freq <- Data_SMD_hypertrophy %>%
  filter(!is.na(freq))

MultiLevelModel_hypertrophy_freq <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_freq,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ freq, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_freq, file = "models/MultiLevelModel_hypertrophy_freq")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_freq$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_freq)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_freq <- 100 * sum(MultiLevelModel_hypertrophy_freq$sigma2) / (sum(MultiLevelModel_hypertrophy_freq$sigma2) + (MultiLevelModel_hypertrophy_freq$k-MultiLevelModel_hypertrophy_freq$p)/sum(diag(P)))
I2bw_hypertrophy_freq <- 100 * MultiLevelModel_hypertrophy_freq$sigma2 / (sum(MultiLevelModel_hypertrophy_freq$sigma2) + (MultiLevelModel_hypertrophy_freq$k-MultiLevelModel_hypertrophy_freq$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_freq <- robust(MultiLevelModel_hypertrophy_freq, Data_SMD_hypertrophy_freq$study)

save(RobuEstMultiLevelModel_hypertrophy_freq, file = "models/RobuEstMultiLevelModel_hypertrophy_freq")

# logCVR
Data_logCVR_freq <- Data_logCVR %>%
  filter(!is.na(freq))

# Strength
MultiLevelModel_logCVR_strength_freq <- rma.mv(yi, V=vi, data=subset(Data_logCVR_freq, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ freq, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_freq, file = "models/MultiLevelModel_logCVR_strength_freq")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_freq, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_freq)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_freq <- 100 * sum(MultiLevelModel_logCVR_strength_freq$sigma2) / (sum(MultiLevelModel_logCVR_strength_freq$sigma2) + (MultiLevelModel_logCVR_strength_freq$k-MultiLevelModel_logCVR_strength_freq$p)/sum(diag(P)))
I2bw_logCVR_strength_freq <- 100 * MultiLevelModel_logCVR_strength_freq$sigma2 / (sum(MultiLevelModel_logCVR_strength_freq$sigma2) + (MultiLevelModel_logCVR_strength_freq$k-MultiLevelModel_logCVR_strength_freq$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_freq <- robust(MultiLevelModel_logCVR_strength_freq, subset(Data_logCVR_freq, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_freq, file = "models/RobuEstMultiLevelModel_logCVR_strength_freq")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_freq <- rma.mv(yi, V=vi, data=subset(Data_logCVR_freq, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ freq, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_freq, file = "models/MultiLevelModel_logCVR_hypertrophy_freq")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_freq, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_freq)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_freq <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_freq$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_freq$sigma2) + (MultiLevelModel_logCVR_hypertrophy_freq$k-MultiLevelModel_logCVR_hypertrophy_freq$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_freq <- 100 * MultiLevelModel_logCVR_hypertrophy_freq$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_freq$sigma2) + (MultiLevelModel_logCVR_hypertrophy_freq$k-MultiLevelModel_logCVR_hypertrophy_freq$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_freq <- robust(MultiLevelModel_logCVR_hypertrophy_freq, subset(Data_logCVR_freq, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_freq, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_freq")

### Number of exercises per workout
# SMD
# Strength
Data_SMD_strength_exercises <- Data_SMD_strength %>%
  filter(!is.na(exercises))

MultiLevelModel_strength_exercises <- rma.mv(yi, V=vi, data=Data_SMD_strength_exercises,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ exercises, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_exercises, file = "models/MultiLevelModel_strength_exercises")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_exercises$vi)
X <- model.matrix(MultiLevelModel_strength_exercises)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_exercises <- 100 * sum(MultiLevelModel_strength_exercises$sigma2) / (sum(MultiLevelModel_strength_exercises$sigma2) + (MultiLevelModel_strength_exercises$k-MultiLevelModel_strength_exercises$p)/sum(diag(P)))
I2bw_strength_exercises <- 100 * MultiLevelModel_strength_exercises$sigma2 / (sum(MultiLevelModel_strength_exercises$sigma2) + (MultiLevelModel_strength_exercises$k-MultiLevelModel_strength_exercises$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_exercises <- robust(MultiLevelModel_strength_exercises, Data_SMD_strength_exercises$study)

save(RobuEstMultiLevelModel_strength_exercises, file = "models/RobuEstMultiLevelModel_strength_exercises")

# Hypertrophy
Data_SMD_hypertrophy_exercises <- Data_SMD_hypertrophy %>%
  filter(!is.na(exercises))

MultiLevelModel_hypertrophy_exercises <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_exercises,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ exercises, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_exercises, file = "models/MultiLevelModel_hypertrophy_exercises")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_exercises$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_exercises)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_exercises <- 100 * sum(MultiLevelModel_hypertrophy_exercises$sigma2) / (sum(MultiLevelModel_hypertrophy_exercises$sigma2) + (MultiLevelModel_hypertrophy_exercises$k-MultiLevelModel_hypertrophy_exercises$p)/sum(diag(P)))
I2bw_hypertrophy_exercises <- 100 * MultiLevelModel_hypertrophy_exercises$sigma2 / (sum(MultiLevelModel_hypertrophy_exercises$sigma2) + (MultiLevelModel_hypertrophy_exercises$k-MultiLevelModel_hypertrophy_exercises$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_exercises <- robust(MultiLevelModel_hypertrophy_exercises, Data_SMD_hypertrophy_exercises$study)

save(RobuEstMultiLevelModel_hypertrophy_exercises, file = "models/RobuEstMultiLevelModel_hypertrophy_exercises")

# logCVR
Data_logCVR_exercises <- Data_logCVR %>%
  filter(!is.na(exercises))

# Strength
MultiLevelModel_logCVR_strength_exercises <- rma.mv(yi, V=vi, data=subset(Data_logCVR_exercises, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ exercises, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_exercises, file = "models/MultiLevelModel_logCVR_strength_exercises")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_exercises, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_exercises)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_exercises <- 100 * sum(MultiLevelModel_logCVR_strength_exercises$sigma2) / (sum(MultiLevelModel_logCVR_strength_exercises$sigma2) + (MultiLevelModel_logCVR_strength_exercises$k-MultiLevelModel_logCVR_strength_exercises$p)/sum(diag(P)))
I2bw_logCVR_strength_exercises <- 100 * MultiLevelModel_logCVR_strength_exercises$sigma2 / (sum(MultiLevelModel_logCVR_strength_exercises$sigma2) + (MultiLevelModel_logCVR_strength_exercises$k-MultiLevelModel_logCVR_strength_exercises$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_exercises <- robust(MultiLevelModel_logCVR_strength_exercises, subset(Data_logCVR_exercises, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_exercises, file = "models/RobuEstMultiLevelModel_logCVR_strength_exercises")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_exercises <- rma.mv(yi, V=vi, data=subset(Data_logCVR_exercises, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ exercises, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_exercises, file = "models/MultiLevelModel_logCVR_hypertrophy_exercises")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_exercises, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_exercises)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_exercises <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_exercises$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_exercises$sigma2) + (MultiLevelModel_logCVR_hypertrophy_exercises$k-MultiLevelModel_logCVR_hypertrophy_exercises$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_exercises <- 100 * MultiLevelModel_logCVR_hypertrophy_exercises$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_exercises$sigma2) + (MultiLevelModel_logCVR_hypertrophy_exercises$k-MultiLevelModel_logCVR_hypertrophy_exercises$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_exercises <- robust(MultiLevelModel_logCVR_hypertrophy_exercises, subset(Data_logCVR_exercises, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_exercises")

### Number of sets per exercise
# SMD
# Strength
Data_SMD_strength_sets_exercise <- Data_SMD_strength %>%
  filter(!is.na(sets_exercise))

MultiLevelModel_strength_sets_exercise <- rma.mv(yi, V=vi, data=Data_SMD_strength_sets_exercise,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ sets_exercise, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_sets_exercise, file = "models/MultiLevelModel_strength_sets_exercise")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_sets_exercise$vi)
X <- model.matrix(MultiLevelModel_strength_sets_exercise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_sets_exercise <- 100 * sum(MultiLevelModel_strength_sets_exercise$sigma2) / (sum(MultiLevelModel_strength_sets_exercise$sigma2) + (MultiLevelModel_strength_sets_exercise$k-MultiLevelModel_strength_sets_exercise$p)/sum(diag(P)))
I2bw_strength_sets_exercise <- 100 * MultiLevelModel_strength_sets_exercise$sigma2 / (sum(MultiLevelModel_strength_sets_exercise$sigma2) + (MultiLevelModel_strength_sets_exercise$k-MultiLevelModel_strength_sets_exercise$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_sets_exercise <- robust(MultiLevelModel_strength_sets_exercise, Data_SMD_strength_sets_exercise$study)

save(RobuEstMultiLevelModel_strength_sets_exercise, file = "models/RobuEstMultiLevelModel_strength_sets_exercise")

# Hypertrophy
Data_SMD_hypertrophy_sets_exercise <- Data_SMD_hypertrophy %>%
  filter(!is.na(sets_exercise))

MultiLevelModel_hypertrophy_sets_exercise <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_sets_exercise,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ sets_exercise, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_sets_exercise, file = "models/MultiLevelModel_hypertrophy_sets_exercise")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_sets_exercise$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_sets_exercise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_sets_exercise <- 100 * sum(MultiLevelModel_hypertrophy_sets_exercise$sigma2) / (sum(MultiLevelModel_hypertrophy_sets_exercise$sigma2) + (MultiLevelModel_hypertrophy_sets_exercise$k-MultiLevelModel_hypertrophy_sets_exercise$p)/sum(diag(P)))
I2bw_hypertrophy_sets_exercise <- 100 * MultiLevelModel_hypertrophy_sets_exercise$sigma2 / (sum(MultiLevelModel_hypertrophy_sets_exercise$sigma2) + (MultiLevelModel_hypertrophy_sets_exercise$k-MultiLevelModel_hypertrophy_sets_exercise$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_sets_exercise <- robust(MultiLevelModel_hypertrophy_sets_exercise, Data_SMD_hypertrophy_sets_exercise$study)

save(RobuEstMultiLevelModel_hypertrophy_sets_exercise, file = "models/RobuEstMultiLevelModel_hypertrophy_sets_exercise")

# logCVR
Data_logCVR_sets_exercise <- Data_logCVR %>%
  filter(!is.na(sets_exercise))

# Strength
MultiLevelModel_logCVR_strength_sets_exercise <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sets_exercise, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ sets_exercise, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_sets_exercise, file = "models/MultiLevelModel_logCVR_strength_sets_exercise")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_sets_exercise, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_sets_exercise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_sets_exercise <- 100 * sum(MultiLevelModel_logCVR_strength_sets_exercise$sigma2) / (sum(MultiLevelModel_logCVR_strength_sets_exercise$sigma2) + (MultiLevelModel_logCVR_strength_sets_exercise$k-MultiLevelModel_logCVR_strength_sets_exercise$p)/sum(diag(P)))
I2bw_logCVR_strength_sets_exercise <- 100 * MultiLevelModel_logCVR_strength_sets_exercise$sigma2 / (sum(MultiLevelModel_logCVR_strength_sets_exercise$sigma2) + (MultiLevelModel_logCVR_strength_sets_exercise$k-MultiLevelModel_logCVR_strength_sets_exercise$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_sets_exercise <- robust(MultiLevelModel_logCVR_strength_sets_exercise, subset(Data_logCVR_sets_exercise, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_sets_exercise, file = "models/RobuEstMultiLevelModel_logCVR_strength_sets_exercise")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_sets_exercise <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sets_exercise, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ sets_exercise, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_sets_exercise, file = "models/MultiLevelModel_logCVR_hypertrophy_sets_exercise")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_sets_exercise, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_sets_exercise)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_sets_exercise <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_sets_exercise$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_sets_exercise$sigma2) + (MultiLevelModel_logCVR_hypertrophy_sets_exercise$k-MultiLevelModel_logCVR_hypertrophy_sets_exercise$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_sets_exercise <- 100 * MultiLevelModel_logCVR_hypertrophy_sets_exercise$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_sets_exercise$sigma2) + (MultiLevelModel_logCVR_hypertrophy_sets_exercise$k-MultiLevelModel_logCVR_hypertrophy_sets_exercise$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise <- robust(MultiLevelModel_logCVR_hypertrophy_sets_exercise, subset(Data_logCVR_sets_exercise, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise")

### Repetitions per exercise
# SMD
# Strength
Data_SMD_strength_reps <- Data_SMD_strength %>%
  filter(!is.na(reps))

MultiLevelModel_strength_reps <- rma.mv(yi, V=vi, data=Data_SMD_strength_reps,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ reps, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_reps, file = "models/MultiLevelModel_strength_reps")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_reps$vi)
X <- model.matrix(MultiLevelModel_strength_reps)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_reps <- 100 * sum(MultiLevelModel_strength_reps$sigma2) / (sum(MultiLevelModel_strength_reps$sigma2) + (MultiLevelModel_strength_reps$k-MultiLevelModel_strength_reps$p)/sum(diag(P)))
I2bw_strength_reps <- 100 * MultiLevelModel_strength_reps$sigma2 / (sum(MultiLevelModel_strength_reps$sigma2) + (MultiLevelModel_strength_reps$k-MultiLevelModel_strength_reps$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_reps <- robust(MultiLevelModel_strength_reps, Data_SMD_strength_reps$study)

save(RobuEstMultiLevelModel_strength_reps, file = "models/RobuEstMultiLevelModel_strength_reps")

# Hypertrophy
Data_SMD_hypertrophy_reps <- Data_SMD_hypertrophy %>%
  filter(!is.na(reps))

MultiLevelModel_hypertrophy_reps <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_reps,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ reps, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_reps, file = "models/MultiLevelModel_hypertrophy_reps")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_reps$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_reps)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_reps <- 100 * sum(MultiLevelModel_hypertrophy_reps$sigma2) / (sum(MultiLevelModel_hypertrophy_reps$sigma2) + (MultiLevelModel_hypertrophy_reps$k-MultiLevelModel_hypertrophy_reps$p)/sum(diag(P)))
I2bw_hypertrophy_reps <- 100 * MultiLevelModel_hypertrophy_reps$sigma2 / (sum(MultiLevelModel_hypertrophy_reps$sigma2) + (MultiLevelModel_hypertrophy_reps$k-MultiLevelModel_hypertrophy_reps$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_reps <- robust(MultiLevelModel_hypertrophy_reps, Data_SMD_hypertrophy_reps$study)

save(RobuEstMultiLevelModel_hypertrophy_reps, file = "models/RobuEstMultiLevelModel_hypertrophy_reps")

# logCVR
Data_logCVR_reps <- Data_logCVR %>%
  filter(!is.na(reps))

# Strength
MultiLevelModel_logCVR_strength_reps <- rma.mv(yi, V=vi, data=subset(Data_logCVR_reps, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ reps, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_reps, file = "models/MultiLevelModel_logCVR_strength_reps")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_reps, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_reps)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_reps <- 100 * sum(MultiLevelModel_logCVR_strength_reps$sigma2) / (sum(MultiLevelModel_logCVR_strength_reps$sigma2) + (MultiLevelModel_logCVR_strength_reps$k-MultiLevelModel_logCVR_strength_reps$p)/sum(diag(P)))
I2bw_logCVR_strength_reps <- 100 * MultiLevelModel_logCVR_strength_reps$sigma2 / (sum(MultiLevelModel_logCVR_strength_reps$sigma2) + (MultiLevelModel_logCVR_strength_reps$k-MultiLevelModel_logCVR_strength_reps$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_reps <- robust(MultiLevelModel_logCVR_strength_reps, subset(Data_logCVR_reps, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_reps, file = "models/RobuEstMultiLevelModel_logCVR_strength_reps")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_reps <- rma.mv(yi, V=vi, data=subset(Data_logCVR_reps, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ reps, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_reps, file = "models/MultiLevelModel_logCVR_hypertrophy_reps")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_reps, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_reps)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_reps <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_reps$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_reps$sigma2) + (MultiLevelModel_logCVR_hypertrophy_reps$k-MultiLevelModel_logCVR_hypertrophy_reps$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_reps <- 100 * MultiLevelModel_logCVR_hypertrophy_reps$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_reps$sigma2) + (MultiLevelModel_logCVR_hypertrophy_reps$k-MultiLevelModel_logCVR_hypertrophy_reps$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_reps <- robust(MultiLevelModel_logCVR_hypertrophy_reps, subset(Data_logCVR_reps, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_reps, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_reps")

### Relative load (%1RM)
# SMD
# Strength
Data_SMD_strength_load <- Data_SMD_strength %>%
  filter(!is.na(load))

MultiLevelModel_strength_load <- rma.mv(yi, V=vi, data=Data_SMD_strength_load,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ load, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_load, file = "models/MultiLevelModel_strength_load")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_load$vi)
X <- model.matrix(MultiLevelModel_strength_load)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_load <- 100 * sum(MultiLevelModel_strength_load$sigma2) / (sum(MultiLevelModel_strength_load$sigma2) + (MultiLevelModel_strength_load$k-MultiLevelModel_strength_load$p)/sum(diag(P)))
I2bw_strength_load <- 100 * MultiLevelModel_strength_load$sigma2 / (sum(MultiLevelModel_strength_load$sigma2) + (MultiLevelModel_strength_load$k-MultiLevelModel_strength_load$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_load <- robust(MultiLevelModel_strength_load, Data_SMD_strength_load$study)

save(RobuEstMultiLevelModel_strength_load, file = "models/RobuEstMultiLevelModel_strength_load")

# Hypertrophy
Data_SMD_hypertrophy_load <- Data_SMD_hypertrophy %>%
  filter(!is.na(load))

MultiLevelModel_hypertrophy_load <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_load,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ load, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_load, file = "models/MultiLevelModel_hypertrophy_load")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_load$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_load)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_load <- 100 * sum(MultiLevelModel_hypertrophy_load$sigma2) / (sum(MultiLevelModel_hypertrophy_load$sigma2) + (MultiLevelModel_hypertrophy_load$k-MultiLevelModel_hypertrophy_load$p)/sum(diag(P)))
I2bw_hypertrophy_load <- 100 * MultiLevelModel_hypertrophy_load$sigma2 / (sum(MultiLevelModel_hypertrophy_load$sigma2) + (MultiLevelModel_hypertrophy_load$k-MultiLevelModel_hypertrophy_load$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_load <- robust(MultiLevelModel_hypertrophy_load, Data_SMD_hypertrophy_load$study)

save(RobuEstMultiLevelModel_hypertrophy_load, file = "models/RobuEstMultiLevelModel_hypertrophy_load")

# logCVR
Data_logCVR_load <- Data_logCVR %>%
  filter(!is.na(load))

# Strength
MultiLevelModel_logCVR_strength_load <- rma.mv(yi, V=vi, data=subset(Data_logCVR_load, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ load, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_load, file = "models/MultiLevelModel_logCVR_strength_load")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_load, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_load)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_load <- 100 * sum(MultiLevelModel_logCVR_strength_load$sigma2) / (sum(MultiLevelModel_logCVR_strength_load$sigma2) + (MultiLevelModel_logCVR_strength_load$k-MultiLevelModel_logCVR_strength_load$p)/sum(diag(P)))
I2bw_logCVR_strength_load <- 100 * MultiLevelModel_logCVR_strength_load$sigma2 / (sum(MultiLevelModel_logCVR_strength_load$sigma2) + (MultiLevelModel_logCVR_strength_load$k-MultiLevelModel_logCVR_strength_load$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_load <- robust(MultiLevelModel_logCVR_strength_load, subset(Data_logCVR_load, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_load, file = "models/RobuEstMultiLevelModel_logCVR_strength_load")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_load <- rma.mv(yi, V=vi, data=subset(Data_logCVR_load, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ load, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_load, file = "models/MultiLevelModel_logCVR_hypertrophy_load")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_load, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_load)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_load <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_load$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_load$sigma2) + (MultiLevelModel_logCVR_hypertrophy_load$k-MultiLevelModel_logCVR_hypertrophy_load$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_load <- 100 * MultiLevelModel_logCVR_hypertrophy_load$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_load$sigma2) + (MultiLevelModel_logCVR_hypertrophy_load$k-MultiLevelModel_logCVR_hypertrophy_load$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_load <- robust(MultiLevelModel_logCVR_hypertrophy_load, subset(Data_logCVR_load, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_load, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_load")

### Trained to task failure?
# SMD
# Strength
Data_SMD_strength_task_failure_y_n <- Data_SMD_strength %>%
  filter(!is.na(task_failure_y_n))

MultiLevelModel_strength_task_failure_y_n <- rma.mv(yi, V=vi, data=Data_SMD_strength_task_failure_y_n,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ task_failure_y_n - 1, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_task_failure_y_n, file = "models/MultiLevelModel_strength_task_failure_y_n")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_task_failure_y_n$vi)
X <- model.matrix(MultiLevelModel_strength_task_failure_y_n)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_task_failure_y_n <- 100 * sum(MultiLevelModel_strength_task_failure_y_n$sigma2) / (sum(MultiLevelModel_strength_task_failure_y_n$sigma2) + (MultiLevelModel_strength_task_failure_y_n$k-MultiLevelModel_strength_task_failure_y_n$p)/sum(diag(P)))
I2bw_strength_task_failure_y_n <- 100 * MultiLevelModel_strength_task_failure_y_n$sigma2 / (sum(MultiLevelModel_strength_task_failure_y_n$sigma2) + (MultiLevelModel_strength_task_failure_y_n$k-MultiLevelModel_strength_task_failure_y_n$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_task_failure_y_n <- robust(MultiLevelModel_strength_task_failure_y_n, Data_SMD_strength_task_failure_y_n$study)

save(RobuEstMultiLevelModel_strength_task_failure_y_n, file = "models/RobuEstMultiLevelModel_strength_task_failure_y_n")

# Hypertrophy
Data_SMD_hypertrophy_task_failure_y_n <- Data_SMD_hypertrophy %>%
  filter(!is.na(task_failure_y_n))

MultiLevelModel_hypertrophy_task_failure_y_n <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_task_failure_y_n,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ task_failure_y_n - 1, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_task_failure_y_n, file = "models/MultiLevelModel_hypertrophy_task_failure_y_n")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_task_failure_y_n$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_task_failure_y_n)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_task_failure_y_n <- 100 * sum(MultiLevelModel_hypertrophy_task_failure_y_n$sigma2) / (sum(MultiLevelModel_hypertrophy_task_failure_y_n$sigma2) + (MultiLevelModel_hypertrophy_task_failure_y_n$k-MultiLevelModel_hypertrophy_task_failure_y_n$p)/sum(diag(P)))
I2bw_hypertrophy_task_failure_y_n <- 100 * MultiLevelModel_hypertrophy_task_failure_y_n$sigma2 / (sum(MultiLevelModel_hypertrophy_task_failure_y_n$sigma2) + (MultiLevelModel_hypertrophy_task_failure_y_n$k-MultiLevelModel_hypertrophy_task_failure_y_n$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_task_failure_y_n <- robust(MultiLevelModel_hypertrophy_task_failure_y_n, Data_SMD_hypertrophy_task_failure_y_n$study)

save(RobuEstMultiLevelModel_hypertrophy_task_failure_y_n, file = "models/RobuEstMultiLevelModel_hypertrophy_task_failure_y_n")

# logCVR
Data_logCVR_task_failure_y_n <- Data_logCVR %>%
  filter(!is.na(task_failure_y_n))

# Strength
MultiLevelModel_logCVR_strength_task_failure_y_n <- rma.mv(yi, V=vi, data=subset(Data_logCVR_task_failure_y_n, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ task_failure_y_n - 1, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_task_failure_y_n, file = "models/MultiLevelModel_logCVR_strength_task_failure_y_n")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_task_failure_y_n, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_task_failure_y_n)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_task_failure_y_n <- 100 * sum(MultiLevelModel_logCVR_strength_task_failure_y_n$sigma2) / (sum(MultiLevelModel_logCVR_strength_task_failure_y_n$sigma2) + (MultiLevelModel_logCVR_strength_task_failure_y_n$k-MultiLevelModel_logCVR_strength_task_failure_y_n$p)/sum(diag(P)))
I2bw_logCVR_strength_task_failure_y_n <- 100 * MultiLevelModel_logCVR_strength_task_failure_y_n$sigma2 / (sum(MultiLevelModel_logCVR_strength_task_failure_y_n$sigma2) + (MultiLevelModel_logCVR_strength_task_failure_y_n$k-MultiLevelModel_logCVR_strength_task_failure_y_n$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n <- robust(MultiLevelModel_logCVR_strength_task_failure_y_n, subset(Data_logCVR_task_failure_y_n, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n, file = "models/RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_task_failure_y_n <- rma.mv(yi, V=vi, data=subset(Data_logCVR_task_failure_y_n, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ task_failure_y_n - 1, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n, file = "models/MultiLevelModel_logCVR_hypertrophy_task_failure_y_n")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_task_failure_y_n, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_task_failure_y_n <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$sigma2) + (MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$k-MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_task_failure_y_n <- 100 * MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$sigma2) + (MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$k-MultiLevelModel_logCVR_hypertrophy_task_failure_y_n$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n <- robust(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n, subset(Data_logCVR_task_failure_y_n, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n")

### Outcome measurement type
# SMD
# Strength
Data_SMD_strength_measure <- Data_SMD_strength %>%
  filter(!is.na(measure))

MultiLevelModel_strength_measure <- rma.mv(yi, V=vi, data=Data_SMD_strength_measure,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | group), mods = ~ measure - 1, method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_strength_measure, file = "models/MultiLevelModel_strength_measure")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_strength_measure$vi)
X <- model.matrix(MultiLevelModel_strength_measure)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_measure <- 100 * sum(MultiLevelModel_strength_measure$sigma2) / (sum(MultiLevelModel_strength_measure$sigma2) + (MultiLevelModel_strength_measure$k-MultiLevelModel_strength_measure$p)/sum(diag(P)))
I2bw_strength_measure <- 100 * MultiLevelModel_strength_measure$sigma2 / (sum(MultiLevelModel_strength_measure$sigma2) + (MultiLevelModel_strength_measure$k-MultiLevelModel_strength_measure$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_measure <- robust(MultiLevelModel_strength_measure, Data_SMD_strength_measure$study)

save(RobuEstMultiLevelModel_strength_measure, file = "models/RobuEstMultiLevelModel_strength_measure")

# Hypertrophy
Data_SMD_hypertrophy_measure <- Data_SMD_hypertrophy %>%
  filter(!is.na(measure))

MultiLevelModel_hypertrophy_measure <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_measure,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | group), mods = ~ measure - 1, method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_hypertrophy_measure, file = "models/MultiLevelModel_hypertrophy_measure")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_SMD_hypertrophy_measure$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_measure)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_measure <- 100 * sum(MultiLevelModel_hypertrophy_measure$sigma2) / (sum(MultiLevelModel_hypertrophy_measure$sigma2) + (MultiLevelModel_hypertrophy_measure$k-MultiLevelModel_hypertrophy_measure$p)/sum(diag(P)))
I2bw_hypertrophy_measure <- 100 * MultiLevelModel_hypertrophy_measure$sigma2 / (sum(MultiLevelModel_hypertrophy_measure$sigma2) + (MultiLevelModel_hypertrophy_measure$k-MultiLevelModel_hypertrophy_measure$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_measure <- robust(MultiLevelModel_hypertrophy_measure, Data_SMD_hypertrophy_measure$study)

save(RobuEstMultiLevelModel_hypertrophy_measure, file = "models/RobuEstMultiLevelModel_hypertrophy_measure")

# logCVR
Data_logCVR_measure <- Data_logCVR %>%
  filter(!is.na(measure))

# Strength
MultiLevelModel_logCVR_strength_measure <- rma.mv(yi, V=vi, data=subset(Data_logCVR_measure, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | group), mods = ~ measure - 1, method="REML",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_measure, file = "models/MultiLevelModel_logCVR_strength_measure")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_measure, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength_measure)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength_measure <- 100 * sum(MultiLevelModel_logCVR_strength_measure$sigma2) / (sum(MultiLevelModel_logCVR_strength_measure$sigma2) + (MultiLevelModel_logCVR_strength_measure$k-MultiLevelModel_logCVR_strength_measure$p)/sum(diag(P)))
I2bw_logCVR_strength_measure <- 100 * MultiLevelModel_logCVR_strength_measure$sigma2 / (sum(MultiLevelModel_logCVR_strength_measure$sigma2) + (MultiLevelModel_logCVR_strength_measure$k-MultiLevelModel_logCVR_strength_measure$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_measure <- robust(MultiLevelModel_logCVR_strength_measure, subset(Data_logCVR_measure, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_measure, file = "models/RobuEstMultiLevelModel_logCVR_strength_measure")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_measure <- rma.mv(yi, V=vi, data=subset(Data_logCVR_measure, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | group), mods = ~ measure - 1, method="REML",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_measure, file = "models/MultiLevelModel_logCVR_hypertrophy_measure")

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR_measure, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy_measure)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy_measure <- 100 * sum(MultiLevelModel_logCVR_hypertrophy_measure$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy_measure$sigma2) + (MultiLevelModel_logCVR_hypertrophy_measure$k-MultiLevelModel_logCVR_hypertrophy_measure$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy_measure <- 100 * MultiLevelModel_logCVR_hypertrophy_measure$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy_measure$sigma2) + (MultiLevelModel_logCVR_hypertrophy_measure$k-MultiLevelModel_logCVR_hypertrophy_measure$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_measure <- robust(MultiLevelModel_logCVR_hypertrophy_measure, subset(Data_logCVR_measure, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_measure, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_measure")

### Collate all moderator estimates for log CVR and SMD

mods_SMD_logCVR <- rbind(data.frame(Outcome = "strength",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_strength$ci.ub, 2),
                                    `I2 study` = round(I2bw_strength, 2)[1],
                                    `I2 group` = round(I2bw_strength, 2)[2]),
                         data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "TESTEX score",
                                      Estimate = round(RobuEstMultiLevelModel_strength_TESTEX$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_TESTEX$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_TESTEX$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_strength_TESTEX, 2)[1],
                                    `I2 group` = round(I2bw_strength_TESTEX, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Age",
                                      Estimate = round(RobuEstMultiLevelModel_strength_age$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_age$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_age$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_age, 2)[1],
                                      `I2 group` = round(I2bw_strength_age, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Proportion Male",
                                      Estimate = round(RobuEstMultiLevelModel_strength_sex_._male$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_sex_._male$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_sex_._male$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_sex_._male, 2)[1],
                                      `I2 group` = round(I2bw_strength_sex_._male, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Weight",
                                      Estimate = round(RobuEstMultiLevelModel_strength_weight$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_weight$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_weight$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_weight, 2)[1],
                                      `I2 group` = round(I2bw_strength_weight, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "BMI",
                                      Estimate = round(RobuEstMultiLevelModel_strength_bmi$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_bmi$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_bmi$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_bmi, 2)[1],
                                      `I2 group` = round(I2bw_strength_bmi, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                      Estimate = round(RobuEstMultiLevelModel_strength_train_status$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_strength_train_status$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_strength_train_status$ci.ub, 2),
                                      `I2 study` = round(I2bw_strength_train_status, 2)[1],
                                      `I2 group` = round(I2bw_strength_train_status, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Healthy Sample", "Clinical Sample"),
                                      Estimate = round(RobuEstMultiLevelModel_strength_healthy_clinical$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_strength_healthy_clinical$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_strength_healthy_clinical$ci.ub, 2),
                                      `I2 study` = round(I2bw_strength_healthy_clinical, 2)[1],
                                      `I2 group` = round(I2bw_strength_healthy_clinical, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                      Estimate = round(RobuEstMultiLevelModel_strength_RT_only$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_strength_RT_only$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_strength_RT_only$ci.ub, 2),
                                      `I2 study` = round(I2bw_strength_RT_only, 2)[1],
                                      `I2 group` = round(I2bw_strength_RT_only, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Duration (weeks)",
                                      Estimate = round(RobuEstMultiLevelModel_strength_weeks$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_weeks$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_weeks$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_weeks, 2)[1],
                                      `I2 group` = round(I2bw_strength_weeks, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Weekly Frequency",
                                      Estimate = round(RobuEstMultiLevelModel_strength_freq$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_freq$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_freq$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_freq, 2)[1],
                                      `I2 group` = round(I2bw_strength_freq, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Number of Exercises",
                                      Estimate = round(RobuEstMultiLevelModel_strength_exercises$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_exercises$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_exercises$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_exercises, 2)[1],
                                      `I2 group` = round(I2bw_strength_exercises, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Sets per Exercise",
                                      Estimate = round(RobuEstMultiLevelModel_strength_sets_exercise$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_sets_exercise$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_sets_exercise$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_sets_exercise, 2)[1],
                                      `I2 group` = round(I2bw_strength_sets_exercise, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Number of Repetitions",
                                      Estimate = round(RobuEstMultiLevelModel_strength_reps$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_reps$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_reps$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_reps, 2)[1],
                                      `I2 group` = round(I2bw_strength_reps, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Load (%1RM)",
                                      Estimate = round(RobuEstMultiLevelModel_strength_load$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_strength_load$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_strength_load$ci.ub, 2)[2],
                                      `I2 study` = round(I2bw_strength_load, 2)[1],
                                      `I2 group` = round(I2bw_strength_load, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                      Estimate = round(RobuEstMultiLevelModel_strength_task_failure_y_n$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_strength_task_failure_y_n$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_strength_task_failure_y_n$ci.ub, 2),
                                      `I2 study` = round(I2bw_strength_task_failure_y_n, 2)[1],
                                      `I2 group` = round(I2bw_strength_task_failure_y_n, 2)[2]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Outcome Measure (12RM)","Outcome Measure (1RM)","Outcome Measure (3RM)","Outcome Measure (5RM)","Outcome Measure (6RM)","Outcome Measure (Isokinetic)", "Outcome Measure (Isometric)"),
                                      Estimate = round(RobuEstMultiLevelModel_strength_measure$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_strength_measure$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_strength_measure$ci.ub, 2),
                                      `I2 study` = round(I2bw_strength_measure, 2)[1],
                                      `I2 group` = round(I2bw_strength_measure, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_TESTEX, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_TESTEX, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_age$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_age, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_age, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_sex_._male, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_sex_._male, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_weight$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_weight, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_weight, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_bmi$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_bmi, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_bmi, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_train_status$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength_train_status, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_train_status, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength_healthy_clinical, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_healthy_clinical, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength_RT_only, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_RT_only, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_weeks$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_weeks, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_weeks, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_freq$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_freq, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_freq, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_exercises$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_exercises, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_exercises, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_sets_exercise, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_sets_exercise, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_reps$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_reps, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_reps, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_load$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_strength_load, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_load, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength_task_failure_y_n, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_task_failure_y_n, 2)[2]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Outcome Measure (12RM)","Outcome Measure (1RM)","Outcome Measure (3RM)","Outcome Measure (5RM)","Outcome Measure (6RM)","Outcome Measure (Isokinetic)", "Outcome Measure (Isometric)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_measure$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_strength_measure, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_strength_measure, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_TESTEX$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_TESTEX, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_TESTEX, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_age$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_age, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_age, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_sex_._male$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_sex_._male, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_sex_._male, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_weight$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_weight, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_weight, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_bmi$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_bmi, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_bmi, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_train_status$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy_train_status, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_train_status, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_healthy_clinical$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy_healthy_clinical, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_healthy_clinical, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_RT_only$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy_RT_only, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_RT_only, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_weeks$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_weeks, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_weeks, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_freq$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_freq, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_freq, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_exercises$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_exercises, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_exercises, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_sets_exercise$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_sets_exercise, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_sets_exercise, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_reps$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_reps, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_reps, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_load$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_hypertrophy_load, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_load, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_task_failure_y_n$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy_task_failure_y_n, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_task_failure_y_n, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Outcome Measure (BIA)","Outcome Measure (Biopsy: Type i)","Outcome Measure (Biopsy: Type ii)","Outcome Measure (Biopsy: Type iia)","Outcome Measure (Biopsy: Type iib)","Outcome Measure (BodPod)", "Outcome Measure (Circumference)","Outcome Measure (CT)","Outcome Measure (DXA)","Outcome Measure (Hydrostatic Weighing)","Outcome Measure (MRI)","Outcome Measure (Skinfold)","Outcome Measure (US)"),
                                    Estimate = round(RobuEstMultiLevelModel_hypertrophy_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_hypertrophy_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_hypertrophy_measure$ci.ub, 2),
                                    `I2 study` = round(I2bw_hypertrophy_measure, 2)[1],
                                    `I2 group` = round(I2bw_hypertrophy_measure, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_TESTEX, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_TESTEX, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_age, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_age, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_sex_._male, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_sex_._male, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_weight, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_weight, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_bmi, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_bmi, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_train_status, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_train_status, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_healthy_clinical, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_healthy_clinical, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_RT_only, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_RT_only, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_weeks, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_weeks, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_freq, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_freq, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_exercises, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_exercises, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_sets_exercise, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_sets_exercise, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_reps, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_reps, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$ci.ub, 2)[2],
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_load, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_load, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_task_failure_y_n, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_task_failure_y_n, 2)[2]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Outcome Measure (BIA)","Outcome Measure (Biopsy: Type i)","Outcome Measure (Biopsy: Type ii)","Outcome Measure (Biopsy: Type iia)","Outcome Measure (Biopsy: Type iib)","Outcome Measure (BodPod)", "Outcome Measure (Circumference)","Outcome Measure (CT)","Outcome Measure (DXA)","Outcome Measure (Hydrostatic Weighing)","Outcome Measure (MRI)","Outcome Measure (Skinfold)","Outcome Measure (US)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$ci.ub, 2),
                                    `I2 study` = round(I2bw_logCVR_hypertrophy_measure, 2)[1],
                                    `I2 group` = round(I2bw_logCVR_hypertrophy_measure, 2)[2])
                           ) 

rownames(mods_SMD_logCVR) <- NULL

save(mods_SMD_logCVR, file = "models/mods_SMD_logCVR")

library(kableExtra)               
library(webshot)

mods_SMD_logCVR <- mods_SMD_logCVR[,-c(1,2)]

knitr::kable(
  mods_SMD_logCVR,
  caption = "Meta-regression models for study and participant characteristics as moderators of SMD and lnCVR effect sizes",
  align = c("l","c")
) %>%
  pack_rows("Strength Outcomes",1, 54, bold = FALSE) %>%
  pack_rows("Hypertrophy Outcomes", 55, 120, bold = FALSE) %>%
  pack_rows("SMD", 1, 27, bold = FALSE) %>%
  pack_rows("lnCVR", 28, 54, bold = FALSE) %>%
  pack_rows("SMD", 55, 87, bold = FALSE) %>%
  pack_rows("lnCVR", 88, 120, bold = FALSE) %>%
  footnote(general = c("SMD = standardised mean difference; ", "lnCVR = log ratio of coefficient of variation; ", "RT = resistance training")
  ) %>%
  row_spec(0, bold = TRUE) %>%
  kable_classic(full_width = FALSE) %>%
  kable_styling() %>%
  save_kable("models/table_mods_SMD_logCVR.html")

webshot("models/table_mods_SMD_logCVR.html", "models/table_mods_SMD_logCVR.pdf")

#### Galaxy plot

Data_SMD_logCVR <- Data_SMD %>%
  select(study, group, es, yi, vi) %>%
  left_join(Data_logCVR, by = c("study", "group", "es")) %>%
  filter(!is.na(yi.x) & !is.na(yi.y)) 

SMD_v_logCVR_plot <- ggplot(data = Data_SMD_logCVR, aes(x=yi.x, y=yi.y, color = outcome)) +
  scale_fill_manual("Outcome", values = alpha(c("#009E73", "#D55E00"),0.5)) +
  scale_color_manual("Outcome", values = alpha(c("#009E73", "#D55E00"),0.5)) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point() +
  labs(x = "Standardised Mean Difference (Positive Values Favour Resistance Training)",
       y = "Log Coefficient of Variation Ratio") +
  scale_y_continuous(limits = c(-8,5), breaks = c(-7.5,-5,-2.5,0,2.5,5)) +
  scale_x_continuous(limits = c(-2,7.5), breaks = c(-1,0,1,2,3,4,5,6)) +
  theme_classic()

save(SMD_v_logCVR_plot, file = "plots/SMD_v_logCVR_plot")

SMD_v_logCVR_plot

ggsave("plots/SMD_v_logCVR_plot.tiff", width = 7.5, height = 5, device = "tiff", dpi = 300)




load(here("models","RobuEstMultiLevelModel_logCVR_strength"))
RobuEstMultiLevelModel_logCVR_strength$beta
