### This script runs all the analyses for the paper

# title: "Meta-Analysis of Variation in Sport and Exercise Science"
# subtitle: "Examples of Application Within Resistance Training Research"
# authors: 
#   - James Steele
# - James Fisher
# - Dave Smith
# - Brad Schoenfeld
# - Yefeng Yang
# - Shinichi Nakagawa

# Open packages
library(metafor)
library(tidyverse)
library(orchaRd)
library(bayestestR)
library(patchwork)
library(europepmc)
library(kableExtra)               
library(webshot)
library(ggtext)

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
                     random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                     control=list(optimizer="optim", optmethod="Nelder-Mead"))

RobuEstMeta_RT_ri <- robust(Meta_RT_ri, Data$study)

z2r_RT <- psych::fisherz2r(RobuEstMeta_RT_ri$b[1])

Data$RT_ri <- ifelse(is.na(Data$RT_ri), z2r_RT, Data$RT_ri)

Data <- escalc(measure = "ZCOR", ri = CON_ri, ni = CON_n, data = Data)

### Note, data is coded with study and arm as having explicit nesting so all random effects are (~ 1 | study, ~ 1 | arm)
Meta_CON_ri <- rma.mv(yi, V=vi, data=Data,
                     slab=paste(label),
                     random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
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
# Following Morris (2008) g_ppc2; see https://www.metafor-project.org/doku.php/analyses:morris2008
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

MultiLevelModel_SMD_strength <- rma.mv(yi, V=vi, data=Data_SMD_strength,
                                         slab=paste(label),
                                         random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength, file = "models/MultiLevelModel_SMD_strength")

### Calculate I^2 
I2_SMD_strength <- i2_ml(MultiLevelModel_SMD_strength)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength <- robust(MultiLevelModel_SMD_strength, Data_SMD_strength$study)

save(RobuEstMultiLevelModel_SMD_strength, file = "models/RobuEstMultiLevelModel_SMD_strength")

### Caterpillar plot 

# Overall estimate
diamond_SMD_strength <- data.frame(x = c(RobuEstMultiLevelModel_SMD_strength$b[1] + (RobuEstMultiLevelModel_SMD_strength$se*1.96),
                            RobuEstMultiLevelModel_SMD_strength$b[1],
                            RobuEstMultiLevelModel_SMD_strength$b[1] - (RobuEstMultiLevelModel_SMD_strength$se*1.96),
                            RobuEstMultiLevelModel_SMD_strength$b[1]),
                      y = c(-15,-25,-15,-5))

# Prediction interval
PI_SMD_strength <- as.data.frame(predict(RobuEstMultiLevelModel_SMD_strength))

# I^2 labels
I2_SMD_strength_lab <- data.frame(level = c("study", "arm", "es"),
                              I2 = I2_SMD_strength[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_SMD_strength <- Data_SMD_strength %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-2,7.5), breaks = c(-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_SMD_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_SMD_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 5, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_SMD_strength_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 5, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_SMD_strength, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_SMD_strength, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standardised Mean Difference (Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        )

save(forest_SMD_strength, file = "plots/forest_SMD_strength")

### Hypertrophy
Data_SMD_hypertrophy <- Data_SMD %>% 
    filter(!is.na(yi) &  outcome == "hypertrophy")

MultiLevelModel_SMD_hypertrophy <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy,
                                   slab=paste(label),
                                   random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy, file = "models/MultiLevelModel_SMD_hypertrophy")

### Calculate I^2 
I2_SMD_hypertrophy <- i2_ml(MultiLevelModel_SMD_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy <- robust(MultiLevelModel_SMD_hypertrophy, Data_SMD_hypertrophy$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_SMD_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_SMD_hypertrophy$b[1] + (RobuEstMultiLevelModel_SMD_hypertrophy$se*1.96),
                            RobuEstMultiLevelModel_SMD_hypertrophy$b[1],
                            RobuEstMultiLevelModel_SMD_hypertrophy$b[1] - (RobuEstMultiLevelModel_SMD_hypertrophy$se*1.96),
                            RobuEstMultiLevelModel_SMD_hypertrophy$b[1]),
                      y = c(-15,-25,-15,-5))

# Prediction interval
PI_SMD_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_SMD_hypertrophy))

# I^2 labels
I2_SMD_hypertrophy_lab <- data.frame(level = c("study", "arm", "es"),
                              I2 = I2_SMD_hypertrophy[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_SMD_hypertrophy <- Data_SMD_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-2,7.5), breaks = c(-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_SMD_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_SMD_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 5, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_SMD_hypertrophy_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 5, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_SMD_hypertrophy, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_SMD_hypertrophy, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standardised Mean Difference (Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_SMD_hypertrophy, file = "plots/forest_SMD_hypertrophy")

### Combine Plots
forest_SMD_plots <- (forest_SMD_strength / forest_SMD_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_SMD_plots, file = "plots/forest_SMD_plots")

forest_SMD_plots

ggsave("plots/forest_SMD_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)


###### Comparison of response ratios between RT and CON ######

### Response ratio difference effect size calculations 
# Following Lajuenesse (2011; 2015) for the interaction effect of a group * time factorial design

Data$lnRR_yi <- log(Data$RT_post_m/Data$RT_pre_m) - log(Data$CON_post_m/Data$CON_pre_m)

Data$lnRR_vi <- (Data$RT_post_sd^2/(Data$RT_post_m^2*Data$RT_n)) + (Data$RT_pre_sd^2/(Data$RT_pre_m^2*Data$RT_n)) + (Data$CON_post_sd^2/(Data$CON_post_m^2*Data$CON_n)) + (Data$CON_pre_sd^2/(Data$CON_pre_m^2*Data$CON_n))

Data_RR <- Data

### Strength
Data_RR_strength <- Data_RR %>% 
  filter(!is.na(lnRR_yi) & !is.na(lnRR_vi) & outcome == "strength")

MultiLevelModel_RR_strength <- rma.mv(lnRR_yi, V=lnRR_vi, data=Data_RR_strength,
                                   slab=paste(label),
                                   random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_RR_strength, file = "models/MultiLevelModel_RR_strength")

### Calculate I^2 
I2_RR_strength <- i2_ml(MultiLevelModel_RR_strength)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_RR_strength <- robust(MultiLevelModel_RR_strength, Data_RR_strength$study)

save(RobuEstMultiLevelModel_RR_strength, file = "models/RobuEstMultiLevelModel_RR_strength")

### Caterpillar plot 

# Overall estimate
diamond_RR_strength <- data.frame(x = c(RobuEstMultiLevelModel_RR_strength$b[1] + (RobuEstMultiLevelModel_RR_strength$se*1.96),
                                     RobuEstMultiLevelModel_RR_strength$b[1],
                                     RobuEstMultiLevelModel_RR_strength$b[1] - (RobuEstMultiLevelModel_RR_strength$se*1.96),
                                     RobuEstMultiLevelModel_RR_strength$b[1]),
                               y = c(-15,-25,-15,-5))

# Prediction interval
PI_RR_strength <- as.data.frame(predict(RobuEstMultiLevelModel_RR_strength))

# I^2 labels
I2_RR_strength_lab <- data.frame(level = c("study", "arm", "es"),
                              I2 = I2_RR_strength[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_RR_strength <- Data_RR_strength %>% 
  mutate(es = factor(es, levels = es[order(lnRR_yi)]),
         se = sqrt(lnRR_vi)) %>%
  ggplot(aes(x=(exp(lnRR_yi)-1)*100, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  # scale_x_continuous(limits = c(-20,40), breaks = c(-20,-10,0,10,20,30,40)) +
  geom_linerange(aes(xmin = (exp(lnRR_yi)-1)*100 - (se*1.96), xmax = (exp(lnRR_yi)-1)*100 + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_RR_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round((exp(pred)-1)*100,2)} [95% Confidence Interval: {round((exp(ci.lb)-1)*100,2)} to {round((exp(ci.ub)-1)*100,2)}]"),
                x = 100, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_RR_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round((exp(pi.lb)-1)*100,2)} to {round((exp(pi.ub)-1)*100,2)}]"),
                x = 100, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_RR_strength_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 100, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_RR_strength, aes(y=-15, yend=-15, x=(exp(pi.lb)-1)*100, xend=(exp(pi.ub)-1)*100),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_RR_strength, aes(x=(exp(x)-1)*100,y=y)) +
  labs(y = "",
       x = "Exponentiated Response Ratio (%; Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_RR_strength, file = "plots/forest_RR_strength")

### Hypertrophy
Data_RR_hypertrophy <- Data_RR %>% 
  filter(!is.na(lnRR_yi) & !is.na(lnRR_vi) & outcome == "hypertrophy")

MultiLevelModel_RR_hypertrophy <- rma.mv(lnRR_yi, V=lnRR_vi, data=Data_RR_hypertrophy,
                                      slab=paste(label),
                                      random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                      control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_RR_hypertrophy, file = "models/MultiLevelModel_RR_hypertrophy")


### Calculate I^2 
I2_RR_hypertrophy <- i2_ml(MultiLevelModel_RR_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_RR_hypertrophy <- robust(MultiLevelModel_RR_hypertrophy, Data_RR_hypertrophy$study)

save(RobuEstMultiLevelModel_RR_hypertrophy, file = "models/RobuEstMultiLevelModel_RR_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_RR_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_RR_hypertrophy$b[1] + (RobuEstMultiLevelModel_RR_hypertrophy$se*1.96),
                                        RobuEstMultiLevelModel_RR_hypertrophy$b[1],
                                        RobuEstMultiLevelModel_RR_hypertrophy$b[1] - (RobuEstMultiLevelModel_RR_hypertrophy$se*1.96),
                                        RobuEstMultiLevelModel_RR_hypertrophy$b[1]),
                                  y = c(-15,-25,-15,-5))

# Prediction interval
PI_RR_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_RR_hypertrophy))

# I^2 labels
I2_RR_hypertrophy_lab <- data.frame(level = c("study", "arm", "es"),
                                 I2 = I2_RR_hypertrophy[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_RR_hypertrophy <- Data_RR_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(lnRR_yi)]),
         se = sqrt(lnRR_vi)) %>%
  ggplot(aes(x=(exp(lnRR_yi)-1)*100, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  # scale_x_continuous(limits = c(-20,40), breaks = c(-20,-10,0,10,20,30,40)) +
  geom_linerange(aes(xmin = (exp(lnRR_yi)-1)*100 - (se*1.96), xmax = (exp(lnRR_yi)-1)*100 + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_RR_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round((exp(pred)-1)*100,2)} [95% Confidence Interval: {round((exp(ci.lb)-1)*100,2)} to {round((exp(ci.ub)-1)*100,2)}]"),
                x = 40, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_RR_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round((exp(pi.lb)-1)*100,2)} to {round((exp(pi.ub)-1)*100,2)}]"),
                x = 40, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_RR_hypertrophy_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 40, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_RR_hypertrophy, aes(y=-15, yend=-15, x=(exp(pi.lb)-1)*100, xend=(exp(pi.ub)-1)*100),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_RR_hypertrophy, aes(x=(exp(x)-1)*100,y=y)) +
  labs(y = "",
       x = "Exponentiated Response Ratio (%; Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_RR_hypertrophy, file = "plots/forest_RR_hypertrophy")

### Combine Plots
forest_RR_plots <- (forest_RR_strength / forest_RR_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_RR_plots, file = "plots/forest_RR_plots")

forest_RR_plots

ggsave("plots/forest_RR_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)


###### Comparison of variance of change between RT and CON ######

### Using the SDir - standard deviation for individual response
Data_SDir <- Data %>%
  filter(!is.na(RT_delta_sd) | !is.na(CON_delta_sd) | CON_delta_sd == 0) %>%
  mutate(SDir = sqrt(pmax(0, RT_delta_sd^2 - CON_delta_sd^2)),
         SDir_se = sqrt(2*(((RT_delta_sd^4)/(RT_n-1))+((CON_delta_sd^4)/(CON_n-1))))) %>%
  filter(!is.na(SDir))

### Strength
Data_SDir_strength <- Data_SDir %>% 
  filter(outcome == "strength")

MultiLevelModel_SDir_strength <- rma.mv(SDir, V=SDir_se^2, data=Data_SDir_strength,
                                         slab=paste(label),
                                         random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SDir_strength, file = "models/MultiLevelModel_SDir_strength")


### Calculate I^2 
I2_SDir_strength <- i2_ml(MultiLevelModel_SDir_strength)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SDir_strength <- robust(MultiLevelModel_SDir_strength, Data_SDir_strength$study)

save(RobuEstMultiLevelModel_SDir_strength, file = "models/RobuEstMultiLevelModel_SDir_strength")

### Caterpillar plot 

# Overall estimate
diamond_SDir_strength <- data.frame(x = c(RobuEstMultiLevelModel_SDir_strength$b[1] + (RobuEstMultiLevelModel_SDir_strength$se*1.96),
                                           RobuEstMultiLevelModel_SDir_strength$b[1],
                                           RobuEstMultiLevelModel_SDir_strength$b[1] - (RobuEstMultiLevelModel_SDir_strength$se*1.96),
                                           RobuEstMultiLevelModel_SDir_strength$b[1]),
                                     y = c(-15,-25,-15,-5))

# Prediction interval
PI_SDir_strength <- as.data.frame(predict(RobuEstMultiLevelModel_SDir_strength))

# I^2 labels
I2_SDir_strength_lab <- data.frame(level = c("study", "arm", "es"),
                                    I2 = I2_SDir_strength[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_SDir_strength <- Data_SDir_strength %>% 
  mutate(es = factor(es, levels = es[order(SDir)])) %>%
  ggplot(aes(x=SDir, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  geom_linerange(aes(xmin = SDir - (SDir_se*1.96), xmax = SDir + (SDir_se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_SDir_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 500000, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_SDir_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 500000, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_SDir_strength_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 500000, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_SDir_strength, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_SDir_strength, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standard Deviation for Individual Responses (Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_SDir_strength, file = "plots/forest_SDir_strength")

### Hypertrophy
Data_SDir_hypertrophy <- Data_SDir %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_SDir_hypertrophy <- rma.mv(SDir, V=SDir_se^2, data=Data_SDir_hypertrophy,
                                            slab=paste(label),
                                            random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                            control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SDir_hypertrophy, file = "models/MultiLevelModel_SDir_hypertrophy")

### Calculate I^2 
I2_SDir_hypertrophy <- i2_ml(MultiLevelModel_SDir_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SDir_hypertrophy <- robust(MultiLevelModel_SDir_hypertrophy, Data_SDir_hypertrophy$study)

save(RobuEstMultiLevelModel_SDir_hypertrophy, file = "models/RobuEstMultiLevelModel_SDir_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_SDir_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_SDir_hypertrophy$b[1] + (RobuEstMultiLevelModel_SDir_hypertrophy$se*1.96),
                                          RobuEstMultiLevelModel_SDir_hypertrophy$b[1],
                                          RobuEstMultiLevelModel_SDir_hypertrophy$b[1] - (RobuEstMultiLevelModel_SDir_hypertrophy$se*1.96),
                                          RobuEstMultiLevelModel_SDir_hypertrophy$b[1]),
                                    y = c(-15,-25,-15,-5))

# Prediction interval
PI_SDir_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_SDir_hypertrophy))

# I^2 labels
I2_SDir_hypertrophy_lab <- data.frame(level = c("study", "arm", "es"),
                                   I2 = I2_SDir_hypertrophy[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_SDir_hypertrophy <- Data_SDir_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(SDir)])) %>%
  ggplot(aes(x=SDir, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  geom_linerange(aes(xmin = SDir - (SDir_se*1.96), xmax = SDir + (SDir_se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_SDir_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 4e+06, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_SDir_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 4e+06, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_SDir_hypertrophy_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 4e+06, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_SDir_hypertrophy, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub), 
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_SDir_hypertrophy, aes(x=x,y=y)) +
  labs(y = "",
       x = "Standard Deviation for Individual Responses (Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_SDir_hypertrophy, file = "plots/forest_SDir_hypertrophy")

### Combine Plots
forest_SDir_plots <- (forest_SDir_strength / forest_SDir_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_SDir_plots, file = "plots/forest_SDir_plots")

forest_SDir_plots

ggsave("plots/forest_SDir_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)


### Using the Log Variability Ratio
Data_logVR <- escalc(measure = "VR", sd1i = RT_delta_sd, sd2i = CON_delta_sd, n1i = RT_n, n2i = CON_n, data = Data)
Data_logVR <- Data_logVR %>% 
    filter(!is.na(yi))

### Strength
Data_logVR_strength <- Data_logVR %>% 
  filter(outcome == "strength")

MultiLevelModel_logVR_strength <- rma.mv(yi, V=vi, data=Data_logVR_strength,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logVR_strength, file = "models/MultiLevelModel_logVR_strength")

### Calculate I^2 
I2_logVR_strength <- i2_ml(MultiLevelModel_logVR_strength)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_strength <- robust(MultiLevelModel_logVR_strength, Data_logVR_strength$study)

save(RobuEstMultiLevelModel_logVR_strength, file = "models/RobuEstMultiLevelModel_logVR_strength")

### Caterpillar plot 

# Overall estimate
diamond_logVR_strength <- data.frame(x = c(RobuEstMultiLevelModel_logVR_strength$b[1] + (RobuEstMultiLevelModel_logVR_strength$se*1.96),
                                          RobuEstMultiLevelModel_logVR_strength$b[1],
                                          RobuEstMultiLevelModel_logVR_strength$b[1] - (RobuEstMultiLevelModel_logVR_strength$se*1.96),
                                          RobuEstMultiLevelModel_logVR_strength$b[1]),
                                    y = c(-15,-25,-15,-5))

# Prediction interval
PI_logVR_strength <- as.data.frame(predict(RobuEstMultiLevelModel_logVR_strength))

# I^2 labels
I2_logVR_strength_lab <- data.frame(level = c("study", "arm", "es"),
                                   I2 = I2_logVR_strength[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_logVR_strength <- Data_logVR_strength %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-4,6), breaks = c(-3,-2,-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_logVR_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 4.5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_logVR_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 4.5, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_logVR_strength_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 4.5, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_logVR_strength, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_logVR_strength, aes(x=x,y=y)) +
  labs(y = "",
       x = "Log Variability Ratio (Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_logVR_strength, file = "plots/forest_logVR_strength")

### Hypertrophy
Data_logVR_hypertrophy <- Data_logVR %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_logVR_hypertrophy <- rma.mv(yi, V=vi, data=Data_logVR_hypertrophy,
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                            control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logVR_hypertrophy, file = "models/MultiLevelModel_logVR_hypertrophy")

### Calculate I^2 
I2_logVR_hypertrophy <- i2_ml(MultiLevelModel_logVR_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_hypertrophy <- robust(MultiLevelModel_logVR_hypertrophy, Data_logVR_hypertrophy$study)

save(RobuEstMultiLevelModel_logVR_hypertrophy, file = "models/RobuEstMultiLevelModel_logVR_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_logVR_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_logVR_hypertrophy$b[1] + (RobuEstMultiLevelModel_logVR_hypertrophy$se*1.96),
                                           RobuEstMultiLevelModel_logVR_hypertrophy$b[1],
                                           RobuEstMultiLevelModel_logVR_hypertrophy$b[1] - (RobuEstMultiLevelModel_logVR_hypertrophy$se*1.96),
                                           RobuEstMultiLevelModel_logVR_hypertrophy$b[1]),
                                     y = c(-15,-25,-15,-5))

# Prediction interval
PI_logVR_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_logVR_hypertrophy))

# I^2 labels
I2_logVR_hypertrophy_lab <- data.frame(level = c("study", "arm", "es"),
                                    I2 = I2_logVR_hypertrophy[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_logVR_hypertrophy <- Data_logVR_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  scale_x_continuous(limits = c(-4,6), breaks = c(-3,-2,-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_logVR_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 4.5, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_logVR_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 4.5, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_logVR_hypertrophy_lab,
            aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                x = 4.5, y = 70), size = 3,
            fill = NA, label.color = NA, # remove background and outline
            label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_logVR_hypertrophy, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_logVR_hypertrophy, aes(x=x,y=y)) +
  labs(y = "",
       x = "Log Variability Ratio (Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_logVR_hypertrophy, file = "plots/forest_logVR_hypertrophy")

### Combine Plots
forest_logVR_plots <- (forest_logVR_strength / forest_logVR_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_logVR_plots, file = "plots/forest_logVR_plots")

forest_logVR_plots

ggsave("plots/forest_logVR_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)

### Combine SMD, RR, SDir, and logVR estimates

SMD_RR_SDir_logVR <- rbind(data.frame(Outcome = "strength",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_strength$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_strength, 2)[2],
                                    `I^2 arm` = round(I2_SMD_strength, 2)[3],
                                    `I^2 effect` = round(I2_SMD_strength, 2)[4]),
                        data.frame(Outcome = "strength",
                                   Model = "RR",
                                   Moderator = "Main model",
                                   Estimate = round((exp(RobuEstMultiLevelModel_RR_strength$b)-1)*100, 2),
                                   Lower = round((exp(RobuEstMultiLevelModel_RR_strength$ci.lb)-1)*100, 2),
                                   Upper = round((exp(RobuEstMultiLevelModel_RR_strength$ci.ub)-1)*100, 2),
                                   `I^2 study` = round(I2_RR_strength, 2)[2],
                                   `I^2 arm` = round(I2_RR_strength, 2)[3],
                                   `I^2 effect` = round(I2_RR_strength, 2)[4]),
                        data.frame(Outcome = "strength",
                                   Model = "SDir",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_SDir_strength$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_SDir_strength$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_SDir_strength$ci.ub, 2),
                                   `I^2 study` = round(I2_SDir_strength, 2)[2],
                                   `I^2 arm` = round(I2_SDir_strength, 2)[3],
                                   `I^2 effect` = round(I2_SDir_strength, 2)[4]),
                        data.frame(Outcome = "strength",
                                   Model = "logVR",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_logVR_strength$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_logVR_strength$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_logVR_strength$ci.ub, 2),
                                   `I^2 study` = round(I2_logVR_strength, 2)[2],
                                   `I^2 arm` = round(I2_logVR_strength, 2)[3],
                                   `I^2 effect` = round(I2_logVR_strength, 2)[4]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "SMD",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy$ci.ub, 2),
                                   `I^2 study` = round(I2_SMD_hypertrophy, 2)[2],
                                   `I^2 arm` = round(I2_SMD_hypertrophy, 2)[3],
                                   `I^2 effect` = round(I2_SMD_hypertrophy, 2)[4]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "RR",
                                   Moderator = "Main model",
                                   Estimate = round((exp(RobuEstMultiLevelModel_RR_hypertrophy$b)-1)*100, 2),
                                   Lower = round((exp(RobuEstMultiLevelModel_RR_hypertrophy$ci.lb)-1)*100, 2),
                                   Upper = round((exp(RobuEstMultiLevelModel_RR_hypertrophy$ci.ub)-1)*100, 2),
                                   `I^2 study` = round(I2_RR_hypertrophy, 2)[2],
                                   `I^2 arm` = round(I2_RR_hypertrophy, 2)[3],
                                   `I^2 effect` = round(I2_RR_hypertrophy, 2)[4]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "SDir",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_SDir_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_SDir_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_SDir_hypertrophy$ci.ub, 2),
                                   `I^2 study` = round(I2_SDir_hypertrophy, 2)[2],
                                   `I^2 arm` = round(I2_SDir_hypertrophy, 2)[3],
                                   `I^2 effect` = round(I2_SDir_hypertrophy, 2)[4]),
                        data.frame(Outcome = "hypertrophy",
                                   Model = "logVR",
                                   Moderator = "Main model",
                                   Estimate = round(RobuEstMultiLevelModel_logVR_hypertrophy$b, 2),
                                   Lower = round(RobuEstMultiLevelModel_logVR_hypertrophy$ci.lb, 2),
                                   Upper = round(RobuEstMultiLevelModel_logVR_hypertrophy$ci.ub, 2),
                                   `I^2 study` = round(I2_logVR_hypertrophy, 2)[2],
                                   `I^2 arm` = round(I2_logVR_hypertrophy, 2)[3],
                                   `I^2 effect` = round(I2_logVR_hypertrophy, 2)[4])
)

save(SMD_RR_SDir_logVR, file = "models/SMD_RR_SDir_logVR")

###### Section - Mean-variance relationships in muscular strength, endurance, and hypertrophy
### We'll look at just the baseline strength data from the Polito et al. data
Data_long_pre_m <- Data %>%
  select(study, arm, es, RT_n, CON_n, RT_pre_m, CON_pre_m, outcome) %>%
  pivot_longer(c(RT_pre_m, CON_pre_m),
                   names_to = "group",
                   values_to = "mean")

Data_long_pre_m$group <- recode(Data_long_pre_m$group, RT_pre_m = "RT", CON_pre_m = "CON")

Data_long_pre_sd <- Data %>%
  select(study, arm, es, RT_n, CON_n, RT_pre_sd, CON_pre_sd, outcome) %>%
  pivot_longer(c(RT_pre_sd, CON_pre_sd),
                   names_to = "group",
                   values_to = "sd")

Data_long_pre_sd$group <- recode(Data_long_pre_sd$group, RT_pre_sd = "RT", CON_pre_sd = "CON")

Data_long_pre <- cbind(Data_long_pre_m, sd = Data_long_pre_sd$sd) 

Data_long_pre_RT <- Data_long_pre %>%
  filter(group == "RT") %>%
  mutate(n = RT_n,
         arm = as.factor(unclass(factor(unlist(arm)))),
         es = as.factor(unclass(factor(unlist(es))))) %>%
  select(study, arm, es, outcome, group, mean, sd, n)

Data_long_pre_CON <- Data_long_pre %>%
  filter(group == "CON") %>%
  distinct(study, outcome, group, mean, sd, .keep_all = TRUE) %>%
  mutate(n = CON_n,
         arm = as.factor(unclass(factor(unlist(arm)))+length(unique(Data_long_pre_RT$arm))),
         es = as.factor(unclass(factor(unlist(es)))+length(unique(Data_long_pre_RT$es)))) %>%
  select(study, arm, es, outcome, group, mean, sd, n)

Data_long_pre <- rbind(Data_long_pre_RT, Data_long_pre_CON) %>%
  filter(!is.na(mean) &
           !is.na(sd))

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

### Comparing a range of models

# Random intercepts only

MultiLevelModel_ri_only_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                                mods = ~ log(mean) + outcome,
                                                method="REML", test="t",
                                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_only_log_mean_variance, file = "models/MultiLevelModel_ri_only_log_mean_variance")

### Calculate I^2 
I2_ri_only_log_mean_variance <- i2_ml(MultiLevelModel_ri_only_log_mean_variance)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_only_log_mean_variance <- robust(MultiLevelModel_ri_only_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_only_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_only_log_mean_variance")

# Random intercepts plus slope for outcome in study

MultiLevelModel_ri_group_slope_study_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                     random = list(~ outcome | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                     mods = ~ log(mean) + outcome,
                                                     method="REML", test="t",
                                                     # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_group_slope_study_log_mean_variance, file = "models/MultiLevelModel_ri_group_slope_study_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance <- robust(MultiLevelModel_ri_group_slope_study_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance")

# Random intercepts plus slope for outcome in study and arm

MultiLevelModel_ri_group_slope_study_arm_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                random = list(~ outcome | study, ~ outcome | arm, ~ 1 | es), struct = "GEN",
                                                mods = ~ log(mean) + outcome,
                                                method="REML", test="t",
                                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
                                            )

save(MultiLevelModel_ri_group_slope_study_arm_log_mean_variance, file = "models/MultiLevelModel_ri_group_slope_study_arm_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance <- robust(MultiLevelModel_ri_group_slope_study_arm_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance")

# Random intercepts plus slope for log mean in study

MultiLevelModel_ri_mean_slope_study_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                                     random = list(~ log(mean) | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                     mods = ~ log(mean) + outcome,
                                                                     method="REML", test="t",
                                                                     # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_log_mean_variance, file = "models/MultiLevelModel_ri_mean_slope_study_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance <- robust(MultiLevelModel_ri_mean_slope_study_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance")

# Random intercepts plus slope for log mean in study and arm

MultiLevelModel_ri_mean_slope_study_arm_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                                random = list(~ log(mean) | study, ~ log(mean) | arm, ~ 1 | es), struct = "GEN",
                                                                mods = ~ log(mean) + outcome,
                                                                method="REML", test="t",
                                                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_arm_log_mean_variance, file = "models/MultiLevelModel_ri_mean_slope_study_arm_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance <- robust(MultiLevelModel_ri_mean_slope_study_arm_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance")

# Random intercepts plus slope for log mean and outcome in study

MultiLevelModel_ri_mean_group_slope_study_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                                    random = list(~ log(mean) + outcome | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                    mods = ~ log(mean) + outcome,
                                                                    method="REML", test="t",
                                                                    # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_log_mean_variance, file = "models/MultiLevelModel_ri_mean_group_slope_study_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance <- robust(MultiLevelModel_ri_mean_group_slope_study_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance")

# Random intercepts plus slope for log mean and group in study and arm

MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_pre,
                                                                      random = list(~ log(mean) + outcome | study, ~ log(mean) + outcome | arm, ~ 1 | es), struct = "GEN",
                                                                      mods = ~ log(mean) + outcome,
                                                                      method="REML", test="t",
                                                                      # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance, file = "models/MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance <- robust(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance, Data_long_pre$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance")

### Compare all models against one another using Bayes Factors i.e., under which model is the data more probable?

BF_variance_models <- bayesfactor_models(RobuEstMultiLevelModel_ri_only_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance,
                                     RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance)

save(BF_variance_models, file = "models/BF_variance_models")

BF_2log <- function(x) (2*x)

BF_variance_models <- as.data.frame(as.matrix(BF_variance_models))  %>%
  mutate_at(1:7, BF_2log) %>%
  rownames_to_column("Denominator") %>%
  rename("Random Intercepts Only" = "RobuEstMultiLevelModel_ri_only_log_mean_variance",
         "Random Intercepts + Random Slope (Outcome) at Study Level" = "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance",
         "Random Intercepts + Random Slope (Outcome) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance",
         "Random Intercepts + Random Slope (log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance",
         "Random Intercepts + Random Slope (log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance",
         "Random Intercepts + Random Slope (Outcome + log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance",
         "Random Intercepts + Random Slope (Outcome + log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance") %>%
  pivot_longer(2:8, names_to = "Numerator", values_to = "logBF")

BF_variance_models$Denominator <-  recode(BF_variance_models$Denominator, 
                                      "RobuEstMultiLevelModel_ri_only_log_mean_variance" = "Random Intercepts Only",
                                      "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_variance" = "Random Intercepts + Random Slope (Outcome) at Study Level",
                                      "RobuEstMultiLevelModel_ri_group_slope_study_arm_log_mean_variance" = "Random Intercepts + Random Slope (Outcome) at Study and Arm Level",
                                      "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_variance" = "Random Intercepts + Random Slope (log Mean) at Study Level",
                                      "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_variance" = "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
                                      "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_variance" = "Random Intercepts + Random Slope (Outcome + log Mean) at Study Level",
                                      "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_variance" = "Random Intercepts + Random Slope (Outcome + log Mean) at Study and Arm Level")

model_mean_variance_pre_model_comparisons <- BF_variance_models %>% 
  mutate(Denominator = factor(Denominator, levels= c( 
    "Random Intercepts Only",
    "Random Intercepts + Random Slope (Outcome) at Study Level",
    "Random Intercepts + Random Slope (Outcome) at Study and Arm Level",
    "Random Intercepts + Random Slope (log Mean) at Study Level",
    "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
    "Random Intercepts + Random Slope (Outcome + log Mean) at Study Level",
    "Random Intercepts + Random Slope (Outcome + log Mean) at Study and Arm Level")),
    Numerator = factor(Numerator, levels= c( 
      "Random Intercepts Only",
      "Random Intercepts + Random Slope (Outcome) at Study Level",
      "Random Intercepts + Random Slope (Outcome) at Study and Arm Level",
      "Random Intercepts + Random Slope (log Mean) at Study Level",
      "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
      "Random Intercepts + Random Slope (Outcome + log Mean) at Study Level",
      "Random Intercepts + Random Slope (Outcome + log Mean) at Study and Arm Level"))) %>%
  ggplot(aes(x=Numerator, y=Denominator, fill=logBF)) +
  geom_tile() +
  geom_raster() +
  geom_text(aes(label = round(logBF,2))) +
  scale_fill_gradient2(low = "#E69F00", mid="white", high = "#56B4E9") +
  scale_y_discrete(limits=rev, labels = function(x) str_wrap(x, width = 25)) +
  scale_x_discrete(position = "top", labels = function(x) str_wrap(x, width = 25)) +
  labs(title = "Testing model specification for log standard deviation models of pre-intervention scores using 2×log(BF)",
       fill = "2×log(BF)",
       caption = "Kass and Raferty (1995) scale:
       -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5))

save(model_mean_variance_pre_model_comparisons, file = "plots/model_mean_variance_pre_model_comparisons")

model_mean_variance_pre_model_comparisons

ggsave("plots/mean_variance_pre_model_comparisons.tiff", width = 15, height = 7.5, device = "tiff", dpi = 300)


### Meta-analytic scatter plot

# get the predicted log values
Data_long_pre <- Data_long_pre %>%
  filter(!is.na(SD_log))

Data_long_pre <- cbind(Data_long_pre, pred = predict(RobuEstMultiLevelModel_ri_only_log_mean_variance)$pred,
                      ci.lb =  predict(RobuEstMultiLevelModel_ri_only_log_mean_variance)$ci.lb,
                      ci.ub =  predict(RobuEstMultiLevelModel_ri_only_log_mean_variance)$ci.ub) %>%
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
  select(study, arm, es, RT_n, CON_n, RT_delta_m, CON_delta_m, outcome) %>%
  pivot_longer(c(RT_delta_m, CON_delta_m),
               names_to = "group",
               values_to = "mean")

Data_long_m$group <- recode(Data_long_m$group, RT_delta_m = "RT", CON_delta_m = "CON")

Data_long_sd <- Data %>%
  select(study, arm, es, RT_n, CON_n, RT_delta_sd, CON_delta_sd, outcome) %>%
  pivot_longer(c(RT_delta_sd, CON_delta_sd),
               names_to = "group",
               values_to = "sd")

Data_long_sd$group <- recode(Data_long_sd$group, RT_delta_sd = "RT", CON_delta_sd = "CON")

Data_long <- cbind(Data_long_m, sd = Data_long_sd$sd) 

Data_long_RT <- Data_long %>%
  filter(group == "RT") %>%
  mutate(n = RT_n,
         arm = as.factor(unclass(factor(unlist(arm)))),
         es = as.factor(unclass(factor(unlist(es))))) %>%
  select(study, arm, es, outcome, group, mean, sd, n)

Data_long_CON <- Data_long %>%
  filter(group == "CON") %>%
  distinct(study, outcome, group, mean, sd, .keep_all = TRUE) %>%
  mutate(n = CON_n,
         arm = as.factor(unclass(factor(unlist(arm)))+length(unique(Data_long_RT$arm))),
         es = as.factor(unclass(factor(unlist(es)))+length(unique(Data_long_RT$es)))) %>%
  select(study, arm, es, outcome, group, mean, sd, n)

Data_long <- rbind(Data_long_RT, Data_long_CON)

# Calculate log SD and variance of log SD
Data_long$SD_log <- log(Data_long$sd) + (1/(2*(Data_long$n-1)))
Data_long$SD_log_vi <- (1/(2*(Data_long$n-1)))

# Plot raw mean and SD
m_sd_strength <- ggplot(subset(Data_long, outcome == "strength"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd, color = group), alpha = 0.2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_y_continuous(limits = c(0,1000)) +
  labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Group") +
  guides(fill = "none") +
  theme_classic() +
  ggtitle("Strength Outcomes")

m_sd_hypertrophy <- ggplot(subset(Data_long, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean, y=sd, color = group), alpha = 0.2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Group") +
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
  geom_point(aes(x=log(mean), y=SD_log, color = group), alpha = 0.2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Group") +
  guides(fill = "none") +
  theme_classic() 

m_sd_hypertrophy_log <- ggplot(subset(Data_long, outcome == "hypertrophy"), aes(x=log(mean), y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=log(mean), y=SD_log, color = group), alpha = 0.2) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Group") +
  guides(fill = "none") +
  theme_classic() 

# Plot together
mean_variance_delta_plots <- ((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A") 

save(mean_variance_delta_plots, file = "plots/mean_variance_delta_plots")

mean_variance_delta_plots

ggsave("plots/mean_variance_delta_plots.tiff", width = 10, height = 7.5, device = "tiff", dpi = 300)

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
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength, file = "models/MultiLevelModel_logCVR_strength")

### Calculate I^2 
I2_logCVR_strength <- i2_ml(MultiLevelModel_logCVR_strength)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength <- robust(MultiLevelModel_logCVR_strength, Data_logCVR_strength$study)

save(RobuEstMultiLevelModel_logCVR_strength, file = "models/RobuEstMultiLevelModel_logCVR_strength")

### Caterpillar plot 

# Overall estimate
diamond_logCVR_strength <- data.frame(x = c(RobuEstMultiLevelModel_logCVR_strength$b[1] + (RobuEstMultiLevelModel_logCVR_strength$se*1.96),
                                           RobuEstMultiLevelModel_logCVR_strength$b[1],
                                           RobuEstMultiLevelModel_logCVR_strength$b[1] - (RobuEstMultiLevelModel_logCVR_strength$se*1.96),
                                           RobuEstMultiLevelModel_logCVR_strength$b[1]),
                                     y = c(-15,-25,-15,-5))

# Prediction interval
PI_logCVR_strength <- as.data.frame(predict(RobuEstMultiLevelModel_logCVR_strength))

# I^2 labels
I2_logCVR_strength_lab <- data.frame(level = c("study", "arm", "es"),
                                    I2 = I2_logCVR_strength[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_logCVR_strength <- Data_logCVR_strength %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  # scale_x_continuous(limits = c(-4,6), breaks = c(-3,-2,-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_logCVR_strength,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 60, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_logCVR_strength,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 60, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_logCVR_strength_lab,
                aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                    x = 60, y = 70), size = 3,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_logCVR_strength, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_logCVR_strength, aes(x=x,y=y)) +
  labs(y = "",
       x = "Log Coefficient of Variability Ratio (Positive Values Favour Resistance Training)",
       title = "Strength Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_logCVR_strength, file = "plots/forest_logCVR_strength")

### Hypertrophy
Data_logCVR_hypertrophy <- Data_logCVR %>% 
  filter(outcome == "hypertrophy")

MultiLevelModel_logCVR_hypertrophy <- rma.mv(yi, V=vi, data=Data_logCVR_hypertrophy,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy, file = "models/MultiLevelModel_logCVR_hypertrophy")

### Calculate I^2 
I2_logCVR_hypertrophy <- i2_ml(MultiLevelModel_logCVR_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy <- robust(MultiLevelModel_logCVR_hypertrophy, Data_logCVR_hypertrophy$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy")

### Caterpillar plot 

# Overall estimate
diamond_logCVR_hypertrophy <- data.frame(x = c(RobuEstMultiLevelModel_logCVR_hypertrophy$b[1] + (RobuEstMultiLevelModel_logCVR_hypertrophy$se*1.96),
                                            RobuEstMultiLevelModel_logCVR_hypertrophy$b[1],
                                            RobuEstMultiLevelModel_logCVR_hypertrophy$b[1] - (RobuEstMultiLevelModel_logCVR_hypertrophy$se*1.96),
                                            RobuEstMultiLevelModel_logCVR_hypertrophy$b[1]),
                                      y = c(-15,-25,-15,-5))

# Prediction interval
PI_logCVR_hypertrophy <- as.data.frame(predict(RobuEstMultiLevelModel_logCVR_hypertrophy))

# I^2 labels
I2_logCVR_hypertrophy_lab <- data.frame(level = c("study", "arm", "es"),
                                     I2 = I2_logCVR_hypertrophy[2:4]) %>%
  pivot_wider(names_from = "level", values_from = "I2")

# Plot
forest_logCVR_hypertrophy <- Data_logCVR_hypertrophy %>% 
  mutate(es = factor(es, levels = es[order(yi)]),
         se = sqrt(vi)) %>%
  ggplot(aes(x=yi, y=fct_rev(es))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size=0.5, alpha=0.25) +
  # scale_x_continuous(limits = c(-4,6), breaks = c(-3,-2,-1,0,1,2,3,4,5,6)) +
  geom_linerange(aes(xmin = yi - (se*1.96), xmax = yi + (se*1.96)), alpha=0.25) +
  geom_text(data = mutate_if(PI_logCVR_hypertrophy,
                             is.numeric, round, 2),
            aes(label = glue::glue("Overall Estimate = {round(pred,2)} [95% Confidence Interval: {round(ci.lb,2)} to {round(ci.ub,2)}]"),
                x = 120, y = 110), hjust = "centre", size = 3) +
  geom_text(data = PI_logCVR_hypertrophy,
            aes(label = glue::glue("[95% Prediction Interval: {round(pi.lb,2)} to {round(pi.ub,2)}]"),
                x = 120, y = 90), hjust = "centre", size = 3) +
  geom_richtext(data = I2_logCVR_hypertrophy_lab,
                aes(label = glue::glue("*I*<sup>2</sup> [Study = {round(study,1)}%; Arm = {round(arm,1)}%; Effect = {round(es,1)}%]"),
                    x = 120, y = 70), size = 3,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt") # remove padding
  ) +
  geom_segment(data = PI_logCVR_hypertrophy, aes(y=-15, yend=-15, x=pi.lb, xend=pi.ub),
               arrow = arrow(length=unit(0.30,"cm"), angle = 90, ends="both", type = "open")) +
  geom_polygon(data=diamond_logCVR_hypertrophy, aes(x=x,y=y)) +
  labs(y = "",
       x = "Log Coefficient of Variability Ratio (Positive Values Favour Resistance Training)",
       title = "Hypertrophy Outcomes") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

save(forest_logCVR_hypertrophy, file = "plots/forest_logCVR_hypertrophy")

### Combine Plots
forest_logCVR_plots <- (forest_logCVR_strength / forest_logCVR_hypertrophy) + 
  plot_annotation(tag_levels = "A")

save(forest_logCVR_plots, file = "plots/forest_logCVR_plots")

forest_logCVR_plots

ggsave("plots/forest_logCVR_plots.tiff", width = 10, height = 10, device = "tiff", dpi = 300)

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

### Comparing a range of models

# Random intercepts only

MultiLevelModel_ri_only_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                                    mods = ~ mean_log + group,
                                                    method="REML", test="t",
                                                    # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_only_log_mean_mod_strength, file = "models/MultiLevelModel_ri_only_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_only_log_mean_mod_strength <- robust(MultiLevelModel_ri_only_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_only_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_only_log_mean_mod_strength")

# Random intercepts plus slope for group in study

MultiLevelModel_ri_group_slope_study_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                                 random = list(~ group | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                 mods = ~ mean_log + group,
                                                                 method="REML", test="t",
                                                                 # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_group_slope_study_log_mean_mod_strength, file = "models/MultiLevelModel_ri_group_slope_study_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength <- robust(MultiLevelModel_ri_group_slope_study_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength")

# Random intercepts plus slope for log mean in study

MultiLevelModel_ri_mean_slope_study_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                                random = list(~ mean_log | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                mods = ~ mean_log + group,
                                                                method="REML", test="t",
                                                                # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_log_mean_mod_strength, file = "models/MultiLevelModel_ri_mean_slope_study_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength <- robust(MultiLevelModel_ri_mean_slope_study_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength")

# Random intercepts plus slope for log mean in study and arm

MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                                    random = list(~ mean_log | study, ~ mean_log | arm, ~ 1 | es), struct = "GEN",
                                                                    mods = ~ mean_log + group,
                                                                    method="REML", test="t",
                                                                    # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength, file = "models/MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength <- robust(MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength")

# Random intercepts plus slope for log mean and group in study

MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                                      random = list(~ mean_log + group | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                      mods = ~ mean_log + group,
                                                                      method="REML", test="t",
                                                                      # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength, file = "models/MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength <- robust(MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength")

# Random intercepts plus slope for log mean and group in study and arm

MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_strength,
                                                                          random = list(~ mean_log + group | study, ~ mean_log | arm, ~ 1 | es), struct = "GEN",
                                                                          mods = ~ mean_log + group,
                                                                          method="REML", test="t",
                                                                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength, file = "models/MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength <- robust(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength, Data_long_strength$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength")

### Compare all models against one another using Bayes Factors i.e., under which model is the data more probable?

BF_mod_strength_models <- bayesfactor_models(RobuEstMultiLevelModel_ri_only_log_mean_mod_strength,
                                         RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength,
                                         RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength,
                                         RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength,
                                         RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength,
                                         RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength)

save(BF_mod_strength_models, file = "models/BF_mod_strength_models")

BF_2log <- function(x) (2*x)

BF_mod_strength_models <- as.data.frame(as.matrix(BF_mod_strength_models))  %>%
  mutate_at(1:6, BF_2log) %>%
  rownames_to_column("Denominator") %>%
  rename("Random Intercepts Only" = "RobuEstMultiLevelModel_ri_only_log_mean_mod_strength",
         "Random Intercepts + Random Slope (group) at Study Level" = "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength",
         "Random Intercepts + Random Slope (log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength",
         "Random Intercepts + Random Slope (log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength",
         "Random Intercepts + Random Slope (group + log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength",
         "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength") %>%
  pivot_longer(2:7, names_to = "Numerator", values_to = "logBF")

BF_mod_strength_models$Denominator <-  recode(BF_mod_strength_models$Denominator, 
                                          "RobuEstMultiLevelModel_ri_only_log_mean_mod_strength" = "Random Intercepts Only",
                                          "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength" = "Random Intercepts + Random Slope (group) at Study Level",
                                          "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength" = "Random Intercepts + Random Slope (log Mean) at Study Level",
                                          "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength" = "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
                                          "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength" = "Random Intercepts + Random Slope (group + log Mean) at Study Level",
                                          "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength" = "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level")

model_mean_mod_strength_model_comparisons <- BF_mod_strength_models %>% 
  mutate(Denominator = factor(Denominator, levels= c( 
    "Random Intercepts Only",
    "Random Intercepts + Random Slope (group) at Study Level",
    "Random Intercepts + Random Slope (log Mean) at Study Level",
    "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
    "Random Intercepts + Random Slope (group + log Mean) at Study Level",
    "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level")),
    Numerator = factor(Numerator, levels= c( 
      "Random Intercepts Only",
      "Random Intercepts + Random Slope (group) at Study Level",
      "Random Intercepts + Random Slope (log Mean) at Study Level",
      "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
      "Random Intercepts + Random Slope (group + log Mean) at Study Level",
      "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level"))) %>%
  ggplot(aes(x=Numerator, y=Denominator, fill=logBF)) +
  geom_tile() +
  geom_raster() +
  geom_text(aes(label = round(logBF,2))) +
  scale_fill_gradient2(low = "#E69F00", mid="white", high = "#56B4E9") +
  scale_y_discrete(limits=rev, labels = function(x) str_wrap(x, width = 25)) +
  scale_x_discrete(position = "top", labels = function(x) str_wrap(x, width = 25)) +
  labs(title = "Testing model specification for log standard deviation models of change scores (strength) using 2×log(BF)",
       fill = "2×log(BF)",
       caption = "Kass and Raferty (1995) scale:
       -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5))

save(model_mean_mod_strength_model_comparisons, file = "plots/model_mean_mod_strength_model_comparisons")

model_mean_mod_strength_model_comparisons

ggsave("plots/mean_mod_strength_model_comparisons.tiff", width = 15, height = 7.5, device = "tiff", dpi = 300)

### Meta-analytic scatter plot

# get the predicted log values
Data_long_strength_log <- cbind(Data_long_strength, predict(RobuEstMultiLevelModel_ri_only_log_mean_mod_strength )) %>%
    mutate(wi = 1/sqrt(SD_log_vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))

model_m_sd_strength_log <- ggplot(Data_long_strength_log, aes(x=mean_log, y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean_log, y=SD_log, color = group, size = size), alpha = 0.2) +
    geom_ribbon(aes(x=mean_log, ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
    geom_line(aes(x=mean_log, y=pred, color = group)) +
    scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Group", shape = "", fill = "") +
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

### Comparing a range of models

# Random intercepts only

MultiLevelModel_ri_only_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                        random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es),
                                                        mods = ~ mean_log + group,
                                                        method="REML", test="t",
                                                        # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_only_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_only_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_ri_only_log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_only_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_only_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy")

# Random intercepts plus slope for group in study

MultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                                     random = list(~ group | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                     mods = ~ mean_log + group,
                                                                     method="REML", test="t",
                                                                     # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_rs_study__log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy")

# Random intercepts plus slope for log mean in study

MultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                                    random = list(~ mean_log | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                    mods = ~ mean_log + group,
                                                                    method="REML", test="t",
                                                                    # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_rs_both__log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy")

# Random intercepts plus slope for log mean in study and arm

MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                                        random = list(~ mean_log | study, ~ mean_log | arm, ~ 1 | es), struct = "GEN",
                                                                        mods = ~ mean_log + group,
                                                                        method="REML", test="t",
                                                                        # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_rs_both__log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy")

# Random intercepts plus slope for log mean and group in study

MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                                          random = list(~ mean_log + group | study, ~ 1 | arm, ~ 1 | es), struct = "GEN",
                                                                          mods = ~ mean_log + group,
                                                                          method="REML", test="t",
                                                                          # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_rs_both__log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy")

# Random intercepts plus slope for log mean and group in study and arm

MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=Data_long_hypertrophy,
                                                                              random = list(~ mean_log + group | study, ~ mean_log | arm, ~ 1 | es), struct = "GEN",
                                                                              mods = ~ mean_log + group,
                                                                              method="REML", test="t",
                                                                              # control=list(optimizer="optim", optmethod="Nelder-Mead")
)

save(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy, file = "models/MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy")

### Calculate I^2 
I2_rs_both__log_mean_mod_hypertrophy <- i2_ml(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy <- robust(MultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy, Data_long_hypertrophy$study)

save(RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy, file = "models/RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy")

### Compare all models against one another using Bayes Factors i.e., under which model is the data more probable?

BF_mod_hypertrophy_models <- bayesfactor_models(RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy,
                                             RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy,
                                             RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy,
                                             RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy,
                                             RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy,
                                             RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy)

save(BF_mod_hypertrophy_models, file = "models/BF_mod_hypertrophy_models")

BF_2log <- function(x) (2*x)

BF_mod_hypertrophy_models <- as.data.frame(as.matrix(BF_mod_hypertrophy_models))  %>%
  mutate_at(1:6, BF_2log) %>%
  rownames_to_column("Denominator") %>%
  rename("Random Intercepts Only" = "RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy",
         "Random Intercepts + Random Slope (group) at Study Level" = "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy",
         "Random Intercepts + Random Slope (log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy",
         "Random Intercepts + Random Slope (log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy",
         "Random Intercepts + Random Slope (group + log Mean) at Study Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy",
         "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level" = "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy") %>%
  pivot_longer(2:7, names_to = "Numerator", values_to = "logBF")

BF_mod_hypertrophy_models$Denominator <-  recode(BF_mod_hypertrophy_models$Denominator, 
                                              "RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy" = "Random Intercepts Only",
                                              "RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy" = "Random Intercepts + Random Slope (group) at Study Level",
                                              "RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy" = "Random Intercepts + Random Slope (log Mean) at Study Level",
                                              "RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy" = "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
                                              "RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy" = "Random Intercepts + Random Slope (group + log Mean) at Study Level",
                                              "RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy" = "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level")

model_mean_mod_hypertrophy_model_comparisons <- BF_mod_hypertrophy_models %>% 
  mutate(Denominator = factor(Denominator, levels= c( 
    "Random Intercepts Only",
    "Random Intercepts + Random Slope (group) at Study Level",
    "Random Intercepts + Random Slope (log Mean) at Study Level",
    "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
    "Random Intercepts + Random Slope (group + log Mean) at Study Level",
    "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level")),
    Numerator = factor(Numerator, levels= c( 
      "Random Intercepts Only",
      "Random Intercepts + Random Slope (group) at Study Level",
      "Random Intercepts + Random Slope (log Mean) at Study Level",
      "Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
      "Random Intercepts + Random Slope (group + log Mean) at Study Level",
      "Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level"))) %>%
  ggplot(aes(x=Numerator, y=Denominator, fill=logBF)) +
  geom_tile() +
  geom_raster() +
  geom_text(aes(label = round(logBF,2))) +
  scale_fill_gradient2(low = "#E69F00", mid="white", high = "#56B4E9") +
  scale_y_discrete(limits=rev, labels = function(x) str_wrap(x, width = 25)) +
  scale_x_discrete(position = "top", labels = function(x) str_wrap(x, width = 25)) +
  labs(title = "Testing model specification for log standard deviation models of change scores (hypertrophy) using 2×log(BF)",
       fill = "2×log(BF)",
       caption = "Kass and Raferty (1995) scale:
       -Inf to 0 = Negative; 0 to 2 = Weak; 2 to 6 = Positive; 6 to 10 = Strong; 10 to +Inf = Very Strong") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5))

save(model_mean_mod_hypertrophy_model_comparisons, file = "plots/model_mean_mod_hypertrophy_model_comparisons")

model_mean_mod_hypertrophy_model_comparisons

ggsave("plots/mean_mod_hypertrophy_model_comparisons.tiff", width = 15, height = 7.5, device = "tiff", dpi = 300)

### Meta-analytic scatter plot

# get the predicted log values
Data_long_hypertrophy_log <- cbind(Data_long_hypertrophy, predict(RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy )) %>%
  mutate(wi = 1/sqrt(SD_log_vi),
         size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi)))

model_m_sd_hypertrophy_log <- ggplot(Data_long_hypertrophy_log, aes(x=mean_log, y=SD_log)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.1) +
  geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
  geom_point(aes(x=mean_log, y=SD_log, color = group, size = size), alpha = 0.2) +
  geom_ribbon(aes(x=mean_log, ymax=ci.ub, ymin=ci.lb, fill = group), alpha = 0.2) +
  geom_line(aes(x=mean_log, y=pred, color = group)) +
  scale_fill_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  scale_color_manual("Group", values = alpha(c("Black", "#E69F00"),0.5)) +
  labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Group", shape = "", fill = "") +
  theme_classic() +
  ggtitle("Hypertrophy Outcomes") +
  guides(size = "none", fill = "none")

model_mean_variance_delta_plots <- (model_m_sd_strength_log | model_m_sd_hypertrophy_log) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = "A")

save(model_mean_variance_delta_plots, file = "plots/model_mean_variance_delta_plots")

model_mean_variance_delta_plots

ggsave("plots/model_mean_variance_delta_plots.tiff", width = 10, height = 5, device = "tiff", dpi = 300)


### Compare estimates from log CVR and each meta-regression model for strength and hypertrophy
model_estimates <- data.frame(Outcome = rep(""),
                 Model = c("lnCVR", 
                           "MLMR: Random Intercepts Only",
                           "MLMR: Random Intercepts + Random Slope (group) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (log Mean) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
                           "MLMR: Random Intercepts + Random Slope (group + log Mean) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level",
                           "lnCVR", 
                           "MLMR: Random Intercepts Only",
                           "MLMR: Random Intercepts + Random Slope (group) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (log Mean) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (log Mean) at Study and Arm Level",
                           "MLMR: Random Intercepts + Random Slope (group + log Mean) at Study Level",
                           "MLMR: Random Intercepts + Random Slope (group + log Mean) at Study and Arm Level"),
                 Estimate = c(RobuEstMultiLevelModel_logCVR_strength$b,
                              RobuEstMultiLevelModel_ri_only_log_mean_mod_strength$b[3], 
                              RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength$b[3], 
                              RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength$b[3],
                              RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength$b[3], 
                              RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength$b[3], 
                              RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength$b[3],
                              RobuEstMultiLevelModel_logCVR_hypertrophy$b,
                              RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy$b[3], 
                               RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy$b[3], 
                               RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy$b[3],
                               RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy$b[3], 
                               RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy$b[3], 
                               RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy$b[3]),
                 `Lower 95% CI` = c(RobuEstMultiLevelModel_logCVR_strength$ci.lb,
                                    RobuEstMultiLevelModel_ri_only_log_mean_mod_strength$ci.lb[3], 
                                    RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength$ci.lb[3], 
                                    RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength$ci.lb[3],
                                    RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength$ci.lb[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength$ci.lb[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength$ci.lb[3],
                                     RobuEstMultiLevelModel_logCVR_hypertrophy$ci.lb,
                                     RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy$ci.lb[3], 
                                     RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy$ci.lb[3], 
                                     RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy$ci.lb[3],
                                     RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy$ci.lb[3], 
                                     RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy$ci.lb[3], 
                                     RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy$ci.lb[3]),
                 `Upper 95% CI` = c(RobuEstMultiLevelModel_logCVR_strength$ci.ub,
                                    RobuEstMultiLevelModel_ri_only_log_mean_mod_strength$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_strength$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_strength$ci.ub[3],
                                    RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_strength$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_strength$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_strength$ci.ub[3],
                                    RobuEstMultiLevelModel_logCVR_hypertrophy$ci.ub,
                                    RobuEstMultiLevelModel_ri_only_log_mean_mod_hypertrophy$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_group_slope_study_log_mean_mod_hypertrophy$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_slope_study_log_mean_mod_hypertrophy$ci.ub[3],
                                    RobuEstMultiLevelModel_ri_mean_slope_study_arm_log_mean_mod_hypertrophy$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_log_mean_mod_hypertrophy$ci.ub[3], 
                                    RobuEstMultiLevelModel_ri_mean_group_slope_study_arm_log_mean_mod_hypertrophy$ci.ub[3])) %>%
  mutate(across(where(is.numeric), round, 2)) 

save(model_estimates, file = "models/group_model_estimates")

knitr::kable(
  model_estimates,
  caption = "Comparison of estimates from model using lnCVR and multilevel meta-regression models of variance~mean with group (RT vs CON) as a predictor",
) %>%
  row_spec(0, bold = TRUE) %>%
  pack_rows(index = c("Strength" = 7, "Hypertrophy" = 7)) %>%
  kable_classic(full_width = FALSE) %>%
  kable_styling() %>%
  save_kable("models/compare_group_model_estimates.html")

webshot::webshot("models/compare_group_model_estimates.html", "models/compare_group_model_estimates.pdf")

####### Summary tables for study/participant characteristics 

### Sample sizes
sample_sizes <- Data %>%
  select(arm, RT_n, CON_n) %>%
  group_by(arm) %>%
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
  pivot_longer(1:8, names_to = "Arm", values_to = "Sample Size")

save(sample_sizes, file = "models/sample_sizes")

### Descriptive Tables
Data_characteristics <- Data %>% 
  distinct(arm, .keep_all = TRUE) %>%
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
                                       RT_only = "RT Intervention Only?",
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

###### Stuff for supplementary materials etc.

###### Small study bias

# Most studies in our field care about detecting treatment effects, not variance
# So we'll explore small study bias first for these across all effects
# We'll create contour enhanced funnel plots for examination of publication bias 

# run that model
Data_SMD <- Data_SMD %>%
    filter(vi > 0) # filter NAs

MultiLevelModel_all <- rma.mv(yi, V=vi, data=Data_SMD,
                                            slab=paste(label),
                                            random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), method="REML")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_all <- robust(MultiLevelModel_all, Data_SMD$study)

### Contour enhanced funnel plot for examination of publication bias
tiff(here::here("plots","funnel_plot.tiff"), units="in", width=10, height=7.5, res=300)

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

# We don't have complete data on moderators for all studies, so we'll look at each separately

### TESTEX score
# SMD
# Strength
Data_SMD_strength_TESTEX <- Data_SMD_strength %>%
  filter(!is.na(TESTEX))

MultiLevelModel_SMD_strength_TESTEX <- rma.mv(yi, V=vi, data=Data_SMD_strength_TESTEX,
                                   slab=paste(label),
                                   random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ TESTEX, method="REML", test="t",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_TESTEX, file = "models/MultiLevelModel_SMD_strength_TESTEX")

### Calculate I^2 
I2_SMD_strength_TESTEX <- i2_ml(MultiLevelModel_SMD_strength_TESTEX)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_TESTEX <- robust(MultiLevelModel_SMD_strength_TESTEX, Data_SMD_strength_TESTEX$study)

save(RobuEstMultiLevelModel_SMD_strength_TESTEX, file = "models/RobuEstMultiLevelModel_SMD_strength_TESTEX")

# Hypertrophy
Data_SMD_hypertrophy_TESTEX <- Data_SMD_hypertrophy %>%
  filter(!is.na(TESTEX))

MultiLevelModel_SMD_hypertrophy_TESTEX <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_TESTEX,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ TESTEX, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_TESTEX, file = "models/MultiLevelModel_SMD_hypertrophy_TESTEX")

### Calculate I^2 
I2_SMD_hypertrophy_TESTEX <- i2_ml(MultiLevelModel_SMD_hypertrophy_TESTEX)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX <- robust(MultiLevelModel_SMD_hypertrophy_TESTEX, Data_SMD_hypertrophy_TESTEX$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX")

# logCVR
Data_logCVR_TESTEX <- Data_logCVR %>%
  filter(!is.na(TESTEX))

# Strength
MultiLevelModel_logCVR_strength_TESTEX <- rma.mv(yi, V=vi, data=subset(Data_logCVR_TESTEX, outcome == "strength"),
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ TESTEX, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_TESTEX, file = "models/MultiLevelModel_logCVR_strength_TESTEX")

### Calculate I^2 
I2_logCVR_strength_TESTEX <- i2_ml(MultiLevelModel_logCVR_strength_TESTEX)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_TESTEX <- robust(MultiLevelModel_logCVR_strength_TESTEX, subset(Data_logCVR_TESTEX, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_TESTEX, file = "models/RobuEstMultiLevelModel_logCVR_strength_TESTEX")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_TESTEX <- rma.mv(yi, V=vi, data=subset(Data_logCVR_TESTEX, outcome == "hypertrophy"),
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ TESTEX, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_TESTEX, file = "models/MultiLevelModel_logCVR_hypertrophy_TESTEX")

### Calculate I^2 
I2_logCVR_hypertrophy_TESTEX <- i2_ml(MultiLevelModel_logCVR_hypertrophy_TESTEX)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX <- robust(MultiLevelModel_logCVR_hypertrophy_TESTEX, subset(Data_logCVR_TESTEX, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX")

### Mean age
# SMD
# Strength
Data_SMD_strength_age <- Data_SMD_strength %>%
  filter(!is.na(age))

MultiLevelModel_SMD_strength_age <- rma.mv(yi, V=vi, data=Data_SMD_strength_age,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ age, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_age, file = "models/MultiLevelModel_SMD_strength_age")

### Calculate I^2 
I2_SMD_strength_age <- i2_ml(MultiLevelModel_SMD_strength_age)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_age <- robust(MultiLevelModel_SMD_strength_age, Data_SMD_strength_age$study)

save(RobuEstMultiLevelModel_SMD_strength_age, file = "models/RobuEstMultiLevelModel_SMD_strength_age")

# Hypertrophy
Data_SMD_hypertrophy_age <- Data_SMD_hypertrophy %>%
  filter(!is.na(age))

MultiLevelModel_SMD_hypertrophy_age <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_age,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ age, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_age, file = "models/MultiLevelModel_SMD_hypertrophy_age")

### Calculate I^2 
I2_SMD_hypertrophy_age <- i2_ml(MultiLevelModel_SMD_hypertrophy_age)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_age <- robust(MultiLevelModel_SMD_hypertrophy_age, Data_SMD_hypertrophy_age$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_age, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_age")

# logCVR
Data_logCVR_age <- Data_logCVR %>%
  filter(!is.na(age))

# Strength
MultiLevelModel_logCVR_strength_age <- rma.mv(yi, V=vi, data=subset(Data_logCVR_age, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ age, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_age, file = "models/MultiLevelModel_logCVR_strength_age")

### Calculate I^2 
I2_logCVR_strength_age <- i2_ml(MultiLevelModel_logCVR_strength_age)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_age <- robust(MultiLevelModel_logCVR_strength_age, subset(Data_logCVR_age, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_age, file = "models/RobuEstMultiLevelModel_logCVR_strength_age")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_age <- rma.mv(yi, V=vi, data=subset(Data_logCVR_age, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ age, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_age, file = "models/MultiLevelModel_logCVR_hypertrophy_age")

### Calculate I^2 
I2_logCVR_hypertrophy_age <- i2_ml(MultiLevelModel_logCVR_hypertrophy_age)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_age <- robust(MultiLevelModel_logCVR_hypertrophy_age, subset(Data_logCVR_age, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_age, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_age")

### Proportion of males
# SMD
# Strength
Data_SMD_strength_sex_._male <- Data_SMD_strength %>%
  filter(!is.na(sex_._male))

MultiLevelModel_SMD_strength_sex_._male <- rma.mv(yi, V=vi, data=Data_SMD_strength_sex_._male,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sex_._male, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_sex_._male, file = "models/MultiLevelModel_SMD_strength_sex_._male")

### Calculate I^2 
I2_SMD_strength_sex_._male <- i2_ml(MultiLevelModel_SMD_strength_sex_._male)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_sex_._male <- robust(MultiLevelModel_SMD_strength_sex_._male, Data_SMD_strength_sex_._male$study)

save(RobuEstMultiLevelModel_SMD_strength_sex_._male, file = "models/RobuEstMultiLevelModel_SMD_strength_sex_._male")

# Hypertrophy
Data_SMD_hypertrophy_sex_._male <- Data_SMD_hypertrophy %>%
  filter(!is.na(sex_._male))

MultiLevelModel_SMD_hypertrophy_sex_._male <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_sex_._male,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sex_._male, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_sex_._male, file = "models/MultiLevelModel_SMD_hypertrophy_sex_._male")

### Calculate I^2 
I2_SMD_hypertrophy_sex_._male <- i2_ml(MultiLevelModel_SMD_hypertrophy_sex_._male)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male <- robust(MultiLevelModel_SMD_hypertrophy_sex_._male, Data_SMD_hypertrophy_sex_._male$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male")

# logCVR
Data_logCVR_sex_._male <- Data_logCVR %>%
  filter(!is.na(sex_._male))

# Strength
MultiLevelModel_logCVR_strength_sex_._male <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sex_._male, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sex_._male, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_sex_._male, file = "models/MultiLevelModel_logCVR_strength_sex_._male")

### Calculate I^2 
I2_logCVR_strength_sex_._male <- i2_ml(MultiLevelModel_logCVR_strength_sex_._male)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_sex_._male <- robust(MultiLevelModel_logCVR_strength_sex_._male, subset(Data_logCVR_sex_._male, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_sex_._male, file = "models/RobuEstMultiLevelModel_logCVR_strength_sex_._male")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_sex_._male <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sex_._male, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sex_._male, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_sex_._male, file = "models/MultiLevelModel_logCVR_hypertrophy_sex_._male")

### Calculate I^2 
I2_logCVR_hypertrophy_sex_._male <- i2_ml(MultiLevelModel_logCVR_hypertrophy_sex_._male)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male <- robust(MultiLevelModel_logCVR_hypertrophy_sex_._male, subset(Data_logCVR_sex_._male, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male")

### Mean weight
# SMD
# Strength
Data_SMD_strength_weight <- Data_SMD_strength %>%
  filter(!is.na(weight))

MultiLevelModel_SMD_strength_weight <- rma.mv(yi, V=vi, data=Data_SMD_strength_weight,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weight, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_weight, file = "models/MultiLevelModel_SMD_strength_weight")

### Calculate I^2 
I2_SMD_strength_weight <- i2_ml(MultiLevelModel_SMD_strength_weight)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_weight <- robust(MultiLevelModel_SMD_strength_weight, Data_SMD_strength_weight$study)

save(RobuEstMultiLevelModel_SMD_strength_weight, file = "models/RobuEstMultiLevelModel_SMD_strength_weight")

# Hypertrophy
Data_SMD_hypertrophy_weight <- Data_SMD_hypertrophy %>%
  filter(!is.na(weight))

MultiLevelModel_SMD_hypertrophy_weight <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_weight,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weight, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_weight, file = "models/MultiLevelModel_SMD_hypertrophy_weight")

### Calculate I^2 
I2_SMD_hypertrophy_weight <- i2_ml(MultiLevelModel_SMD_hypertrophy_weight)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_weight <- robust(MultiLevelModel_SMD_hypertrophy_weight, Data_SMD_hypertrophy_weight$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_weight, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_weight")

# logCVR
Data_logCVR_weight <- Data_logCVR %>%
  filter(!is.na(weight))

# Strength
MultiLevelModel_logCVR_strength_weight <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weight, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weight, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_weight, file = "models/MultiLevelModel_logCVR_strength_weight")

### Calculate I^2 
I2_logCVR_strength_weight <- i2_ml(MultiLevelModel_logCVR_strength_weight)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_weight <- robust(MultiLevelModel_logCVR_strength_weight, subset(Data_logCVR_weight, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_weight, file = "models/RobuEstMultiLevelModel_logCVR_strength_weight")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_weight <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weight, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weight, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_weight, file = "models/MultiLevelModel_logCVR_hypertrophy_weight")

### Calculate I^2 
I2_logCVR_hypertrophy_weight <- i2_ml(MultiLevelModel_logCVR_hypertrophy_weight)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_weight <- robust(MultiLevelModel_logCVR_hypertrophy_weight, subset(Data_logCVR_weight, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_weight, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_weight")

### Mean BMI
# SMD
# Strength
Data_SMD_strength_bmi <- Data_SMD_strength %>%
  filter(!is.na(bmi))

MultiLevelModel_SMD_strength_bmi <- rma.mv(yi, V=vi, data=Data_SMD_strength_bmi,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ bmi, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_bmi, file = "models/MultiLevelModel_SMD_strength_bmi")

### Calculate I^2 
I2_SMD_strength_bmi <- i2_ml(MultiLevelModel_SMD_strength_bmi)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_bmi <- robust(MultiLevelModel_SMD_strength_bmi, Data_SMD_strength_bmi$study)

save(RobuEstMultiLevelModel_SMD_strength_bmi, file = "models/RobuEstMultiLevelModel_SMD_strength_bmi")

# Hypertrophy
Data_SMD_hypertrophy_bmi <- Data_SMD_hypertrophy %>%
  filter(!is.na(bmi))

MultiLevelModel_SMD_hypertrophy_bmi <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_bmi,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ bmi, method="REML", test="t",
                                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
                                             )

save(MultiLevelModel_SMD_hypertrophy_bmi, file = "models/MultiLevelModel_SMD_hypertrophy_bmi")

### Calculate I^2 
I2_SMD_hypertrophy_bmi <- i2_ml(MultiLevelModel_SMD_hypertrophy_bmi)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_bmi <- robust(MultiLevelModel_SMD_hypertrophy_bmi, Data_SMD_hypertrophy_bmi$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_bmi, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_bmi")

# logCVR
Data_logCVR_bmi <- Data_logCVR %>%
  filter(!is.na(bmi))

# Strength
MultiLevelModel_logCVR_strength_bmi <- rma.mv(yi, V=vi, data=subset(Data_logCVR_bmi, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ bmi, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_bmi, file = "models/MultiLevelModel_logCVR_strength_bmi")

### Calculate I^2 
I2_logCVR_strength_bmi <- i2_ml(MultiLevelModel_logCVR_strength_bmi)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_bmi <- robust(MultiLevelModel_logCVR_strength_bmi, subset(Data_logCVR_bmi, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_bmi, file = "models/RobuEstMultiLevelModel_logCVR_strength_bmi")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_bmi <- rma.mv(yi, V=vi, data=subset(Data_logCVR_bmi, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ bmi, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_bmi, file = "models/MultiLevelModel_logCVR_hypertrophy_bmi")

### Calculate I^2 
I2_logCVR_hypertrophy_bmi <- i2_ml(MultiLevelModel_logCVR_hypertrophy_bmi)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_bmi <- robust(MultiLevelModel_logCVR_hypertrophy_bmi, subset(Data_logCVR_bmi, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_bmi")

### Trained vs untrained
# SMD
# Strength
Data_SMD_strength_train_status <- Data_SMD_strength %>%
  filter(!is.na(train_status))

MultiLevelModel_SMD_strength_train_status <- rma.mv(yi, V=vi, data=Data_SMD_strength_train_status,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ train_status - 1, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_train_status, file = "models/MultiLevelModel_SMD_strength_train_status")

### Calculate I^2 
I2_SMD_strength_train_status <- i2_ml(MultiLevelModel_SMD_strength_train_status)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_train_status <- robust(MultiLevelModel_SMD_strength_train_status, Data_SMD_strength_train_status$study)

save(RobuEstMultiLevelModel_SMD_strength_train_status, file = "models/RobuEstMultiLevelModel_SMD_strength_train_status")

# Hypertrophy
Data_SMD_hypertrophy_train_status <- Data_SMD_hypertrophy %>%
  filter(!is.na(train_status))

MultiLevelModel_SMD_hypertrophy_train_status <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_train_status,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ train_status - 1, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_train_status, file = "models/MultiLevelModel_SMD_hypertrophy_train_status")

### Calculate I^2 
I2_SMD_hypertrophy_train_status <- i2_ml(MultiLevelModel_SMD_hypertrophy_train_status)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_train_status <- robust(MultiLevelModel_SMD_hypertrophy_train_status, Data_SMD_hypertrophy_train_status$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_train_status, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_train_status")

# logCVR
Data_logCVR_train_status <- Data_logCVR %>%
  filter(!is.na(train_status))

# Strength
MultiLevelModel_logCVR_strength_train_status <- rma.mv(yi, V=vi, data=subset(Data_logCVR_train_status, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ train_status - 1, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_train_status, file = "models/MultiLevelModel_logCVR_strength_train_status")

### Calculate I^2 
I2_logCVR_strength_train_status <- i2_ml(MultiLevelModel_logCVR_strength_train_status)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_train_status <- robust(MultiLevelModel_logCVR_strength_train_status, subset(Data_logCVR_train_status, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_train_status, file = "models/RobuEstMultiLevelModel_logCVR_strength_train_status")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_train_status <- rma.mv(yi, V=vi, data=subset(Data_logCVR_train_status, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ train_status - 1, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_train_status, file = "models/MultiLevelModel_logCVR_hypertrophy_train_status")

### Calculate I^2 
I2_logCVR_hypertrophy_train_status <- i2_ml(MultiLevelModel_logCVR_hypertrophy_train_status)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_train_status <- robust(MultiLevelModel_logCVR_hypertrophy_train_status, subset(Data_logCVR_train_status, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_train_status")

### Healthy vs clinical sample
# SMD
# Strength
Data_SMD_strength_healthy_clinical <- Data_SMD_strength %>%
  filter(!is.na(healthy_clinical))

MultiLevelModel_SMD_strength_healthy_clinical <- rma.mv(yi, V=vi, data=Data_SMD_strength_healthy_clinical,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ healthy_clinical - 1, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_healthy_clinical, file = "models/MultiLevelModel_SMD_strength_healthy_clinical")

### Calculate I^2 
I2_SMD_strength_healthy_clinical <- i2_ml(MultiLevelModel_SMD_strength_healthy_clinical)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_healthy_clinical <- robust(MultiLevelModel_SMD_strength_healthy_clinical, Data_SMD_strength_healthy_clinical$study)

save(RobuEstMultiLevelModel_SMD_strength_healthy_clinical, file = "models/RobuEstMultiLevelModel_SMD_strength_healthy_clinical")

# Hypertrophy
Data_SMD_hypertrophy_healthy_clinical <- Data_SMD_hypertrophy %>%
  filter(!is.na(healthy_clinical))

MultiLevelModel_SMD_hypertrophy_healthy_clinical <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_healthy_clinical,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ healthy_clinical - 1, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_healthy_clinical, file = "models/MultiLevelModel_SMD_hypertrophy_healthy_clinical")

### Calculate I^2 
I2_SMD_hypertrophy_healthy_clinical <- i2_ml(MultiLevelModel_SMD_hypertrophy_healthy_clinical)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical <- robust(MultiLevelModel_SMD_hypertrophy_healthy_clinical, Data_SMD_hypertrophy_healthy_clinical$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical")

# logCVR
Data_logCVR_healthy_clinical <- Data_logCVR %>%
  filter(!is.na(healthy_clinical))

# Strength
MultiLevelModel_logCVR_strength_healthy_clinical <- rma.mv(yi, V=vi, data=subset(Data_logCVR_healthy_clinical, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ healthy_clinical - 1, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_healthy_clinical, file = "models/MultiLevelModel_logCVR_strength_healthy_clinical")

### Calculate I^2 
I2_logCVR_strength_healthy_clinical <- i2_ml(MultiLevelModel_logCVR_strength_healthy_clinical)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_healthy_clinical <- robust(MultiLevelModel_logCVR_strength_healthy_clinical, subset(Data_logCVR_healthy_clinical, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical, file = "models/RobuEstMultiLevelModel_logCVR_strength_healthy_clinical")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_healthy_clinical <- rma.mv(yi, V=vi, data=subset(Data_logCVR_healthy_clinical, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ healthy_clinical - 1, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_healthy_clinical, file = "models/MultiLevelModel_logCVR_hypertrophy_healthy_clinical")

### Calculate I^2 
I2_logCVR_hypertrophy_healthy_clinical <- i2_ml(MultiLevelModel_logCVR_hypertrophy_healthy_clinical)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical <- robust(MultiLevelModel_logCVR_hypertrophy_healthy_clinical, subset(Data_logCVR_healthy_clinical, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical")

### Only RT intervention (y), or additional e.g., supplements, aerobic (n)
# SMD
# Strength
Data_SMD_strength_RT_only <- Data_SMD_strength %>%
  filter(!is.na(RT_only))

MultiLevelModel_SMD_strength_RT_only <- rma.mv(yi, V=vi, data=Data_SMD_strength_RT_only,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ RT_only - 1, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_RT_only, file = "models/MultiLevelModel_SMD_strength_RT_only")

### Calculate I^2 
I2_SMD_strength_RT_only <- i2_ml(MultiLevelModel_SMD_strength_RT_only)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_RT_only <- robust(MultiLevelModel_SMD_strength_RT_only, Data_SMD_strength_RT_only$study)

save(RobuEstMultiLevelModel_SMD_strength_RT_only, file = "models/RobuEstMultiLevelModel_SMD_strength_RT_only")

# Hypertrophy
Data_SMD_hypertrophy_RT_only <- Data_SMD_hypertrophy %>%
  filter(!is.na(RT_only))

MultiLevelModel_SMD_hypertrophy_RT_only <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_RT_only,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ RT_only - 1, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_RT_only, file = "models/MultiLevelModel_SMD_hypertrophy_RT_only")

### Calculate I^2 
I2_SMD_hypertrophy_RT_only <- i2_ml(MultiLevelModel_SMD_hypertrophy_RT_only)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_RT_only <- robust(MultiLevelModel_SMD_hypertrophy_RT_only, Data_SMD_hypertrophy_RT_only$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_RT_only, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_RT_only")

# logCVR
Data_logCVR_RT_only <- Data_logCVR %>%
  filter(!is.na(RT_only))

# Strength
MultiLevelModel_logCVR_strength_RT_only <- rma.mv(yi, V=vi, data=subset(Data_logCVR_RT_only, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ RT_only - 1, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_RT_only, file = "models/MultiLevelModel_logCVR_strength_RT_only")

### Calculate I^2 
I2_logCVR_strength_RT_only <- i2_ml(MultiLevelModel_logCVR_strength_RT_only)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_RT_only <- robust(MultiLevelModel_logCVR_strength_RT_only, subset(Data_logCVR_RT_only, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_RT_only, file = "models/RobuEstMultiLevelModel_logCVR_strength_RT_only")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_RT_only <- rma.mv(yi, V=vi, data=subset(Data_logCVR_RT_only, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ RT_only - 1, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_RT_only, file = "models/MultiLevelModel_logCVR_hypertrophy_RT_only")

### Calculate I^2 
I2_logCVR_hypertrophy_RT_only <- i2_ml(MultiLevelModel_logCVR_hypertrophy_RT_only)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only <- robust(MultiLevelModel_logCVR_hypertrophy_RT_only, subset(Data_logCVR_RT_only, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only")

### Duration of intervention in weeks
# SMD
# Strength
Data_SMD_strength_weeks <- Data_SMD_strength %>%
  filter(!is.na(weeks))

MultiLevelModel_SMD_strength_weeks <- rma.mv(yi, V=vi, data=Data_SMD_strength_weeks,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weeks, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_weeks, file = "models/MultiLevelModel_SMD_strength_weeks")

### Calculate I^2 
I2_SMD_strength_weeks <- i2_ml(MultiLevelModel_SMD_strength_weeks)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_weeks <- robust(MultiLevelModel_SMD_strength_weeks, Data_SMD_strength_weeks$study)

save(RobuEstMultiLevelModel_SMD_strength_weeks, file = "models/RobuEstMultiLevelModel_SMD_strength_weeks")

# Hypertrophy
Data_SMD_hypertrophy_weeks <- Data_SMD_hypertrophy %>%
  filter(!is.na(weeks))

MultiLevelModel_SMD_hypertrophy_weeks <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_weeks,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weeks, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_weeks, file = "models/MultiLevelModel_SMD_hypertrophy_weeks")

### Calculate I^2 
I2_SMD_hypertrophy_weeks <- i2_ml(MultiLevelModel_SMD_hypertrophy_weeks)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_weeks <- robust(MultiLevelModel_SMD_hypertrophy_weeks, Data_SMD_hypertrophy_weeks$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_weeks, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_weeks")

# logCVR
Data_logCVR_weeks <- Data_logCVR %>%
  filter(!is.na(weeks))

# Strength
MultiLevelModel_logCVR_strength_weeks <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weeks, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weeks, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_weeks, file = "models/MultiLevelModel_logCVR_strength_weeks")

### Calculate I^2 
I2_logCVR_strength_weeks <- i2_ml(MultiLevelModel_logCVR_strength_weeks)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_weeks <- robust(MultiLevelModel_logCVR_strength_weeks, subset(Data_logCVR_weeks, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_weeks, file = "models/RobuEstMultiLevelModel_logCVR_strength_weeks")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_weeks <- rma.mv(yi, V=vi, data=subset(Data_logCVR_weeks, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ weeks, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_weeks, file = "models/MultiLevelModel_logCVR_hypertrophy_weeks")

### Calculate I^2 
I2_logCVR_hypertrophy_weeks <- i2_ml(MultiLevelModel_logCVR_hypertrophy_weeks)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_weeks <- robust(MultiLevelModel_logCVR_hypertrophy_weeks, subset(Data_logCVR_weeks, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_weeks")

### Weekly frequency in days
# SMD
# Strength
Data_SMD_strength_freq <- Data_SMD_strength %>%
  filter(!is.na(freq))

MultiLevelModel_SMD_strength_freq <- rma.mv(yi, V=vi, data=Data_SMD_strength_freq,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ freq, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_freq, file = "models/MultiLevelModel_SMD_strength_freq")

### Calculate I^2 
I2_SMD_strength_freq <- i2_ml(MultiLevelModel_SMD_strength_freq)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_freq <- robust(MultiLevelModel_SMD_strength_freq, Data_SMD_strength_freq$study)

save(RobuEstMultiLevelModel_SMD_strength_freq, file = "models/RobuEstMultiLevelModel_SMD_strength_freq")

# Hypertrophy
Data_SMD_hypertrophy_freq <- Data_SMD_hypertrophy %>%
  filter(!is.na(freq))

MultiLevelModel_SMD_hypertrophy_freq <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_freq,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ freq, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_freq, file = "models/MultiLevelModel_SMD_hypertrophy_freq")

### Calculate I^2 
I2_SMD_hypertrophy_freq <- i2_ml(MultiLevelModel_SMD_hypertrophy_freq)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_freq <- robust(MultiLevelModel_SMD_hypertrophy_freq, Data_SMD_hypertrophy_freq$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_freq, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_freq")

# logCVR
Data_logCVR_freq <- Data_logCVR %>%
  filter(!is.na(freq))

# Strength
MultiLevelModel_logCVR_strength_freq <- rma.mv(yi, V=vi, data=subset(Data_logCVR_freq, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ freq, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_freq, file = "models/MultiLevelModel_logCVR_strength_freq")

### Calculate I^2 
I2_logCVR_strength_freq <- i2_ml(MultiLevelModel_logCVR_strength_freq)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_freq <- robust(MultiLevelModel_logCVR_strength_freq, subset(Data_logCVR_freq, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_freq, file = "models/RobuEstMultiLevelModel_logCVR_strength_freq")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_freq <- rma.mv(yi, V=vi, data=subset(Data_logCVR_freq, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ freq, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_freq, file = "models/MultiLevelModel_logCVR_hypertrophy_freq")

### Calculate I^2 
I2_logCVR_hypertrophy_freq <- i2_ml(MultiLevelModel_logCVR_hypertrophy_freq)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_freq <- robust(MultiLevelModel_logCVR_hypertrophy_freq, subset(Data_logCVR_freq, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_freq, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_freq")

### Number of exercises per workout
# SMD
# Strength
Data_SMD_strength_exercises <- Data_SMD_strength %>%
  filter(!is.na(exercises))

MultiLevelModel_SMD_strength_exercises <- rma.mv(yi, V=vi, data=Data_SMD_strength_exercises,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ exercises, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_exercises, file = "models/MultiLevelModel_SMD_strength_exercises")

### Calculate I^2 
I2_SMD_strength_exercises <- i2_ml(MultiLevelModel_SMD_strength_exercises)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_exercises <- robust(MultiLevelModel_SMD_strength_exercises, Data_SMD_strength_exercises$study)

save(RobuEstMultiLevelModel_SMD_strength_exercises, file = "models/RobuEstMultiLevelModel_SMD_strength_exercises")

# Hypertrophy
Data_SMD_hypertrophy_exercises <- Data_SMD_hypertrophy %>%
  filter(!is.na(exercises))

MultiLevelModel_SMD_hypertrophy_exercises <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_exercises,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ exercises, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_exercises, file = "models/MultiLevelModel_SMD_hypertrophy_exercises")

### Calculate I^2 
I2_SMD_hypertrophy_exercises <- i2_ml(MultiLevelModel_SMD_hypertrophy_exercises)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_exercises <- robust(MultiLevelModel_SMD_hypertrophy_exercises, Data_SMD_hypertrophy_exercises$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_exercises, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_exercises")

# logCVR
Data_logCVR_exercises <- Data_logCVR %>%
  filter(!is.na(exercises))

# Strength
MultiLevelModel_logCVR_strength_exercises <- rma.mv(yi, V=vi, data=subset(Data_logCVR_exercises, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ exercises, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_exercises, file = "models/MultiLevelModel_logCVR_strength_exercises")

### Calculate I^2 
I2_logCVR_strength_exercises <- i2_ml(MultiLevelModel_logCVR_strength_exercises)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_exercises <- robust(MultiLevelModel_logCVR_strength_exercises, subset(Data_logCVR_exercises, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_exercises, file = "models/RobuEstMultiLevelModel_logCVR_strength_exercises")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_exercises <- rma.mv(yi, V=vi, data=subset(Data_logCVR_exercises, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ exercises, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_exercises, file = "models/MultiLevelModel_logCVR_hypertrophy_exercises")

### Calculate I^2 
I2_logCVR_hypertrophy_exercises <- i2_ml(MultiLevelModel_logCVR_hypertrophy_exercises)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_exercises <- robust(MultiLevelModel_logCVR_hypertrophy_exercises, subset(Data_logCVR_exercises, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_exercises")

### Number of sets per exercise
# SMD
# Strength
Data_SMD_strength_sets_exercise <- Data_SMD_strength %>%
  filter(!is.na(sets_exercise))

MultiLevelModel_SMD_strength_sets_exercise <- rma.mv(yi, V=vi, data=Data_SMD_strength_sets_exercise,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sets_exercise, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_sets_exercise, file = "models/MultiLevelModel_SMD_strength_sets_exercise")

### Calculate I^2 
I2_SMD_strength_sets_exercise <- i2_ml(MultiLevelModel_SMD_strength_sets_exercise)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_sets_exercise <- robust(MultiLevelModel_SMD_strength_sets_exercise, Data_SMD_strength_sets_exercise$study)

save(RobuEstMultiLevelModel_SMD_strength_sets_exercise, file = "models/RobuEstMultiLevelModel_SMD_strength_sets_exercise")

# Hypertrophy
Data_SMD_hypertrophy_sets_exercise <- Data_SMD_hypertrophy %>%
  filter(!is.na(sets_exercise))

MultiLevelModel_SMD_hypertrophy_sets_exercise <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_sets_exercise,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sets_exercise, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_sets_exercise, file = "models/MultiLevelModel_SMD_hypertrophy_sets_exercise")

### Calculate I^2 
I2_SMD_hypertrophy_sets_exercise <- i2_ml(MultiLevelModel_SMD_hypertrophy_sets_exercise)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise <- robust(MultiLevelModel_SMD_hypertrophy_sets_exercise, Data_SMD_hypertrophy_sets_exercise$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise")

# logCVR
Data_logCVR_sets_exercise <- Data_logCVR %>%
  filter(!is.na(sets_exercise))

# Strength
MultiLevelModel_logCVR_strength_sets_exercise <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sets_exercise, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sets_exercise, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_sets_exercise, file = "models/MultiLevelModel_logCVR_strength_sets_exercise")

### Calculate I^2 
I2_logCVR_strength_sets_exercise <- i2_ml(MultiLevelModel_logCVR_strength_sets_exercise)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_sets_exercise <- robust(MultiLevelModel_logCVR_strength_sets_exercise, subset(Data_logCVR_sets_exercise, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_sets_exercise, file = "models/RobuEstMultiLevelModel_logCVR_strength_sets_exercise")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_sets_exercise <- rma.mv(yi, V=vi, data=subset(Data_logCVR_sets_exercise, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ sets_exercise, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_sets_exercise, file = "models/MultiLevelModel_logCVR_hypertrophy_sets_exercise")

### Calculate I^2 
I2_logCVR_hypertrophy_sets_exercise <- i2_ml(MultiLevelModel_logCVR_hypertrophy_sets_exercise)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise <- robust(MultiLevelModel_logCVR_hypertrophy_sets_exercise, subset(Data_logCVR_sets_exercise, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise")

### Repetitions per exercise
# SMD
# Strength
Data_SMD_strength_reps <- Data_SMD_strength %>%
  filter(!is.na(reps))

MultiLevelModel_SMD_strength_reps <- rma.mv(yi, V=vi, data=Data_SMD_strength_reps,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ reps, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_reps, file = "models/MultiLevelModel_SMD_strength_reps")

### Calculate I^2 
I2_SMD_strength_reps <- i2_ml(MultiLevelModel_SMD_strength_reps)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_reps <- robust(MultiLevelModel_SMD_strength_reps, Data_SMD_strength_reps$study)

save(RobuEstMultiLevelModel_SMD_strength_reps, file = "models/RobuEstMultiLevelModel_SMD_strength_reps")

# Hypertrophy
Data_SMD_hypertrophy_reps <- Data_SMD_hypertrophy %>%
  filter(!is.na(reps))

MultiLevelModel_SMD_hypertrophy_reps <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_reps,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ reps, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_reps, file = "models/MultiLevelModel_SMD_hypertrophy_reps")

### Calculate I^2 
I2_SMD_hypertrophy_reps <- i2_ml(MultiLevelModel_SMD_hypertrophy_reps)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_reps <- robust(MultiLevelModel_SMD_hypertrophy_reps, Data_SMD_hypertrophy_reps$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_reps, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_reps")

# logCVR
Data_logCVR_reps <- Data_logCVR %>%
  filter(!is.na(reps))

# Strength
MultiLevelModel_logCVR_strength_reps <- rma.mv(yi, V=vi, data=subset(Data_logCVR_reps, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ reps, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_reps, file = "models/MultiLevelModel_logCVR_strength_reps")

### Calculate I^2 
I2_logCVR_strength_reps <- i2_ml(MultiLevelModel_logCVR_strength_reps)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_reps <- robust(MultiLevelModel_logCVR_strength_reps, subset(Data_logCVR_reps, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_reps, file = "models/RobuEstMultiLevelModel_logCVR_strength_reps")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_reps <- rma.mv(yi, V=vi, data=subset(Data_logCVR_reps, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ reps, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_reps, file = "models/MultiLevelModel_logCVR_hypertrophy_reps")

### Calculate I^2 
I2_logCVR_hypertrophy_reps <- i2_ml(MultiLevelModel_logCVR_hypertrophy_reps)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_reps <- robust(MultiLevelModel_logCVR_hypertrophy_reps, subset(Data_logCVR_reps, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_reps, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_reps")

### Relative load (%1RM)
# SMD
# Strength
Data_SMD_strength_load <- Data_SMD_strength %>%
  filter(!is.na(load))

MultiLevelModel_SMD_strength_load <- rma.mv(yi, V=vi, data=Data_SMD_strength_load,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ load, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_load, file = "models/MultiLevelModel_SMD_strength_load")

### Calculate I^2 
I2_SMD_strength_load <- i2_ml(MultiLevelModel_SMD_strength_load)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_load <- robust(MultiLevelModel_SMD_strength_load, Data_SMD_strength_load$study)

save(RobuEstMultiLevelModel_SMD_strength_load, file = "models/RobuEstMultiLevelModel_SMD_strength_load")

# Hypertrophy
Data_SMD_hypertrophy_load <- Data_SMD_hypertrophy %>%
  filter(!is.na(load))

MultiLevelModel_SMD_hypertrophy_load <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_load,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ load, method="REML", test="t",
                                             # control=list(optimizer="optim", optmethod="Nelder-Mead")
                                             )

save(MultiLevelModel_SMD_hypertrophy_load, file = "models/MultiLevelModel_SMD_hypertrophy_load")

### Calculate I^2 
I2_SMD_hypertrophy_load <- i2_ml(MultiLevelModel_SMD_hypertrophy_load)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_load <- robust(MultiLevelModel_SMD_hypertrophy_load, Data_SMD_hypertrophy_load$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_load, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_load")

# logCVR
Data_logCVR_load <- Data_logCVR %>%
  filter(!is.na(load))

# Strength
MultiLevelModel_logCVR_strength_load <- rma.mv(yi, V=vi, data=subset(Data_logCVR_load, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ load, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_load, file = "models/MultiLevelModel_logCVR_strength_load")

### Calculate I^2 
I2_logCVR_strength_load <- i2_ml(MultiLevelModel_logCVR_strength_load)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_load <- robust(MultiLevelModel_logCVR_strength_load, subset(Data_logCVR_load, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_load, file = "models/RobuEstMultiLevelModel_logCVR_strength_load")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_load <- rma.mv(yi, V=vi, data=subset(Data_logCVR_load, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ load, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_load, file = "models/MultiLevelModel_logCVR_hypertrophy_load")

### Calculate I^2 
I2_logCVR_hypertrophy_load <- i2_ml(MultiLevelModel_logCVR_hypertrophy_load)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_load <- robust(MultiLevelModel_logCVR_hypertrophy_load, subset(Data_logCVR_load, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_load, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_load")

### Trained to task failure?
# SMD
# Strength
Data_SMD_strength_task_failure_y_n <- Data_SMD_strength %>%
  filter(!is.na(task_failure_y_n))

MultiLevelModel_SMD_strength_task_failure_y_n <- rma.mv(yi, V=vi, data=Data_SMD_strength_task_failure_y_n,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ task_failure_y_n - 1, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_task_failure_y_n, file = "models/MultiLevelModel_SMD_strength_task_failure_y_n")

### Calculate I^2 
I2_SMD_strength_task_failure_y_n <- i2_ml(MultiLevelModel_SMD_strength_task_failure_y_n)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_task_failure_y_n <- robust(MultiLevelModel_SMD_strength_task_failure_y_n, Data_SMD_strength_task_failure_y_n$study)

save(RobuEstMultiLevelModel_SMD_strength_task_failure_y_n, file = "models/RobuEstMultiLevelModel_SMD_strength_task_failure_y_n")

# Hypertrophy
Data_SMD_hypertrophy_task_failure_y_n <- Data_SMD_hypertrophy %>%
  filter(!is.na(task_failure_y_n))

MultiLevelModel_SMD_hypertrophy_task_failure_y_n <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_task_failure_y_n,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ task_failure_y_n - 1, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_task_failure_y_n, file = "models/MultiLevelModel_SMD_hypertrophy_task_failure_y_n")

### Calculate I^2 
I2_SMD_hypertrophy_task_failure_y_n <- i2_ml(MultiLevelModel_SMD_hypertrophy_task_failure_y_n)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n <- robust(MultiLevelModel_SMD_hypertrophy_task_failure_y_n, Data_SMD_hypertrophy_task_failure_y_n$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n")

# logCVR
Data_logCVR_task_failure_y_n <- Data_logCVR %>%
  filter(!is.na(task_failure_y_n))

# Strength
MultiLevelModel_logCVR_strength_task_failure_y_n <- rma.mv(yi, V=vi, data=subset(Data_logCVR_task_failure_y_n, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ task_failure_y_n - 1, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_task_failure_y_n, file = "models/MultiLevelModel_logCVR_strength_task_failure_y_n")

### Calculate I^2 
I2_logCVR_strength_task_failure_y_n <- i2_ml(MultiLevelModel_logCVR_strength_task_failure_y_n)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n <- robust(MultiLevelModel_logCVR_strength_task_failure_y_n, subset(Data_logCVR_task_failure_y_n, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n, file = "models/RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_task_failure_y_n <- rma.mv(yi, V=vi, data=subset(Data_logCVR_task_failure_y_n, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ task_failure_y_n - 1, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n, file = "models/MultiLevelModel_logCVR_hypertrophy_task_failure_y_n")

### Calculate I^2 
I2_logCVR_hypertrophy_task_failure_y_n <- i2_ml(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n <- robust(MultiLevelModel_logCVR_hypertrophy_task_failure_y_n, subset(Data_logCVR_task_failure_y_n, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n")

### Outcome measurement type
# SMD
# Strength
Data_SMD_strength_measure <- Data_SMD_strength %>%
  filter(!is.na(measure))

MultiLevelModel_SMD_strength_measure <- rma.mv(yi, V=vi, data=Data_SMD_strength_measure,
                                          slab=paste(label),
                                          random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ measure - 1, method="REML", test="t",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_strength_measure, file = "models/MultiLevelModel_SMD_strength_measure")

### Calculate I^2 
I2_SMD_strength_measure <- i2_ml(MultiLevelModel_SMD_strength_measure)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_strength_measure <- robust(MultiLevelModel_SMD_strength_measure, Data_SMD_strength_measure$study)

save(RobuEstMultiLevelModel_SMD_strength_measure, file = "models/RobuEstMultiLevelModel_SMD_strength_measure")

# Hypertrophy
Data_SMD_hypertrophy_measure <- Data_SMD_hypertrophy %>%
  filter(!is.na(measure))

MultiLevelModel_SMD_hypertrophy_measure <- rma.mv(yi, V=vi, data=Data_SMD_hypertrophy_measure,
                                             slab=paste(label),
                                             random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ measure - 1, method="REML", test="t",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_SMD_hypertrophy_measure, file = "models/MultiLevelModel_SMD_hypertrophy_measure")

### Calculate I^2 
I2_SMD_hypertrophy_measure <- i2_ml(MultiLevelModel_SMD_hypertrophy_measure)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_SMD_hypertrophy_measure <- robust(MultiLevelModel_SMD_hypertrophy_measure, Data_SMD_hypertrophy_measure$study)

save(RobuEstMultiLevelModel_SMD_hypertrophy_measure, file = "models/RobuEstMultiLevelModel_SMD_hypertrophy_measure")

# logCVR
Data_logCVR_measure <- Data_logCVR %>%
  filter(!is.na(measure))

# Strength
MultiLevelModel_logCVR_strength_measure <- rma.mv(yi, V=vi, data=subset(Data_logCVR_measure, outcome == "strength"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ measure - 1, method="REML", test="t",
                                                 control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_strength_measure, file = "models/MultiLevelModel_logCVR_strength_measure")

### Calculate I^2 
I2_logCVR_strength_measure <- i2_ml(MultiLevelModel_logCVR_strength_measure)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength_measure <- robust(MultiLevelModel_logCVR_strength_measure, subset(Data_logCVR_measure, outcome == "strength")$study)

save(RobuEstMultiLevelModel_logCVR_strength_measure, file = "models/RobuEstMultiLevelModel_logCVR_strength_measure")

# Hypertrophy
MultiLevelModel_logCVR_hypertrophy_measure <- rma.mv(yi, V=vi, data=subset(Data_logCVR_measure, outcome == "hypertrophy"),
                                                    slab=paste(label),
                                                    random = list(~ 1 | study, ~ 1 | arm, ~ 1 | es), mods = ~ measure - 1, method="REML", test="t",
                                                    control=list(optimizer="optim", optmethod="Nelder-Mead"))

save(MultiLevelModel_logCVR_hypertrophy_measure, file = "models/MultiLevelModel_logCVR_hypertrophy_measure")

### Calculate I^2 
I2_logCVR_hypertrophy_measure <- i2_ml(MultiLevelModel_logCVR_hypertrophy_measure)

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy_measure <- robust(MultiLevelModel_logCVR_hypertrophy_measure, subset(Data_logCVR_measure, outcome == "hypertrophy")$study)

save(RobuEstMultiLevelModel_logCVR_hypertrophy_measure, file = "models/RobuEstMultiLevelModel_logCVR_hypertrophy_measure")

### Collate all moderator estimates for log CVR and SMD

mods_SMD_logCVR <- rbind(data.frame(Outcome = "strength",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_strength$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_strength, 2)[2],
                                    `I^2 arm` = round(I2_SMD_strength, 2)[3],
                                    `I^2 effect` = round(I2_SMD_strength, 2)[4]),
                         data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "TESTEX score",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_TESTEX$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_TESTEX$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_TESTEX$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_strength_TESTEX, 2)[2],
                                    `I^2 arm` = round(I2_SMD_strength_TESTEX, 2)[3],
                                    `I^2 effect` = round(I2_SMD_strength_TESTEX, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Age",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_age$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_age$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_age$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_age, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_age, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_age, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Proportion Male",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_sex_._male$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_sex_._male$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_sex_._male$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_sex_._male, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_sex_._male, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_sex_._male, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Weight",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_weight$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_weight$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_weight$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_weight, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_weight, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_weight, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "BMI",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_bmi$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_bmi$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_bmi$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_bmi, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_bmi, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_bmi, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_train_status$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_train_status$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_train_status$ci.ub, 2),
                                      `I^2 study` = round(I2_SMD_strength_train_status, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_train_status, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_train_status, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Healthy Sample", "Clinical Sample"),
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_healthy_clinical$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_healthy_clinical$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_healthy_clinical$ci.ub, 2),
                                      `I^2 study` = round(I2_SMD_strength_healthy_clinical, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_healthy_clinical, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_healthy_clinical, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_RT_only$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_RT_only$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_RT_only$ci.ub, 2),
                                      `I^2 study` = round(I2_SMD_strength_RT_only, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_RT_only, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_RT_only, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Duration (weeks)",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_weeks$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_weeks$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_weeks$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_weeks, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_weeks, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_weeks, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Weekly Frequency",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_freq$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_freq$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_freq$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_freq, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_freq, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_freq, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Number of Exercises",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_exercises$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_exercises$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_exercises$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_exercises, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_exercises, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_exercises, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Sets per Exercise",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_sets_exercise$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_sets_exercise$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_sets_exercise$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_sets_exercise, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_sets_exercise, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_sets_exercise, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Number of Repetitions",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_reps$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_reps$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_reps$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_reps, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_reps, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_reps, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = "Load (%1RM)",
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_load$b, 2)[2],
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_load$ci.lb, 2)[2],
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_load$ci.ub, 2)[2],
                                      `I^2 study` = round(I2_SMD_strength_load, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_load, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_load, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_task_failure_y_n$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_task_failure_y_n$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_task_failure_y_n$ci.ub, 2),
                                      `I^2 study` = round(I2_SMD_strength_task_failure_y_n, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_task_failure_y_n, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_task_failure_y_n, 2)[4]),
                           data.frame(Outcome = "strength",
                                      Model = "SMD",
                                      Moderator = c("Outcome Measure (12RM)","Outcome Measure (1RM)","Outcome Measure (3RM)","Outcome Measure (5RM)","Outcome Measure (6RM)","Outcome Measure (Isokinetic)", "Outcome Measure (Isometric)"),
                                      Estimate = round(RobuEstMultiLevelModel_SMD_strength_measure$b, 2),
                                      Lower = round(RobuEstMultiLevelModel_SMD_strength_measure$ci.lb, 2),
                                      Upper = round(RobuEstMultiLevelModel_SMD_strength_measure$ci.ub, 2),
                                      `I^2 study` = round(I2_SMD_strength_measure, 2)[2],
                                      `I^2 arm` = round(I2_SMD_strength_measure, 2)[3],
                                      `I^2 effect` = round(I2_SMD_strength_measure, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_TESTEX$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_TESTEX, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_TESTEX, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_TESTEX, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_age$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_age, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_age, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_age, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_sex_._male$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_sex_._male, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_sex_._male, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_sex_._male, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_weight$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_weight, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_weight, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_weight, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_bmi$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_bmi, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_bmi, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_bmi, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_train_status$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength_train_status, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_train_status, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_train_status, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_healthy_clinical$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength_healthy_clinical, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_healthy_clinical, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_healthy_clinical, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_RT_only$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength_RT_only, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_RT_only, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_RT_only, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_weeks$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_weeks, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_weeks, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_weeks, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_freq$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_freq, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_freq, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_freq, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_exercises$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_exercises, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_exercises, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_exercises, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_sets_exercise$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_sets_exercise, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_sets_exercise, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_sets_exercise, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_reps$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_reps, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_reps, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_reps, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_load$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_strength_load, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_load, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_load, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_task_failure_y_n$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength_task_failure_y_n, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_task_failure_y_n, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_task_failure_y_n, 2)[4]),
                         data.frame(Outcome = "strength",
                                    Model = "logCVR",
                                    Moderator = c("Outcome Measure (12RM)","Outcome Measure (1RM)","Outcome Measure (3RM)","Outcome Measure (5RM)","Outcome Measure (6RM)","Outcome Measure (Isokinetic)", "Outcome Measure (Isometric)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_strength_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_strength_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_strength_measure$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_strength_measure, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_strength_measure, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_strength_measure, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_TESTEX$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_TESTEX, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_TESTEX, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_TESTEX, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_age$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_age, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_age, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_age, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_sex_._male$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_sex_._male, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_sex_._male, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_sex_._male, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_weight$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_weight, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_weight, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_weight, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_bmi$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_bmi, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_bmi, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_bmi, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_train_status$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy_train_status, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_train_status, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_train_status, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_healthy_clinical$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy_healthy_clinical, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_healthy_clinical, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_healthy_clinical, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_RT_only$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy_RT_only, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_RT_only, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_RT_only, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_weeks$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_weeks, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_weeks, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_weeks, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_freq$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_freq, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_freq, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_freq, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_exercises$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_exercises, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_exercises, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_exercises, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_sets_exercise$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_sets_exercise, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_sets_exercise, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_sets_exercise, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_reps$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_reps, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_reps, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_reps, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_load$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_SMD_hypertrophy_load, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_load, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_load, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_task_failure_y_n$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy_task_failure_y_n, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_task_failure_y_n, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_task_failure_y_n, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "SMD",
                                    Moderator = c("Outcome Measure (BIA)","Outcome Measure (Biopsy: Type i)","Outcome Measure (Biopsy: Type ii)","Outcome Measure (Biopsy: Type iia)","Outcome Measure (Biopsy: Type iib)","Outcome Measure (BodPod)", "Outcome Measure (Circumference)","Outcome Measure (CT)","Outcome Measure (DXA)","Outcome Measure (Hydrostatic Weighing)","Outcome Measure (MRI)","Outcome Measure (Skinfold)","Outcome Measure (US)"),
                                    Estimate = round(RobuEstMultiLevelModel_SMD_hypertrophy_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_SMD_hypertrophy_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_SMD_hypertrophy_measure$ci.ub, 2),
                                    `I^2 study` = round(I2_SMD_hypertrophy_measure, 2)[2],
                                    `I^2 arm` = round(I2_SMD_hypertrophy_measure, 2)[3],
                                    `I^2 effect` = round(I2_SMD_hypertrophy_measure, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Main model",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "TESTEX score",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_TESTEX$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_TESTEX, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_TESTEX, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_TESTEX, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Age",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_age$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_age, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_age, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_age, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Proportion Male",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sex_._male$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_sex_._male, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_sex_._male, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_sex_._male, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Weight",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weight$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_weight, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_weight, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_weight, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "BMI",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_bmi$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_bmi, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_bmi, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_bmi, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Training Status (trained)","Training Status (untrained)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_train_status$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy_train_status, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_train_status, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_train_status, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Healthy Sample", "Clinical Sample"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_healthy_clinical$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy_healthy_clinical, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_healthy_clinical, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_healthy_clinical, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("RT + Adjuvant Intervention", "RT Only Intervention"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_RT_only$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy_RT_only, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_RT_only, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_RT_only, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Duration (weeks)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_weeks$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_weeks, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_weeks, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_weeks, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Weekly Frequency",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_freq$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_freq, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_freq, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_freq, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Number of Exercises",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_exercises$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_exercises, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_exercises, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_exercises, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Sets per Exercise",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_sets_exercise$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_sets_exercise, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_sets_exercise, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_sets_exercise, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Number of Repetitions",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_reps$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_reps, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_reps, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_reps, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = "Load (%1RM)",
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$b, 2)[2],
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$ci.lb, 2)[2],
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_load$ci.ub, 2)[2],
                                    `I^2 study` = round(I2_logCVR_hypertrophy_load, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_load, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_load, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Task Failure (No)", "Task Failure (Y)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_task_failure_y_n$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy_task_failure_y_n, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_task_failure_y_n, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_task_failure_y_n, 2)[4]),
                         data.frame(Outcome = "hypertrophy",
                                    Model = "logCVR",
                                    Moderator = c("Outcome Measure (BIA)","Outcome Measure (Biopsy: Type i)","Outcome Measure (Biopsy: Type ii)","Outcome Measure (Biopsy: Type iia)","Outcome Measure (Biopsy: Type iib)","Outcome Measure (BodPod)", "Outcome Measure (Circumference)","Outcome Measure (CT)","Outcome Measure (DXA)","Outcome Measure (Hydrostatic Weighing)","Outcome Measure (MRI)","Outcome Measure (Skinfold)","Outcome Measure (US)"),
                                    Estimate = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$b, 2),
                                    Lower = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$ci.lb, 2),
                                    Upper = round(RobuEstMultiLevelModel_logCVR_hypertrophy_measure$ci.ub, 2),
                                    `I^2 study` = round(I2_logCVR_hypertrophy_measure, 2)[2],
                                    `I^2 arm` = round(I2_logCVR_hypertrophy_measure, 2)[3],
                                    `I^2 effect` = round(I2_logCVR_hypertrophy_measure, 2)[4])
                           ) 

rownames(mods_SMD_logCVR) <- NULL

save(mods_SMD_logCVR, file = "models/mods_SMD_logCVR")

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
  select(study, arm, es, yi, vi) %>%
  left_join(Data_logCVR, by = c("study", "arm", "es")) %>%
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
