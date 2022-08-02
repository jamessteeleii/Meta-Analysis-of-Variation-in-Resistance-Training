### Refit everything in brms aswell ###


### To save trying to get hold of all the data for these studies, conduct a sensitivity analysis
# those where it was possible to directly calculate the change score standard deviation vs those where it is estimated

# Need to sort out SMD calculations for those with pre-post, and those with pre-change

# also need to now add in paired logVR and logCVR calculations

            ### All models to run
            
            # Main treatment effects - both as d/g and as RR (on % improvement scale)
            # logVR for delta sd
            # But mean-variance relationship can be an issue
            # Show example using the reps to failure dataset
            # Show relationships between m/sd and log(m)/log(sd) for deltas
            # logCVR for delta sd
            # Model log(sd) with log(mean)+group as moderators and random slopes by log(mean) with random intercepts for study/group/es and condition
            # add size to points on meta-scatter
            # logVR and logCVR for post_m (assume randomisation deals with baseline imbalance, but maybe check with pre-scores as moderator)
            # Introduce meta-regression of intervention characteristics on variance

library(metafor)
library(scales)
library(tidyverse)
library(patchwork)
library(rms)
library(lspline)

# Read csv as data frame into environment - Note: change source address
Data <- read.csv("C:/Users/James/Dropbox/Research/Response Heterogeneity - RT/RT Extracted Data.csv")

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

# Impute median within group pre-post correlations assumptions to studies
# We'll filter to > -1 < 1 as there are some odd values likely due to misreporting or miscalcuations in original studies
ri <- c(Data$RT_ri, Data$CON_ri)
ri <- subset(ri, ri >-1 & ri <1)
Data$ri_avg <- as.numeric(strrep(median(ri, na.rm = TRUE), 1))

# Estimate change score difference SD where only pre-post data available
Data$RT_delta_sd <- replmiss(Data$RT_delta_sd, with(Data, sqrt(RT_pre_sd^2 + RT_post_sd^2 - (2*ri_avg*RT_pre_sd*RT_post_sd))))
Data$CON_delta_sd <- replmiss(Data$CON_delta_sd, with(Data, sqrt(CON_pre_sd^2 + CON_post_sd^2 - (2*ri_avg*CON_pre_sd*CON_post_sd))))

### Group by design for comparative treatment standardised effect size calculations 

# Between studies 
Data_between_std <- subset(Data, study_design == "between")
Data_between_std$SD_pool <- sqrt(((Data_between_std$RT_n - 1)*Data_between_std$RT_pre_sd^2 + (Data_between_std$CON_n - 1)*Data_between_std$CON_pre_sd^2) / (Data_between_std$RT_n + Data_between_std$CON_n - 2))

Data_between_std_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                   m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=ri_avg, data = Data_between_std)
Data_between_std_CON <- escalc(measure="SMCR", m1i=CON_post_m, 
                          m2i=CON_pre_m, sd1i=SD_pool, ni=CON_n, ri=ri_avg, data = Data_between_std)

Data_between_std$yi <- (Data_between_std_RT$yi - Data_between_std_CON$yi)
Data_between_std$vi <- (Data_between_std_RT$vi + Data_between_std_CON$vi)


# Within participant studies
Data_within_std <- subset(Data, study_design == "within")
Data_within_std$SD_pool <- (((Data_within_std$RT_n - 1)*Data_within_std$RT_pre_sd) + ((Data_within_std$CON_n - 1)*Data_within_std$CON_pre_sd)) / (Data_within_std$RT_n + Data_within_std$CON_n - 2)

Data_within_std_RT <- escalc(measure="SMCR", m1i=RT_post_m, 
                          m2i=RT_pre_m, sd1i=SD_pool, ni=RT_n, ri=ri_avg, data = Data_within_std)
Data_within_std_CON <- escalc(measure="SMCR", m1i=CON_pre_m, 
                            m2i=CON_post_m, sd1i=SD_pool, ni=CON_n, ri=ri_avg, data = Data_within_std)

Data_within_std$yi <- (Data_within_std_RT$yi - Data_within_std_CON$yi)
Data_within_std$vi <- (Data_within_std_RT$vi + Data_within_std_CON$vi)

# recombine standardised effect size 
Data <- rbind(Data_between_std,Data_within_std)

### Strength
Data_strength <- Data %>% 
    filter(!is.na(yi) &  outcome == "strength")

MultiLevelModel_strength <- rma.mv(yi, V=vi, data=Data_strength,
                                         slab=paste(label),
                                         random = list(~ 1 | study/group/es), method="REML",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_strength$vi)
X <- model.matrix(MultiLevelModel_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength <- 100 * sum(MultiLevelModel_strength$sigma2) / (sum(MultiLevelModel_strength$sigma2) + (MultiLevelModel_strength$k-MultiLevelModel_strength$p)/sum(diag(P)))
I2bw_strength <- 100 * MultiLevelModel_strength$sigma2 / (sum(MultiLevelModel_strength$sigma2) + (MultiLevelModel_strength$k-MultiLevelModel_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength <- robust(MultiLevelModel_strength, Data_strength$study)

### create plot
par(mar=c(4,4,1,2))

### 2x1 panel
par(mfrow=c(2,1))

forest(RobuEstMultiLevelModel_strength$yi, RobuEstMultiLevelModel_strength$vi,
       alim = c(-4,8),
       at = c(-1,-0.5,0,0.5,1,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0),
       xlim=c(-4.5,8),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_strength$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Hedges g (Postive values indicate greater treatment response in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_strength, -1, mlab="", annotate=FALSE, cex=1, addpred = TRUE)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-4.5, -1, pos=4, cex=0.75, bquote(paste("Strength Outcomes (Q = ",
                                           
                                           .(formatC(RobuEstMultiLevelModel_strength$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_strength$k - RobuEstMultiLevelModel_strength$p),
                                           
                                           ", p = ", .(formatC(RobuEstMultiLevelModel_strength$QEp, digits=2, format="f")),
                                           
                                           ", ", I^2," = ", .(formatC(I2_strength, digits=2, format="f")), "%)")))


### Hypertrophy
Data_hypertrophy <- Data %>% 
    filter(!is.na(yi) &  outcome == "hypertrophy")

MultiLevelModel_hypertrophy <- rma.mv(yi, V=vi, data=Data_hypertrophy,
                                   slab=paste(label),
                                   random = list(~ 1 | study/group/es), method="REML",
                                   control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_hypertrophy$vi)
X <- model.matrix(MultiLevelModel_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy <- 100 * sum(MultiLevelModel_hypertrophy$sigma2) / (sum(MultiLevelModel_hypertrophy$sigma2) + (MultiLevelModel_hypertrophy$k-MultiLevelModel_hypertrophy$p)/sum(diag(P)))
I2bw_hypertrophy <- 100 * MultiLevelModel_hypertrophy$sigma2 / (sum(MultiLevelModel_hypertrophy$sigma2) + (MultiLevelModel_hypertrophy$k-MultiLevelModel_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy <- robust(MultiLevelModel_hypertrophy, Data_hypertrophy$study)

### create plot
par(mar=c(4,4,1,2))

forest(RobuEstMultiLevelModel_hypertrophy$yi, RobuEstMultiLevelModel_hypertrophy$vi,
       alim = c(-4,8),
       at = c(-1,-0.5,0,0.5,1,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0),
       xlim=c(-4.5,8),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_hypertrophy$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Hedges g (Postive values indicate greater treatment response in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_hypertrophy, -1, mlab="", annotate=FALSE, cex=1, addpred = TRUE)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-4.5, -1, pos=4, cex=0.75, bquote(paste("Hypertrophy Outcomes (Q = ",
                                           
                                           .(formatC(RobuEstMultiLevelModel_hypertrophy$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_hypertrophy$k - RobuEstMultiLevelModel_hypertrophy$p),
                                           
                                           ", p = ", .(formatC(RobuEstMultiLevelModel_hypertrophy$QEp, digits=2, format="f")),
                                           
                                           ", ", I^2," = ", .(formatC(I2_hypertrophy, digits=2, format="f")), "%)")))


################### SDs - compare delta SDs for RT vs CON to determine if 'true' interindividual variation in response exists


# Log Variability Ratios
Data_logVR <- escalc(measure = "VR", sd1i = RT_delta_sd, sd2i = CON_delta_sd, n1i = RT_n, n2i = CON_n, data = Data)
Data_logVR <- Data_logVR %>% 
    filter(!is.na(yi))

### Strength
MultiLevelModel_logVR_strength <- rma.mv(yi, V=vi, data=subset(Data_logVR, outcome == "strength"),
                                          slab=paste(label),
                                          random = list(~ 1 | study/group/es), method="REML",
                                         control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logVR, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logVR_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logVR_strength <- 100 * sum(MultiLevelModel_logVR_strength$sigma2) / (sum(MultiLevelModel_logVR_strength$sigma2) + (MultiLevelModel_logVR_strength$k-MultiLevelModel_logVR_strength$p)/sum(diag(P)))
I2bw_logVR_strength <- 100 * MultiLevelModel_logVR_strength$sigma2 / (sum(MultiLevelModel_logVR_strength$sigma2) + (MultiLevelModel_logVR_strength$k-MultiLevelModel_logVR_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_strength <- robust(MultiLevelModel_logVR_strength, subset(Data_logVR, outcome == "strength")$study)

### create plot
par(mar=c(4,4,1,2))

### 2x1 panel
par(mfrow=c(2,1))

forest(RobuEstMultiLevelModel_logVR_strength$yi, RobuEstMultiLevelModel_logVR_strength$vi,
       # alim = c(-4,3),
       # at = c(-3.0,-2.5,-2.0,-1.5,-1,-0.5,0,0.5,1,1.5,2.0,2.5),
       xlim=c(-4,2),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_logVR_strength$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Log Variability Ratio (Postive values indicate greater variability in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_logVR_strength, -0, mlab="", annotate=FALSE, cex=1)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-4, -0, pos=4, cex=0.75, bquote(paste("Strength Outcomes (Q = ",
                                           
                                           .(formatC(RobuEstMultiLevelModel_logVR_strength$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_logVR_strength$k - RobuEstMultiLevelModel_logVR_strength$p),
                                           
                                           ", p = ", .(formatC(RobuEstMultiLevelModel_logVR_strength$QEp, digits=2, format="f")),
                                           
                                           ", ", I^2," = ", .(formatC(I2_logVR_strength, digits=2, format="f")), "%)")))


### Hypertrophy
MultiLevelModel_logVR_hypertrophy <- rma.mv(yi, V=vi, data=subset(Data_logVR, outcome == "hypertrophy"),
                                                 slab=paste(label),
                                                 random = list(~ 1 | study/group/es), method="REML",
                                            control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logVR, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logVR_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logVR_hypertrophy <- 100 * sum(MultiLevelModel_logVR_hypertrophy$sigma2) / (sum(MultiLevelModel_logVR_hypertrophy$sigma2) + (MultiLevelModel_logVR_hypertrophy$k-MultiLevelModel_logVR_hypertrophy$p)/sum(diag(P)))
I2bw_logVR_hypertrophy <- 100 * MultiLevelModel_logVR_hypertrophy$sigma2 / (sum(MultiLevelModel_logVR_hypertrophy$sigma2) + (MultiLevelModel_logVR_hypertrophy$k-MultiLevelModel_logVR_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logVR_hypertrophy <- robust(MultiLevelModel_logVR_hypertrophy, subset(Data_logVR, outcome == "hypertrophy")$study)

### create plot
par(mar=c(4,4,1,2))

forest(RobuEstMultiLevelModel_logVR_hypertrophy$yi, RobuEstMultiLevelModel_logVR_hypertrophy$vi,
       # alim = c(-4,3),
       # at = c(-3.0,-2.5,-2.0,-1.5,-1,-0.5,0,0.5,1,1.5,2.0,2.5),
       xlim=c(-4,2),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_logVR_hypertrophy$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Log Variability Ratio (Postive values indicate greater variability in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_logVR_hypertrophy, -0, mlab="", annotate=FALSE, cex=1)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-4, -0, pos=4, cex=0.75, bquote(paste("Hypertrophy Outcomes (Q = ",
                                           
                                           .(formatC(RobuEstMultiLevelModel_logVR_hypertrophy$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_logVR_hypertrophy$k - RobuEstMultiLevelModel_logVR_hypertrophy$p),
                                           
                                           ", p = ", .(formatC(RobuEstMultiLevelModel_logVR_hypertrophy$QEp, digits=2, format="f")),
                                           
                                           ", ", I^2," = ", .(formatC(I2_logVR_hypertrophy, digits=2, format="f")), "%)")))

# Log Coefficient of Variation Ratios

# Add missing deltas
Data$RT_delta_m <- replmiss(Data$RT_delta_m, with(Data, RT_post_m - RT_pre_m))
Data$CON_delta_m <- replmiss(Data$CON_delta_m, with(Data, CON_post_m - CON_pre_m))


# First, we check to see what the relationship between mean and variance for change scores looks like
deltas_m <- Data %>%
    select(study, group, es, RT_n, CON_n, RT_delta_m, CON_delta_m, outcome) %>%
    pivot_longer(c(RT_delta_m, CON_delta_m),
                 names_to = "key",
                 values_to = "mean")

deltas_m$key <- recode(deltas_m$key, RT_delta_m = "RT", CON_delta_m = "CON")

deltas_sd <- Data %>%
    select(study, group, es, RT_n, CON_n, RT_delta_sd, CON_delta_sd, outcome) %>%
    pivot_longer(c(RT_delta_sd, CON_delta_sd),
                 names_to = "key",
                 values_to = "sd")

deltas_sd$key <- recode(deltas_sd$key, RT_delta_sd = "RT", CON_delta_sd = "CON")

deltas <- cbind(deltas_m, sd = deltas_sd$sd) 

deltas_RT <- deltas %>%
    filter(key == "RT") %>%
    mutate(n = RT_n) %>%
    select(study, group, es, outcome, key, mean, sd, n)

deltas_CON <- deltas %>%
    filter(key == "CON") %>%
    distinct() %>%
    mutate(n = CON_n) %>%
    select(study, group, es, outcome, key, mean, sd, n)

deltas <- rbind(deltas_RT, deltas_CON)

# Calculate log SD and variance of log SD
deltas$SD_log <- log(deltas$sd) + (1/(2*(deltas$n-1)))
deltas$SD_log_vi <- (1/(2*(deltas$n-1)))

# Let's plot and look at the raw unweighted mean-variance relationships first
# We use add a simple linear regression with a knot at zero
# (don't worry, we meta-analyse properly later accounting for dependence of effects)

m_sd_strength <- ggplot(subset(deltas, outcome == "strength"), aes(x=mean, y=sd)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean, y=sd, color = key), alpha = 0.2) +
    geom_smooth(aes(color = key, fill = key), formula = y ~ lspline(x, 1), method = "lm") +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_y_continuous(limits = c(0,1000)) +
    labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Condition") +
    guides(fill = "none") +
    theme_classic() +
    ggtitle("Strength Outcomes")

m_sd_hypertrophy <- ggplot(subset(deltas, outcome == "hypertrophy"), aes(x=mean, y=sd)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=mean, y=sd, color = key), alpha = 0.2) +
    geom_smooth(aes(color = key, fill = key), formula = y ~ lspline(x, 0), method = "lm") +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Mean Change Score", y = "Standard Deviation of the Change Score", color = "Condition") +
    guides(fill = "none") +
    theme_classic() +
    ggtitle("Hypertrophy Outcomes")

(m_sd_strength / m_sd_hypertrophy)

# Make positive as the relationship between mean and variance seems to roughly hold in both directions
deltas$mean <- ifelse(deltas$mean < 0, deltas$mean * -1, deltas$mean)

# Now let's examine raw unweighted log mean - log variance relationships

m_sd_strength_log <- ggplot(subset(deltas, outcome == "strength"), aes(x=log(mean), y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=log(mean), y=SD_log, color = key), alpha = 0.2) +
    geom_smooth(aes(color = key, fill = key), method = "lm") +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition") +
    guides(fill = "none") +
    theme_classic() 

m_sd_hypertrophy_log <- ggplot(subset(deltas, outcome == "hypertrophy"), aes(x=log(mean), y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
    geom_vline(aes(xintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 1), alpha = 0.1, lty = "dashed") +
    geom_point(aes(x=log(mean), y=SD_log, color = key), alpha = 0.2) +
    geom_smooth(aes(color = key, fill = key), method = "lm") +
    scale_fill_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    scale_color_manual("Condition", values = alpha(c("Black", "#E69F00"),0.5)) +
    labs(x = "Log Mean Change Score", y = "Log Standard Deviation of the Change Score", color = "Condition") +
    guides(fill = "none") +
    theme_classic() 


# Plot together
((m_sd_strength / m_sd_strength_log) | (m_sd_hypertrophy / m_sd_hypertrophy_log)) + plot_layout(guides = 'collect')

### There seems to be a relationship, so we'll produce some models using the log CVR 
# Note, we'll truncate the forest plots x-axes because some of the CIs are relatively wide.

### Log CVR models
Data_logCVR <- escalc(measure = "CVR", m1i = RT_delta_m, m2i = CON_delta_m,
                      sd1i = RT_delta_sd, sd2i = CON_delta_sd, n1i = RT_n, n2i = CON_n, data = Data)
Data_logCVR <- Data_logCVR %>% 
    filter(!is.na(yi))

### Strength
MultiLevelModel_logCVR_strength <- rma.mv(yi, V=vi, data=subset(Data_logCVR, outcome == "strength"),
                                          slab=paste(label),
                                          random = list(~ 1 | study/group/es), method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR, outcome == "strength")$vi)
X <- model.matrix(MultiLevelModel_logCVR_strength)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_strength <- 100 * sum(MultiLevelModel_logCVR_strength$sigma2) / (sum(MultiLevelModel_logCVR_strength$sigma2) + (MultiLevelModel_logCVR_strength$k-MultiLevelModel_logCVR_strength$p)/sum(diag(P)))
I2bw_logCVR_strength <- 100 * MultiLevelModel_logCVR_strength$sigma2 / (sum(MultiLevelModel_logCVR_strength$sigma2) + (MultiLevelModel_logCVR_strength$k-MultiLevelModel_logCVR_strength$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_strength <- robust(MultiLevelModel_logCVR_strength, subset(Data_logCVR, outcome == "strength")$study)

### create plot
par(mar=c(4,4,1,2))

### 2x1 panel
par(mfrow=c(2,1))

forest(RobuEstMultiLevelModel_logCVR_strength$yi, RobuEstMultiLevelModel_logCVR_strength$vi,
       olim = c(-20,20),
       # at = c(-3.0,-2.5,-2.0,-1.5,-1,-0.5,0,0.5,1,1.5,2.0,2.5),
       xlim=c(-4,2),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_logCVR_strength$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Log Coefficient of Variability Ratio (Postive values indicate greater variability in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_logCVR_strength, -1, mlab="", annotate=FALSE, cex=1)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-20, -1, pos=4, cex=0.75, bquote(paste("Strength Outcomes (Q = ",
                                            
                                            .(formatC(RobuEstMultiLevelModel_logCVR_strength$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_logCVR_strength$k - RobuEstMultiLevelModel_logCVR_strength$p),
                                            
                                            ", p = ", .(formatC(RobuEstMultiLevelModel_logCVR_strength$QEp, digits=2, format="f")),
                                            
                                            ", ", I^2," = ", .(formatC(I2_logCVR_strength, digits=2, format="f")), "%)")))


### Hypertrophy
MultiLevelModel_logCVR_hypertrophy <- rma.mv(yi, V=vi, data=subset(Data_logCVR, outcome == "hypertrophy"),
                                             slab=paste(label),
                                             random = list(~ 1 | study/group/es), method="REML",
                                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/subset(Data_logCVR, outcome == "hypertrophy")$vi)
X <- model.matrix(MultiLevelModel_logCVR_hypertrophy)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_logCVR_hypertrophy <- 100 * sum(MultiLevelModel_logCVR_hypertrophy$sigma2) / (sum(MultiLevelModel_logCVR_hypertrophy$sigma2) + (MultiLevelModel_logCVR_hypertrophy$k-MultiLevelModel_logCVR_hypertrophy$p)/sum(diag(P)))
I2bw_logCVR_hypertrophy <- 100 * MultiLevelModel_logCVR_hypertrophy$sigma2 / (sum(MultiLevelModel_logCVR_hypertrophy$sigma2) + (MultiLevelModel_logCVR_hypertrophy$k-MultiLevelModel_logCVR_hypertrophy$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_logCVR_hypertrophy <- robust(MultiLevelModel_logCVR_hypertrophy, subset(Data_logCVR, outcome == "hypertrophy")$study)

### create plot
par(mar=c(4,4,1,2))

forest(RobuEstMultiLevelModel_logCVR_hypertrophy$yi, RobuEstMultiLevelModel_logCVR_hypertrophy$vi,
       olim = c(-20,20),
       # at = c(-3.0,-2.5,-2.0,-1.5,-1,-0.5,0,0.5,1,1.5,2.0,2.5),
       xlim=c(-4,2),       ### adjust horizontal plot region limits
       # steps = 11, 
       order = RobuEstMultiLevelModel_logCVR_hypertrophy$yi,        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       
       cex.lab=0.75, cex.axis=0.75,   ### increase size of x-axis title/labels
       lty=c("solid","blank"),
       xlab = "Log Coefficient of Variability Ratio (Postive values indicate greater variability in RT compared to CON)")

### add summary polygon at bottom and text
addpoly(RobuEstMultiLevelModel_logCVR_hypertrophy, -1, mlab="", annotate=FALSE, cex=1)

## add text with Q-value, dfs, p-value, and I^2 statistic
text(-20, -1, pos=4, cex=0.75, bquote(paste("Hypertrophy Outcomes (Q = ",
                                            
                                            .(formatC(RobuEstMultiLevelModel_logCVR_hypertrophy$QE, digits=2, format="f")), ", df = ", .(RobuEstMultiLevelModel_logCVR_hypertrophy$k - RobuEstMultiLevelModel_logCVR_hypertrophy$p),
                                            
                                            ", p = ", .(formatC(RobuEstMultiLevelModel_logCVR_hypertrophy$QEp, digits=2, format="f")),
                                            
                                            ", ", I^2," = ", .(formatC(I2_logCVR_hypertrophy, digits=2, format="f")), "%)")))

### The log CVR models have far less heterogeneity it seems 

### Let's fit some additional models though addressing some of the possible limitations of log CVR
# See "Limitations of lnCVR and an alternative approach" in https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12309 

### Strength
deltas_strength <- deltas %>%
    filter(!is.na(sd) &
               !is.na(mean)) %>%
    filter(!is.infinite(SD_log) &
               !is.infinite(SD_log_vi)) %>%
    mutate(mean_log = log(mean)) %>%
    filter(!is.na(mean_log) &
               !is.infinite(mean_log)) %>%
    filter(outcome == "strength")

MultiLevelModel_deltas_strength <- rma.mv(SD_log, V=SD_log_vi, data=deltas_strength,
                                          random = list(~ 1 | study/group/es, ~1 | key),
                                          mods = ~ mean_log + key,
                                          method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))


### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_deltas_strength <- robust(MultiLevelModel_deltas_strength, deltas_strength$study)


### Meta-analytic scatter plot

# get the predicted log values
deltas_strength_log <- cbind(deltas_strength, predict(RobuEstMultiLevelModel_deltas_strength)) %>%
    mutate(wi = 1/sqrt(SD_log_vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

model_m_sd_strength_log <- ggplot(deltas_strength_log, aes(x=mean_log, y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
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
deltas_hypertrophy <- deltas %>%
    filter(!is.na(sd) &
               !is.na(mean)) %>%
    filter(!is.infinite(SD_log) &
               !is.infinite(SD_log_vi)) %>%
    mutate(mean_log = log(mean)) %>%
    filter(!is.na(mean_log) &
               !is.infinite(mean_log)) %>%
    filter(outcome == "hypertrophy")

MultiLevelModel_deltas_hypertrophy <- rma.mv(SD_log, V=SD_log_vi, data=deltas_hypertrophy,
                                          random = list(~ 1 | study/group/es, ~1 | key),
                                          mods = ~ mean_log + key,
                                          method="REML",
                                          control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_deltas_hypertrophy <- robust(MultiLevelModel_deltas_hypertrophy, deltas_hypertrophy$study)


### Meta-analytic scatter plot

# get the predicted log values
deltas_hypertrophy_log <- cbind(deltas_hypertrophy, predict(RobuEstMultiLevelModel_deltas_hypertrophy)) %>%
    mutate(wi = 1/sqrt(SD_log_vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

model_m_sd_hypertrophy_log <- ggplot(deltas_hypertrophy_log, aes(x=mean_log, y=SD_log)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
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

((model_m_sd_strength_log) |
        (model_m_sd_hypertrophy_log)) + plot_layout(guides = 'collect')

# Most studies in our field care about detecting treatment effects, not variance
# So we'll explore small study bias first for these across all effects
# We'll create contour enhanced funnel plots for examination of publication bias 

# run that model
Data <- Data %>%
    filter(vi > 0) # filter NAs

MultiLevelModel_all <- rma.mv(yi, V=vi, data=Data,
                                            slab=paste(label),
                                            random = list(~ 1 | study/group/es), method="REML")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_all <- robust(MultiLevelModel_all, Data$study)

### 1x3 panel
par(mfrow=c(1,3))

### produce funnel plot
funnelplotdata <- funnel(RobuEstMultiLevelModel_all, yaxis = "sei", 
                         main="Funnel plot of all treatment effects",
                         xlab = "Between Condition Treatment Effect Comparison (Hedge's g; Postive values favour RT)",
                         level=c(90, 95, 99), 
                         shade=c("white", "lightgray", "darkgray"), col = alpha(0.01), back = "white",
                         refline=0, legend=TRUE) # Contour-enhanced funnel plot

# Plot over points
with(funnelplotdata, points(x, y, col = alpha("black",0.5), pch = 19))


# Let's also have a little look at the log VR and log CVR and add them to the plot

MultiLevelModel_all_logVR <- rma.mv(yi, V=vi, data=Data_logVR,
                                     slab=paste(label),
                                     random = list(~ 1 | study/group/es), method="REML")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_all_logVR <- robust(MultiLevelModel_all_logVR, Data_logVR$study)

funnel(RobuEstMultiLevelModel_all_logVR,
       main="Funnel plot of all log VR effects",
       xlab = "Log Variability Ratio (Postive values indicate greater variability in RT compared to CON)",
       back = "white", col = alpha("black",0.5))

MultiLevelModel_all_logCVR <- rma.mv(yi, V=vi, data=Data_logCVR,
                                    slab=paste(label),
                                    random = list(~ 1 | study/group/es), method="REML")

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_all_logCVR <- robust(MultiLevelModel_all_logCVR, Data_logCVR$study)

funnel(RobuEstMultiLevelModel_all_logCVR,
       main="Funnel plot of all log CVR effects",
       xlab = "Log Coefficient of Variability Ratio (Postive values indicate greater variability in RT compared to CON)",
       back = "white", col = alpha("black",0.5))

################ Outputs

### summarise all descriptives
sample_sizes <- Data %>%
        filter(study_design == "between") %>%
        select(group, RT_n, CON_n) %>%
        group_by(group) %>%
        summarise(RT_n = max(RT_n), 
                  CON_n = max(CON_n)) %>%
    summarise(n_all_RT = sum(RT_n),
              n_min_RT = min(RT_n),
              n_median_RT = median(RT_n),
              n_max_RT = max(RT_n),
              n_all_CON = sum(CON_n, na.rm =TRUE), 
              n_min_CON = min(CON_n, na.rm =TRUE),
              n_median_CON = median(CON_n, na.rm =TRUE),
              n_max_CON = max(CON_n, na.rm =TRUE))

### Descriptive Tables
library(gtsummary)
library(flextable)

Data_characteristics <- Data %>% 
    distinct(group, .keep_all = TRUE) %>%
    select(age,
           sex_._male,
           bmi,
           train_status,
           weeks,
           freq,
           exercises,
           sets_exercise,
           reps,
           load,
           task_failure_y_n,
           measure
           )

Data_characteristics %>% tbl_summary() %>% as_flex_table() %>% save_as_docx(path = "Table_characteristics.docx")

# Descriptive characteristics
sample_sizes

### Primary Outcomes ###

# Log variability ratio
RobuEstMultiLevelModel_strength
I2_strength

RobuEstMultiLevelModel_hypertrophy
I2_hypertrophy

RobuEstMultiLevelModel_logVR_strength
I2_logVR_strength

RobuEstMultiLevelModel_logVR_hypertrophy
I2_logVR_hypertrophy

RobuEstMultiLevelModel_logCVR_strength
I2_logCVR_strength

RobuEstMultiLevelModel_logCVR_hypertrophy
I2_logCVR_hypertrophy



############## Exploring moderators of treatment effects and variances

### Duration
Data_strength_duration <- Data_strength %>%
    filter(!is.na(weeks))

MultiLevelModel_strength_duration <- rma.mv(yi, V=vi, data=Data_strength_duration,
                                                 random = list(~ 1 | study/group/es), mods = ~ log(weeks),
                                                 method="REML", control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_strength_duration$vi)
X <- model.matrix(MultiLevelModel_strength_duration)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_strength_duration <- 100 * sum(MultiLevelModel_strength_duration$sigma2) / (sum(MultiLevelModel_strength_duration$sigma2) + (MultiLevelModel_strength_duration$k-MultiLevelModel_strength_duration$p)/sum(diag(P)))
I2bw_strength_duration <- 100 * MultiLevelModel_strength_duration$sigma2 / (sum(MultiLevelModel_strength_duration$sigma2) + (MultiLevelModel_strength_duration$k-MultiLevelModel_strength_duration$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_strength_duration <- robust(MultiLevelModel_strength_duration, Data_strength_duration$study)

### Meta-analytic scatter plot

# get the predicted values
Data_strength_duration <- cbind(Data_strength_duration, predict(RobuEstMultiLevelModel_strength_duration)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

meta_reg_strength_duration <- Data_strength_duration %>%
    ggplot(aes(x=log(weeks), y=pred)) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_ribbon(aes(x=log(weeks), ymax=ci.ub, ymin=ci.lb), alpha = 0.2) +
    geom_line(aes(x=log(weeks), y=pred)) +
    geom_point(aes(x=log(weeks), y=yi, size = size), alpha = 0.2) +
    labs(x = "Intervention Length (Weeks)", y = "Hedge's g (Postive values indicate greater treatment response in RT compared to CON)", size = "") +
    theme_classic() +
    guides(size = "none") +
    theme(
    # legend.text = element_text(size = 6),
    #       legend.title = element_text(size = 8),
          axis.title = element_text(size = 8),
    #       axis.text = element_text(size = 6),                    plot.title = element_text(size = 10) 
    ) +
    ggtitle("Strength Outcomes") 

Data_hypertrophy_duration <- Data_hypertrophy %>%
    filter(!is.na(weeks))

MultiLevelModel_hypertrophy_duration <- rma.mv(yi, V=vi, data=Data_hypertrophy_duration,
                                            random = list(~ 1 | study/group/es), mods = ~ log(weeks),
                                            method="REML", control=list(optimizer="optim", optmethod="Nelder-Mead"))

### Calculate I^2 for see "http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate"
W <- diag(1/Data_hypertrophy_duration$vi)
X <- model.matrix(MultiLevelModel_hypertrophy_duration)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_hypertrophy_duration <- 100 * sum(MultiLevelModel_hypertrophy_duration$sigma2) / (sum(MultiLevelModel_hypertrophy_duration$sigma2) + (MultiLevelModel_hypertrophy_duration$k-MultiLevelModel_hypertrophy_duration$p)/sum(diag(P)))
I2bw_hypertrophy_duration <- 100 * MultiLevelModel_hypertrophy_duration$sigma2 / (sum(MultiLevelModel_hypertrophy_duration$sigma2) + (MultiLevelModel_hypertrophy_duration$k-MultiLevelModel_hypertrophy_duration$p)/sum(diag(P)))

### Calculate robust estimate from multi-level model
RobuEstMultiLevelModel_hypertrophy_duration <- robust(MultiLevelModel_hypertrophy_duration, Data_hypertrophy_duration$study)

### Meta-analytic scatter plot

# get the predicted values
Data_hypertrophy_duration <- cbind(Data_hypertrophy_duration, predict(RobuEstMultiLevelModel_hypertrophy_duration)) %>%
    mutate(wi = 1/sqrt(vi),
           size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

meta_reg_hypertrophy_duration <- Data_hypertrophy_duration %>%
    ggplot(aes(x=log(weeks), y=pred)) +
    geom_vline(aes(xintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_hline(aes(yintercept = 0), alpha = 0.1, lty = "dashed") +
    geom_ribbon(aes(x=log(weeks), ymax=ci.ub, ymin=ci.lb), alpha = 0.2) +
    geom_line(aes(x=log(weeks), y=pred)) +
    geom_point(aes(x=log(weeks), y=yi, size = size), alpha = 0.2) +
    labs(x = "Intervention Length (Weeks)", y = "Hedge's g (Postive values indicate greater treatment response in RT compared to CON)", size = "") +
    theme_classic() +
    guides(size = "none") +
    theme(
        # legend.text = element_text(size = 6),
        #       legend.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        #       axis.text = element_text(size = 6),                    plot.title = element_text(size = 10) 
    ) +
    ggtitle("Hypertrophy Outcomes") 


meta_reg_strength_duration / meta_reg_hypertrophy_duration
