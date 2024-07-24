library(rstatix)
library(afex)
library(car)
library(broom)
library(emmeans)
library(stringr)
library(lmerTest)
library(tidyverse)
library(MOTE)

# Import and prepare data 

data <- read_csv("FilipStachnik2021_replicationdata.csv") 

# Convert to long data set

long_data <- data %>%
  rename(plac_1rm = `PLAC 1RM`,
         caf3_1rm = `CAF-3 1RM`,
         caf6_1rm = `CAF-6	 1RM`) %>%
  select(id, plac_1rm, caf3_1rm, caf6_1rm) %>%
  pivot_longer(cols = c("plac_1rm", "caf3_1rm", "caf6_1rm"),
               names_to = "condition",
               values_to = "one_rep_max") 

long_data$condition <- as.factor(long_data$condition)

# Descriptives ------------

desc <- long_data %>% 
  group_by(condition) %>%
  summarize(count = n (),
            overall_mean = mean(one_rep_max,na.rm=TRUE),
            overall_sd = sd(one_rep_max, na.rm = TRUE))
desc

# Repeated Measures ANOVA -----

##afex::aov_4(continuous_var ~ group_var + (RM_var|id_var)

anova_results <- afex::aov_4(one_rep_max ~ (condition|id), 
                                   data = long_data,
                                   anova_table = list(es = "pes")) # partial eta squared
anova_results

summary(anova_results)

## Resolving assumptions --------

### Normality test 

long_data %>% 
  dplyr::group_by(condition) %>% 
  rstatix::shapiro_test(one_rep_max) # shapiro-wilk test on individual groups

norm <- performance::check_normality(anova_results)
plot(norm)
plot(norm, type = "qq")

### Outliers check

long_data %>%
  group_by(condition) %>%
  identify_outliers(one_rep_max)

## Plots

## violin

long_data %>% 
  ggplot(aes(condition, one_rep_max)) +  
  geom_violin(fill = "gray") +
  geom_boxplot(width = .07,
               fill = "white") +
  geom_jitter(position = position_jitter(0.21)) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 18,
               color = "red",
               size = 5) +
  theme_bw()

## Individual qq plots 

long_data %>% 
  ggplot(aes(sample = one_rep_max)) +    
  geom_qq() +                               
  stat_qq_line() +                          
  facet_wrap(~ condition,                   # Panel by group
             labeller = label_both) +    
  theme_bw()

# Normal dataset, remove outlier --------------- 

normal_long_data <- long_data %>% 
  filter(id!= 3) 

## Normal Descriptives ------------

normal_desc <- normal_long_data %>%
  group_by(condition) %>%
  summarize(count = n (),
            overall_mean = mean(one_rep_max,na.rm=TRUE),
            overall_sd = sd(one_rep_max, na.rm = TRUE))
normal_desc

## Repeated Measures ANOVA -----

##afex::aov_4(continuous_var ~ group_var + (RM_var|id_var)

normal_anova_results <- afex::aov_4(one_rep_max ~ (condition|id), 
                             data = normal_long_data,
                             anova_table = list(es = "pes")) # partial eta squared
normal_anova_results

summary(normal_anova_results)

## Post hoc ------------

normal_data_emm <- normal_anova_results %>% 
  emmeans::emmeans(~ condition, model = "multivariate")
normal_data_emm

normal_posthocresults <- pairs(normal_data_emm, adjust = "Tukey") %>% # original study uses Tukey
  broom::tidy(conf.int = T)
normal_posthocresults

## Resolving assumptions --------

### Normality test 

normal_long_data %>% 
  dplyr::group_by(condition) %>% 
  rstatix::shapiro_test(one_rep_max) # shapiro-wilk test on individual groups

norm <- performance::check_normality(normal_anova_results)
plot(norm)
plot(norm, type = "qq")

### Outliers check

normal_long_data %>%
  group_by(condition) %>%
  identify_outliers(one_rep_max)

## Plots

## violin

normal_long_data %>% 
  ggplot(aes(condition, one_rep_max)) +  
  geom_violin(fill = "gray") +
  geom_boxplot(width = .07,
               fill = "white") +
  geom_jitter(position = position_jitter(0.21)) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 18,
               color = "red",
               size = 5) +
  theme_bw()

## Individual qq plots 

normal_long_data %>% 
  ggplot(aes(sample = one_rep_max)) +    
  geom_qq() +                               
  stat_qq_line() +                          
  facet_wrap(~ condition,                   # Panel by group
             labeller = label_both) +    
  theme_bw()

# Replication Effect Size ----------

pes_rep_ci <- eta.F(
  dfm = normal_anova_results$anova_table$`num Df`,
  dfe = normal_anova_results$anova_table$`den Df`,
  Fvalue = normal_anova_results$anova_table$F,
  a = 0.05) %>%
  as.data.frame() %>%
  select(eta, etalow, etahigh) %>%
  mutate(study_id = c("Replication study")) # add identifier
pes_rep_ci


# Original values ------

orig_values <- data.frame(
  plac_mean = 40.48,
    plac_sd = 9.21,
    caf3_mean = 41.68,
    caf3_sd = 8.98,
    caf6_mean = 42.98,
    caf6_sd = 8.79,
    f_val = 14.74,
  df1 = 2,
  df2 = 19,
  p_val = 0.0099 # conservative estimate
) 


## Calculate original partial eta squared ------

  pes_orig_ci <- eta.F(
    dfm = orig_values$df1,
    dfe = orig_values$df2,
    Fvalue = orig_values$f_val,
    a = 0.05) %>%
    as.data.frame() %>%
    select(eta, etalow, etahigh) %>%
    mutate(study_id = c("Original study")) # add identifier
  pes_orig_ci

# Replication test -----

pes_rep = normal_anova_results$anova_table$pes
df_rep = 14
pes_ori = pes_orig_ci$eta
df_ori = orig_values$df2

rho_ori = 2*sqrt(pes_ori)-1
rho_rep = 2*sqrt(pes_rep)-1

rep_test = TOSTER::compare_cor(r1 = rho_ori,
                               df1 = df_ori,
                               r2 = rho_rep,
                               df2 = df_rep,
                               alternative = "greater")
rep_test


