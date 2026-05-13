# 1.1 Load data ----------------------------------------------------------------
my_data <- read.csv("https://github.com/katjia/ISPPD_TTE_workshop/raw/main/dat.csv")

# Display data with descriptive column names
my_data_display <- my_data
names(my_data_display)[names(my_data_display) == "time.inf.since.dose.1"] <- "time.inf.since.dose.1 \n (time to infection from Dose 1)"
names(my_data_display)[names(my_data_display) == "days.btw.doses.12"] <- "days.btw.doses.12 \n (time to Dose 2 from Dose 1)"

# Print title
title <- "                Time to infection and time to Dose 2 among those who received Dose 1                  "
separator <- paste(rep("=", nchar(title)), collapse = "")
cat(separator, "\n")
cat(title, "\n")
cat(separator, "\n")
head(my_data_display, 15)

# 1.2 Package setup ------------------------------------------------------------
install_missing_pkgs <- function(pkgs) {
  new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, dependencies = c("Depends","Imports"))
  }
}

load_pkgs <- function(pkgs) {
  invisible(lapply(pkgs, library, character.only = TRUE))
}

pkgs <- c(
  "dplyr",
  "tidyr",
  "data.table",
  "survival",
  "purrr",
  "future",
  "future.apply",
  "patchwork",
  "ggplot2"
)

install_missing_pkgs(pkgs)
load_pkgs(pkgs)

# 2 Cloning & censoring --------------------------------------------------------
# 2.1 Cloning
# 👉 At the index date (i.e., day of receiving Dose 1), 
# we assign each person to all protocols that are compatible with their observed data.
# 👥 We create two "clones" of each person under two protocols:
# 1. 2-dose protocol: received two doses
# 2. 1-dose protocol: received only dose

# 2.2 Censoring 
# Individuals are **censored** when they deviate from the assigned protocol.

# 🥾 BOOTSTRAPPING is used to construct confidence intervals
# Bootstrap resampling: Randomly sample individuals WITH REPLACEMENT n_boot times
# to estimate uncertainty in cumulative incidence curves.
# Each iteration applies the full clone-censor-weight procedure.
# We then calculate 95% CI from the 2.5th and 97.5th percentiles across iterations.

n_boot <- 500 

# Each uncensored individual is followed up until 240 days after the index date (Dose 1).

followup_period <- 240

# Create placeholders
all_cases_2_doses <- list()
all_cases_1_dose <- list()

all_res_2_doses <- list()
all_res_1_dose <- list()

for(i in 1:n_boot){

# Set a seed for reproducibility
seed <- 788 + i

set.seed(seed) 

target_pop <- my_data[sample(1:dim(my_data)[1], replace=TRUE), ] 

# 2.1-2.2 (a) Cloning & censoring for the 2-dose schedule-----------------------

df_2_doses <- target_pop %>% 
  
  ### Step 1 for the 2-dose schedule ###
  
  # *---*---*---*---* KEY code for cloning & censoring *---*---*---*---*---*
  
  # Calculate censoring time for each clone under each protocol
  dplyr::mutate(censoring_time = case_when(
    
    #---Criteria 1---#
    # If never received dose 2: follow until Day 28
    is.na(days.btw.doses.12) ~ 28,
    
    #---Criteria 2---#                                                                
    # If received dose 2: follow until end of follow-up
    !is.na(days.btw.doses.12) ~ followup_period))

  # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*

    ### Create a function for Steps 2-4 ###
    compute_metrics <- function(df) {
    
      # Step 2: Determine if the outcome occurred BEFORE censoring for each clone
      df <- df %>%
        dplyr::mutate(outcome_before_censoring = case_when(
          # Outcome_before_censoring is 1 if the outcome occurred before 
          # an individual was censored
          time.inf.since.dose.1 < censoring_time ~ 1,
          
          # Otherwise, the individual did not have an outcome, or the individual
          # had an outcome after censoring_time 
          # For either case, outcome_before_censoring = 0.
          is.na(df$time.inf.since.dose.1) | df$time.inf.since.dose.1 >= df$censoring_time ~ 0
        ))
      
      # Step 3: Determine the duration of follow up for each clone
      # The follow-up duration "dur_followup" is the time from dose 1 until either: 
      # - censoring_time (censoring due to non-adherence or end-of-follow-up) OR
      # - time.inf.since.dose.1 (outcome occurred before censoring)
      df$dur_followup <- df$censoring_time
      
      df$dur_followup[!is.na(df$outcome_before_censoring) & df$outcome_before_censoring == 1L] <-
        df$time.inf.since.dose.1[!is.na(df$outcome_before_censoring) & df$outcome_before_censoring == 1L]
      
      ### Step 4: Create censoring indicator variable
      # This indicator is used in the Cox model for censoring weights
      # censored = 1: The clone was censored due to the protocol non-adherence or end of follow up
      # censored = 0: otherwise
      df$censored <- as.integer(df$outcome_before_censoring == 0L)
      
      df
    }
    
    df_2_doses <- compute_metrics(df_2_doses)

    # Label protocol group
    df_2_doses$protocol <- "2 doses"

# 2.1-2.2 (b) Cloning for the 1-dose schedule-----------------------------------
df_1_dose <- target_pop %>% 
      
    ### Step 1 for the 1-dose schedule ###
      
    # *---*---*---*---* KEY code for cloning & censoring *---*---*---*---*---*
      
    # Calculate censoring time for each clone
      dplyr::mutate(censoring_time = case_when(
        
        #---Criteria 1---#
        # If never received dose 2: follow until end of follow-up 
        is.na(days.btw.doses.12) ~ followup_period,
        
        #---Criteria 2---#                                                                
        # If received dose 2: until Day 28 (day of receiving dose 2)
        !is.na(days.btw.doses.12) ~ 28))
    
    # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
    
    ### Steps 2-4 ###
    df_1_dose <- compute_metrics(df_1_dose)
    
    # Label protocol group
    df_1_dose$protocol <- "1 dose"

# 3 Calculate censor weights ---------------------------------------------------
  # After cloning, we censor clones who deviate from the assigned protocol.
  # 🎯 However, our goal is to construct a hypothetical population in which 
  # nobody is censored because everybody followed their assigned protocol.
  # To do this, we upweight uncensored individuals by giving them a weight equal 
  # to the inverse of their probability of being uncensored. That is, people who 
  # are censored transfer their weights to those who are uncensored.

  # We do so in three steps:
    
  # 3 (a) Calculate censor weights for 2-dose protocol -------------------------
    
    # *---*---*---*---*---* KEY code for calculating weights *---*---*---*---*---*
    
    # STEP 1: 
    # Fit a Cox proportional hazard (PH) model, for which the "event" is CENSORING, 
    # to estimate the probability of remaining uncensored over time.
    
      coxph.censor <- survival::coxph(survival::Surv(dur_followup, censored) ~ 1, data = df_2_doses)
    
    # where
    #   - dur_followup: time until censoring (non-adherence or end-of-follow-up) or outcome
    #   - censored: 1 if censored, 0 if had outcome
    #   - ~ 1: null model (no covariates), estimates baseline censoring hazard only  
    
    # Subset the clone population to just those who had the outcome because 
    # total population size is fixed and we do not need to estimate weights for 
    # the non-cases.
    cases <- df_2_doses[df_2_doses$outcome_before_censoring==1,]
      
    # STEP 2: Predict the probability of remaining uncensored at each person’s event time
    # (date of the outcome).
    # For a null Cox model, use baseline hazard
    bh <- survival::basehaz(coxph.censor, centered = FALSE)
    expected_vals <- approx(bh$time, bh$hazard, xout = cases$dur_followup, rule = 2)$y
    cases_2_doses <- cases %>%
      dplyr::mutate(.fitted = expected_vals,
                    prob = exp(-.fitted))
      
    # Order the subset data by the event time 
    cases_2_doses <- cases_2_doses[order(cases_2_doses$dur_followup),]
      
    # STEP 3: Compute the inverse probability censoring weights
    cases_2_doses$wt <- 1/cases_2_doses$prob
      
    # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
      
    # Create time data 
    cases_2_doses <- cases_2_doses %>% 
      subset(., select = c(ID, wt, dur_followup, protocol))
      
    cases_2_doses <- cbind(cases_2_doses, 
                           sim=i)
      
    cases_2_doses <- data.frame(dur_followup=seq(1, followup_period, 1)) %>% 
        left_join(., cases_2_doses, by = "dur_followup")
    
    # *---*---* KEY code for calculating cumulative incidence curve *---*---*
    
    # At each time t: Cumulative Incidence(t) = Cumulative sum of weighted cases / Total population size
    # This accounts for censoring by upweighting cases who were uncensored
    # This gives us our cumulative incidence curve
    cases_2_doses$risk <- cumsum(tidyr::replace_na(cases_2_doses$wt, 0)) / nrow(df_2_doses)
    
    # where:
    # - wt = weight assigned to each uncensored case
    # - cumsum(...) = cumulative sum over time
    # - nrow(df_2_doses) = total population size (20,000)
    
    # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
    
    res_2_doses <- cases_2_doses[which(cases_2_doses$dur_followup <= followup_period),] %>%
        dplyr::group_by(dur_followup) %>%
        dplyr::summarise(cumrisk = max(risk, na.rm=TRUE))
      
    # Add the protocol column
    res_2_doses$protocol <- "2 doses"
    cases_2_doses$protocol <- "2 doses"
      
    # Note iteration 
    res_2_doses$sim <- i
      
    # Combine results across all simulations
    all_res_2_doses <- rbind(all_res_2_doses, res_2_doses)
      
    # Save the patient-level weights
    all_cases_2_doses <- rbind(all_cases_2_doses, cases_2_doses)

  # 3 (b) Censoring weights for 1-dose protocol --------------------------------
      
    # *---*---*---*---*---* KEY code for calculating weights *---*---*---*---*---*
    
    # STEP 1
    # Fit a Cox proportional hazard (PH) model, for which the "event" is CENSORING, 
    # to estimate the probability of remaining uncensored over time.
    
    coxph.censor <- survival::coxph(survival::Surv(dur_followup, censored) ~ 1, data = df_1_dose)
      
    # Subset the clone population to just those who had the outcome
    cases <- df_1_dose[df_1_dose$outcome_before_censoring==1,]
      
    # STEP 2: Predict the probability of remaining uncensored at each person’s event time
    # (date of the outcome).
    # For a null Cox model, use baseline hazard
    bh <- survival::basehaz(coxph.censor, centered = FALSE)
    expected_vals <- approx(bh$time, bh$hazard, xout = cases$dur_followup, rule = 2)$y
    cases_1_dose <- cases %>%
        dplyr::mutate(.fitted = expected_vals,
                      prob = exp(-.fitted))
      
    # Order the subset data by the event time 
    cases_1_dose <- cases_1_dose[order(cases_1_dose$dur_followup),]
      
    # STEP 3: Compute the inverse probability censoring weights
    cases_1_dose$wt <- 1/cases_1_dose$prob
      
    # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
      
    # Create time data 
    cases_1_dose <- cases_1_dose %>% 
        subset(., select = c(ID, wt, dur_followup, protocol))
      
    cases_1_dose <- cbind(cases_1_dose, 
                             sim=i)
      
    cases_1_dose <- data.frame(dur_followup=seq(1, followup_period, 1)) %>% 
        left_join(., cases_1_dose, by = "dur_followup")
      
    # *---*---* KEY code for calculating cumulative incidence curve *---*---*
    
    cases_1_dose$risk <- cumsum(tidyr::replace_na(cases_1_dose$wt, 0)) / nrow(df_1_dose)
    
    # *---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*
    
    res_1_dose <- cases_1_dose[which(cases_1_dose$dur_followup <= followup_period),] %>%
        dplyr::group_by(dur_followup) %>%
        dplyr::summarise(cumrisk = max(risk, na.rm=TRUE))
      
    # Add the protocol column
    res_1_dose$protocol <- "1 dose"
    cases_1_dose$protocol <- "1 dose"
      
    # Note iteration 
    res_1_dose$sim <- i
      
    # Combine results across all simulations
    all_res_1_dose <- rbind(all_res_1_dose, res_1_dose)
      
    # Save the patient-level weights
    all_cases_1_dose <- rbind(all_cases_1_dose, cases_1_dose)
      
    # clean up
    rm(df_2_doses)
    rm(df_1_dose)
    
    rm(res_2_doses)
    rm(res_1_dose)
    
    rm(cases_2_doses)
    rm(cases_1_dose)
}

# 4. Visualize the weights -----------------------------------------------------
all_cases <- rbind(all_cases_2_doses, 
                   all_cases_1_dose)

# Remove rows with NA weights (time points with no outcomes)
all_cases_clean <- all_cases %>% 
  filter(!is.na(wt))

# Randomly select 3 simulations
set.seed(123)
available_sims <- unique(all_cases_clean$sim)
selected_sims <- sample(available_sims, 3)

# Filter to only the selected simulations
all_cases_selected <- all_cases_clean %>%
  filter(sim %in% selected_sims)

# Visualize the distribution of the inverse probability weights
# 2x3 grid: upper row for only 1 dose, lower row for 2 doses, columns for each sim
ggplot(all_cases_selected, aes(x = wt)) +
  geom_histogram(fill = "steelblue", color = "white", binwidth = 0.1) +
  labs(
    title = "Distribution of Weights across Selected Simulations",
    x = "Weight (1/Probability of remaining uncensored)",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  facet_grid(protocol ~ sim, labeller = labeller(sim = function(x) paste0("sim=", x))) +
  scale_x_continuous(breaks = c(1, 2)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        panel.spacing.x = unit(1.5, "lines"))

# 5. Plot cumulative incidence -------------------------------------------------
all_data <- rbind(all_res_2_doses,
                  all_res_1_dose)

# Calculate median and 95% confidence intervals by region
summary_data <- all_data %>%
  group_by(protocol, dur_followup) %>%
  summarize(
    median = quantile(cumrisk, 0.5, na.rm = T),
    lb = quantile(cumrisk, 0.025, na.rm = T),
    ub = quantile(cumrisk, 0.975, na.rm = T),
    .groups = "drop"
  )

ggplot(
  summary_data,
  aes(x = dur_followup, y = median, color = protocol, fill = protocol)
) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Days since dose 1",
    y = "Cumulative incidence",
    color = " ",
    fill  = " "
  ) +
  scale_color_manual(
    values = c("2 doses" = "dodgerblue", 
               "1 dose" = "hotpink"),
    labels = c("2 doses" = "2 doses", 
               "1 dose" = "1 dose")
  ) +
  scale_fill_manual(
    values = c("2 doses" = "dodgerblue", 
               "1 dose" = "hotpink"),
    labels = c("2 doses" = "2 doses", 
               "1 dose" = "1 dose")
  ) +
  scale_y_continuous(labels = scales::label_percent(
    accuracy = 0.1
  )) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) 

# 6. Plot risk ratio -----------------------------------------------------------
RR_time_series <- all_data %>%
  select(dur_followup, sim, protocol, cumrisk) %>%
  pivot_wider(
    names_from  = protocol,
    values_from = cumrisk
  ) %>%
  mutate(
    RR = `1 dose`/`2 doses`) %>%
  group_by(dur_followup) %>%
  summarize(
    RR_median = quantile(RR, 0.5, na.rm=T),
    RR_lb = quantile(RR, 0.025, na.rm=T),
    RR_ub = quantile(RR, 0.975, na.rm=T),
    .groups = "drop"
  )

ggplot(
  RR_time_series %>%
    select(
      dur_followup,
      median = RR_median,
      lb     = RR_lb,
      ub     = RR_ub
    ),
  aes(x = dur_followup, y = median)
) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(
    x = "Days since dose 1",
    y = "Risk ratio (1 dose vs 2 doses)"
  ) +
  scale_y_continuous(limits = c(NA, 3)) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
