# 1.1 Load data -----------------------------------------------------------------
my_data <- read.csv("https://github.com/katjia/ISPPD_TTE_workshop/raw/main/dat.csv")
head(my_data)

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
# We create "clones" of each individual under two protocols:
# 1. 2-dose protocol: received two doses
# 2. 1-dose protocol: received only dose
# 
# Cloning: Each individual is assigned to both protocols at the time they receive dose 1 (index date). 
# Censoring: Individuals are censored when they deviate from the assigned protocol. 

# BOOTSTRAPPING to construct confidence intervals
# Bootstrap resampling: Randomly sample participants WITH REPLACEMENT n_boot times
# to estimate uncertainty in cumulative incidence curves.
# Each iteration applies the full clone-censor-weight procedure.
# We then calculate 95% CI from the 2.5th and 97.5th percentiles across iterations.
n_boot <- 500 # 1000

# We follow up each individual until 240 days after dose 1
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

bootm <- my_data[sample(1:dim(my_data)[1], replace=TRUE), ] 

# 2.1 Cloning for the 2-dose schedule---------------------------------

df_2_doses <- bootm %>% 
  
  ### Step 1 for the 2-dose schedule ###
  # Calculate censoring time for each clone
  # CENSORING RULE for 2-dose protocol:
  # - Participants who never get dose 2 are CENSORED at day 28 (they didn't follow protocol)
  # - Participants who get dose 2 are followed until the end of the study period
  dplyr::mutate(censoring_time = case_when(
    
    #---Criteria 1---#
    # If never received dose 2: follow until Day 28
    is.na(days.btw.doses.12) ~ 28,
    
    #---Criteria 2---#                                                                
    # If received dose 2: follow until end of follow-up
    !is.na(days.btw.doses.12) ~ followup_period))

    ### Create a function for Steps 2-4 ###
    compute_metrics <- function(df) {
    
      # Step 2: Determine if the outcome occurred BEFORE censoring for each clone
      df <- df %>%
        dplyr::mutate(outcome_before_censoring = case_when(
          # Outcome_before_censoring is 1 if the outcome occurred before 
          # an individual was censored
          time.inf < censoring_time ~ 1,
          
          # Otherwise, the individual did not have an outcome, or the individual
          # had an outcome after censoring_time 
          # For either case, outcome_before_censoring = 0.
          is.na(df$time.inf) | df$time.inf >= df$censoring_time ~ 0
        ))
      
      # Step 3: Determine the duration of follow up for each clone
      # The follow-up duration "dur_followup" is the time from dose 1 until either: 
      # - censoring_time (censoring occurs due to non-adherence to the specified protocol) OR
      # - time.inf (outcome occurred before censoring)
      df$dur_followup <- df$censoring_time
      
      df$dur_followup[!is.na(df$outcome_before_censoring) & df$outcome_before_censoring == 1L] <-
        df$time.inf[!is.na(df$outcome_before_censoring) & df$outcome_before_censoring == 1L]
      
      ### Step 4: Create censoring indicator variable
      # This indicator is used in the Cox model for censoring weights
      # censored = 1: The clone was censored due to the protocol non-adherence or the end of follow up
      # censored = 0: otherwise
      df$censored <- as.integer(df$outcome_before_censoring == 0L)
      
      df
    }
    
    df_2_doses <- compute_metrics(df_2_doses)

    # Label protocol group
    df_2_doses$protocol <- "2 doses"

# 2.2 Cloning for the 1-dose schedule----------------------------
df_1_dose <- bootm %>% 
      
      ### Step 1 for the 1-dose schedule ###
      # Calculate censoring time for each clone
      # CENSORING RULE for 1-dose protocol:
      # - Those who never get dose 2 are followed until end of study (they followed protocol)
      # - Those who get dose 2 are CENSORED at day 28 (they didn't follow the 1-dose protocol)
      dplyr::mutate(censoring_time = case_when(
        
        #---Criteria 1---#
        # If never received dose 2: follow until end of follow-up 
        is.na(days.btw.doses.12) ~ followup_period,
        
        #---Criteria 2---#                                                                
        # If received dose 2: until Day 28 (day of receiving dose 2)
        !is.na(days.btw.doses.12) ~ 28))
    
    ### Steps 2-4 ###
    df_1_dose <- compute_metrics(df_1_dose)
    
    # Label protocol group
    df_1_dose$protocol <- "1 dose"

    # Define the outcome at dur_followup under each protocol:
    # clone.outcome=1 if a clone had an outcome before censoring under each protocol
    # clone.outcome=0 otherwise.
    df_2_doses$clone.outcome <- as.integer(df_2_doses$outcome_before_censoring == 1)
    df_1_dose$clone.outcome <- as.integer(df_1_dose$outcome_before_censoring == 1)
    
# 3 Compute the inverse probability of censoring weight (IPCW) for each clone under each protocol --------------

    # We use the Cox proportional hazard model to estimate the probability of remaining recensored over time
    # We fit a Cox proportional hazards model where the "event" is CENSORING (not outcome)
    # Model: Surv(dur_followup, censored) ~ 1
    #   - dur_followup: time until censoring (non-adherence or EOFU) or outcome
    #   - censored: 1 if censored, 0 if had outcome
    #   - ~ 1: null model (no covariates), estimates baseline censoring hazard only

  # 3.1 Censoring weights for 2-dose protocol ----------------------------------
    
      coxph.censor <- survival::coxph(survival::Surv(dur_followup, censored) ~ 1, data = df_2_doses)
      
      # Subset the clone population to just those who had the outcome
      # because those who did not have the outcome do not contribute to the calculation
      # of the cumulative risk of outcome.
      cases <- df_2_doses[df_2_doses$clone.outcome==1,] 
      
      # Predict the probability of remaining uncensored at each person’s event time
      # (date of the outcome).
      # For a null Cox model, use baseline hazard
      bh <- survival::basehaz(coxph.censor, centered = FALSE)
      expected_vals <- approx(bh$time, bh$hazard, xout = cases$dur_followup, rule = 2)$y
      cases_2_doses <- cases %>%
        dplyr::mutate(.fitted = expected_vals,
                      prob = exp(-.fitted))
      
      # Order the subset data by the event time 
      cases_2_doses <- cases_2_doses[order(cases_2_doses$dur_followup),]
      
      # Compute the inverse probability censoring weights
      cases_2_doses$wt <- 1/cases_2_doses$prob
      
      # Create time data 
      cases_2_doses <- cases_2_doses %>% 
        subset(., select = c(ID, wt, dur_followup, protocol))
      
      cases_2_doses <- cbind(cases_2_doses, 
                             sim=i)
      
      cases_2_doses <- data.frame(dur_followup=seq(1, followup_period, 1)) %>% 
        left_join(., cases_2_doses, by = "dur_followup")
      
      # At each time t: Cumulative Incidence(t) = Cumulative sum of weights / Total population size 
      # This accounts for censoring by upweighting cases who were at risk of being censored
      # This gives us our cumulative incidence curve
      cases_2_doses$risk <- cumsum(tidyr::replace_na(cases_2_doses$wt, 0)) / nrow(df_2_doses)
      
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

      # 3.2 Censoring weights for 1-dose protocol ------------------------------
      
      coxph.censor <- survival::coxph(survival::Surv(dur_followup, censored) ~ 1, data = df_1_dose)
      
      # Subset the clone population to just those who had the outcome
      # because those who did not have the outcome do not contribute to the calculation
      # of the cumulative risk of outcome.
      cases <- df_1_dose[df_1_dose$clone.outcome==1,] 
      
      # Predict the probability of remaining uncensored at each person’s event time
      # (date of the outcome).
      # For a null Cox model, use baseline hazard
      bh <- survival::basehaz(coxph.censor, centered = FALSE)
      expected_vals <- approx(bh$time, bh$hazard, xout = cases$dur_followup, rule = 2)$y
      cases_1_dose <- cases %>%
        dplyr::mutate(.fitted = expected_vals,
                      prob = exp(-.fitted))
      
      # Order the subset data by the event time 
      cases_1_dose <- cases_1_dose[order(cases_1_dose$dur_followup),]
      
      # Compute the inverse probability censoring weights
      cases_1_dose$wt <- 1/cases_1_dose$prob
      
      # Create time data 
      cases_1_dose <- cases_1_dose %>% 
        subset(., select = c(ID, wt, dur_followup, protocol))
      
      cases_1_dose <- cbind(cases_1_dose, 
                             sim=i)
      
      cases_1_dose <- data.frame(dur_followup=seq(1, followup_period, 1)) %>% 
        left_join(., cases_1_dose, by = "dur_followup")
      
      # At each time t: Cumulative Incidence(t) = Cumulative sum of weights / Total population size 
      # This accounts for censoring by upweighting cases who were at risk of being censored
      # This gives us our cumulative incidence curve
      cases_1_dose$risk <- cumsum(tidyr::replace_na(cases_1_dose$wt, 0)) / nrow(df_1_dose)
      
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

# 4. Visualize the distribution of weights -------------------------------------
all_cases <- rbind(all_cases_2_doses, 
                   all_cases_1_dose)

# Remove rows with NA weights (time points with no events)
all_cases_clean <- all_cases %>% 
  filter(!is.na(wt))

# Randomly select 3 simulations
set.seed(123)
available_sims <- unique(all_cases_clean$sim)
selected_sims <- sample(available_sims, 3)

# Filter to only the selected simulations
all_cases_selected <- all_cases_clean %>%
  filter(sim %in% selected_sims)

# Calculate 1st and 99th percentiles per protocol and sim for the vertical lines
percentile_data <- all_cases_selected %>%
  group_by(protocol, sim) %>%
  summarise(
    p01 = quantile(wt, 0.01, na.rm = TRUE),
    p99 = quantile(wt, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

# Visualize the distribution of the inverse probability weights
# 2x3 grid: upper row for 1 dose, lower row for 2 doses, columns for each sim
ggplot(all_cases_selected, aes(x = wt)) +
  geom_histogram(fill = "steelblue", color = "white", binwidth = 0.1) +
  labs(
    title = "Distribution of Weights across Selected Simulations",
    x = "Weight (1/Probability of remaining uncensored)",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  facet_grid(protocol ~ sim, labeller = labeller(sim = function(x) paste0("sim=", x))) +
  geom_vline(data = percentile_data, aes(xintercept = p01), color = "blue", linetype = "dashed") +
  geom_vline(data = percentile_data, aes(xintercept = p99), color = "red", linetype = "dashed") +
  geom_text(data = percentile_data, aes(x = p01, label = "1%"), 
            y = Inf, vjust = 1.4, hjust = 1, color = "blue", size = 3) +
  geom_text(data = percentile_data, aes(x = p99, label = "99%"), 
            y = Inf, vjust = 1.4, hjust = -0.1, color = "red", size = 3) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5),
        panel.spacing.x = unit(1.5, "lines"))

# 5. Plot cumulative incidence -------------------------------------------------
all_data <- rbind(all_res_2_doses,
                  all_res_1_dose)

# Calculate median and 95% credible intervals by region
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

# 6. Plot relative risk -------------------------------------------------
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
  scale_y_continuous(limits = c(NA, 10)) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
