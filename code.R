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

# 2.1 Cloning & censoring ------------------------------------------------------
# We specify 100 bootstrap samples to get the confidence interval
n_boot <- 100 # 1000
boot_workers <- max(1, parallel::detectCores() - 1)
# Set a seed for reproducibility
boot_seed <- 788
end.eval <- 240

# main function
func.create.clone <- function(df, range_min, range_max) {
  df %>%
    mutate(
      clone_follow_end = case_when(
        # Received dose 2 within the protocol schedule: follow until the earliest of end of follow-up or infection
        !is.na(days.btw.doses.12) & days.btw.doses.12 >= range_min & days.btw.doses.12 <= range_max ~ pmin(days.btw.dose1.eofu, time.inf, na.rm = TRUE),
        
        # Received dose 2 earlier than the protocol schedule: follow until the earliest of dose 2 or infection
        !is.na(days.btw.doses.12) & days.btw.doses.12 < range_min ~ pmin(days.btw.doses.12, time.inf, na.rm = TRUE),
        
        # Never received dose 2 or received it later than the protocol schedule: follow until the earliest of range_max, infection, or end of follow-up
        is.na(days.btw.doses.12) | days.btw.doses.12 > range_max ~ pmin(range_max, time.inf, days.btw.dose1.eofu, na.rm = TRUE)
      ),
      clone_event = case_when(
        is.na(time.inf) ~ 0,
        
        # Infection before dose 2 and within follow-up
        !is.na(days.btw.doses.12) & time.inf < days.btw.doses.12 & time.inf <= clone_follow_end ~ 1,
        
        # Infection and dose 2 on the same day — count only if dose 2 was within the protocol schedule
        !is.na(days.btw.doses.12) & time.inf == days.btw.doses.12 & days.btw.doses.12 >= range_min & days.btw.doses.12 <= range_max ~ 1,
        
        # Infection after dose 2 but within follow-up
        !is.na(days.btw.doses.12) & time.inf > days.btw.doses.12 & time.inf <= clone_follow_end ~ 1,
        
        # Infection occurred and participant never received dose 2
        is.na(days.btw.doses.12) & time.inf <= clone_follow_end ~ 1, TRUE ~ 0
      ),
      clone_censor = ifelse(clone_event == 0 & clone_follow_end != days.btw.dose1.eofu, 1, 0)
    ) %>%
    rename(time = clone_follow_end)
}

# helper
extract_single_strata_var <- function(fm) {
  fm_chr <- paste(deparse(fm), collapse = " ")
  m <- regexpr("strata\\s*\\(\\s*(`[^`]+`|[A-Za-z\\.][A-Za-z0-9\\._]*)\\s*\\)", fm_chr, perl = TRUE)
  if (m[1] < 0) {
    return(NULL)
  }
  s <- regmatches(fm_chr, m)
  s_var <- sub("^.*strata\\s*\\(\\s*(`?)([^`)]+)\\1\\s*\\).*$", "\\2", s, perl = TRUE)
  trimws(s_var)
}

# Cox --------
# main function
run_cox <- function(fm, clone, cases, time_col = "time", mode2 = c("expected", "risk")) {
  mode2 <- match.arg(mode2)
  
  analysis_title <- paste("Cox:", paste(deparse(fm), collapse = " "))
  
  dat <- clone
  fm_local <- stats::as.formula(fm)
  environment(fm_local) <- environment()
  
  fit <- survival::coxph(fm_local, data = dat, model = FALSE, x = FALSE, y = FALSE)
  
  s_var <- extract_single_strata_var(fm_local)
  
  if (!is.null(s_var)) {
    expected <- predict(fit, newdata = cases, type = "expected")
  } else {
    bh <- survival::basehaz(fit, centered = FALSE)
    H0t <- approx(bh$time, bh$hazard, xout = cases[[time_col]], rule = 2)$y
    co <- coef(fit)
    if (mode2 == "expected") {
      if (length(co)) {
        expected <- predict(fit, newdata = cases, type = "expected")
      } else {
        expected <- H0t
      }
    } else {
      if (length(co)) {
        risk <- predict(fit, newdata = cases, type = "risk")
        expected <- H0t * risk
      } else {
        expected <- H0t
      }
    }
  }
  
  cases$prob <- exp(-expected)
  fit <- remove_environment_attr(fit)
  list(analysis_title = analysis_title, fit = fit, cases = cases)
}

# helper
remove_environment_attr <- function(fit) {
  if (!is.null(fit$terms)) {
    attr(fit$terms, ".Environment") <- NULL
  }
  if (!is.null(fit$formula)) {
    environment(fit$formula) <- emptyenv()
  }
  fit
}

# wrapper 1
est.ccw.cumrisk <- function(df, range_min, range_max, fm, full_res = FALSE, mode = "cox") {
  fm2 <- stats::as.formula(fm)
  environment(fm2) <- baseenv()
  
  clone <- func.create.clone(df, range_min, range_max)
  cases <- clone[clone$clone_event == 1, ] %>%
    select(-clone_event, -time.inf)
  
  output <- run_cox(fm2, clone, cases, mode2 = "expected")
  
  cases <- output$cases
  data.table::setorder(cases, time)
  cases$risk <- cumsum(1 / cases$prob) / nrow(df)
  
  output$res <- cases[cases$time <= end.eval, ] %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(cumrisk = max(risk, na.rm = TRUE), .groups = "drop") %>%
    complete(time = 1:end.eval) %>%
    fill(cumrisk, .direction = "down")
  
  output$cases <- cases %>%
    select(-clone_censor)
  output$range_min <- range_min
  output$range_max <- range_max
  if (full_res) {
    output$clone <- clone
  }
  output
}

# wrapper 2
run.ccw.analysis <- function(df, fm, folder_name, filename_base, protocol_mat, mode, plot_title = NULL, full_res = FALSE) {
  df <- df %>% mutate(dose2_sch = case_when(
    is.na(days.btw.doses.12) ~ protocol_mat$label[3],
    days.btw.doses.12 >= protocol_mat$min[1] & days.btw.doses.12 <= protocol_mat$max[1] ~ protocol_mat$label[1],
    days.btw.doses.12 >= protocol_mat$min[2] & days.btw.doses.12 <= protocol_mat$max[2] ~ protocol_mat$label[2]
  ))
  setDT(df)
  res.list <- lapply(seq_along(protocol_mat$protocol), function(i) {
    est.ccw.cumrisk(df, protocol_mat$min[i], protocol_mat$max[i], fm, full_res = full_res, mode = mode)
  })
  
  names(res.list) <- protocol_mat$label
  
  if (n_boot > 1) {
    out_dir <- file.path(folder_name, "bootstrap", filename_base)
    if (full_res) dir.create(out_dir, recursive = TRUE)
    
    fm_boot <- stats::as.formula(fm)
    environment(fm_boot) <- baseenv()
    
    future::plan("multisession", workers = boot_workers)
    on.exit(future::plan("sequential"), add = TRUE)
    
    globals_list <- list(
      df = df, protocol_mat = protocol_mat, fm_boot = fm_boot, full_res = full_res, mode = mode,
      end.eval = end.eval, func.create.clone = func.create.clone, est.ccw.cumrisk = est.ccw.cumrisk,
      run_cox = run_cox, extract_single_strata_var = extract_single_strata_var,
      remove_environment_attr = remove_environment_attr
    )
    
    boots_res <- future.apply::future_lapply(
      X = seq_len(n_boot),
      FUN = function(boot_id) {
        idx <- sample.int(nrow(df), size = nrow(df), replace = TRUE)
        df_b <- data.table::copy(df[idx, , drop = FALSE])
        data.table::setDT(df_b)
        
        res_b <- lapply(seq_along(protocol_mat$protocol), function(i) {
          est.ccw.cumrisk(df_b, protocol_mat$min[i], protocol_mat$max[i], fm_boot, full_res = FALSE, mode = mode)$res
        })
        names(res_b) <- protocol_mat$label
        
        #    if (full_res) saveRDS(list(idx = idx, res = res_b), file.path(out_dir, sprintf("boot-b%04d.rds", boot_id)))
        
        res_b
      },
      future.seed = boot_seed,
      future.globals = globals_list,
      future.packages = c("dplyr", "tidyr", "data.table", "survival", "survminer", "tibble", "purrr", "zoo", "splines")
    )
    
    boot_long <- purrr::imap_dfr(boots_res, function(res_b, boot_id) {
      dplyr::bind_rows(lapply(names(res_b), function(lbl) res_b[[lbl]] %>% dplyr::mutate(label = lbl, boot = as.integer(boot_id))))
    })
    
    ci_long <- boot_long %>%
      dplyr::group_by(label, time) %>%
      dplyr::summarise(
        cumrisk_lwr = stats::quantile(cumrisk, 0.025, na.rm = TRUE),
        cumrisk_med = stats::median(cumrisk, na.rm = TRUE),
        cumrisk_upr = stats::quantile(cumrisk, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    for (lbl in names(res.list)) {
      res.list[[lbl]]$res <- res.list[[lbl]]$res %>%
        dplyr::left_join(
          ci_long %>% dplyr::filter(label == lbl),
          by = "time"
        )
    }
  }
  
  res.list
}

# wrapper 3
run.ccw.analysis_with_formulas <- function(formulas, protocol_mat, folder_name, mode = "cox", full_res = FALSE, sample_prop = 1) {
  
  dat <- read_csv("https://github.com/katjia/ISPPD_TTE_workshop/raw/main/dat.csv") %>%
    slice_sample(prop = sample_prop)
  results_all <- list()
  for (fm_name in names(formulas)) {
    plot_title <- paste("Censoring model covariates:", paste(all.vars(formulas[[fm_name]][[3]]), collapse = ", "))
    results_all[[fm_name]] <- run.ccw.analysis(dat, formulas[[fm_name]], folder_name, paste("cumrisk", fm_name, mode, sep = "-"), protocol_mat, mode, plot_title = plot_title, full_res = full_res)
  }
  
  dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    results_all,
    file = file.path(folder_name, paste0(gsub("\\.", "", paste0("res_", mode)), ".rds"))
  )
  
  invisible(NULL)
}

# Specification -------

set.seed(123456)

# Read Data
data_date <- "20260424"

# CHANGE HERE
# Replace with your specific Raw URL
github_url <- "https://github.com/katjia/ISPPD_TTE_workshop/raw/main/data-gen_run_details.RData"
load(url(github_url))

## Analysis Protocols
protocol_analysis_fixed <- list(
  protocol = protocol_data$protocol[2:4],
  label = protocol_data$protocol[2:4],
  min = protocol_data$min[2:4],
  max = c(protocol_data$max[2], protocol_data$max[3], Inf),
  col = c("dodgerblue", "goldenrod", "hotpink")[1:3]
)
protocol_analysis_fixed$label[protocol_analysis_fixed$protocol == "recommended"] <- "shorter"
protocol_analysis_fixed$label[protocol_analysis_fixed$protocol == "allowable"] <- "two doses"
protocol_analysis_fixed$label[protocol_analysis_fixed$protocol == "late"] <- "1 dose"
protocol_analysis_fixed$label_fig <- tools::toTitleCase(protocol_analysis_fixed$label)
protocol_analysis_fixed$label_fig[protocol_analysis_fixed$protocol == "recommended"] <- "Shorter Interval"
protocol_analysis_fixed$label_fig[protocol_analysis_fixed$protocol == "allowable"] <- "two doses"
protocol_analysis_fixed$label_fig[protocol_analysis_fixed$protocol == "late"] <- "1 dose"
protocol_analysis_fixed$min <- round((protocol_analysis_fixed$min + protocol_analysis_fixed$max) / 2)
protocol_analysis_fixed$max <- round((protocol_analysis_fixed$min + protocol_analysis_fixed$max) / 2)
rm(protocol_data)

rm(all_nodes, data_seed, dose1_prop_0, dose1_prop_1, group1_prop, infection_log_combined, node_vax_counts, initial_SEIR, phase_lookup, vax_detailed_final, vax_phase_group_plan)

## Formula List
formulas_without_group_month <- list(
  inf = as.formula("survival::Surv(time, clone_censor) ~ prior.inf"),
  no_cov = as.formula("survival::Surv(time, clone_censor) ~ 1")
)

xlab_from_vax <- "Days since first dose"

ylab_cumrisk <- "Cumulative risk, %"

ylim_vec_cumrisk <- c(0, 0.285)
ylim_vec_rr <- c(0.5, 25)
scale_y_continuous_percent <- function(ylim_vec = c(0, NA), breaks = waiver(), n.breaks = NULL, expand_mult = c(0, 0)) {
  scale_y_continuous(limits = ylim_vec, breaks = breaks, n.breaks = n.breaks, labels = function(y) paste0(y * 100), expand = expansion(mult = expand_mult))
}

scale_y_continuous_rr_log <- function(ylim_vec = ylim_vec_rr) {
  list(
    scale_y_continuous(trans = "log", breaks = 2^seq(-1, 5, 1), limits = ylim_vec, expand = expansion(mult = c(0, 0))),
    labs(x = xlab_from_vax, y = "Risk ratio (vs. shorter)")
  )
}

scale_x_continuous_days <- function(break_by = 60) {
  scale_x_continuous(breaks = seq(0, 1000, break_by), minor_breaks = seq(0, 1000, break_by / 2), expand = expansion(mult = c(0, 0)))
}

# Run -------
set.seed(123456)

single_time_noGrp_analysis <- function(prefix, comment, dose1_prop_char, folder_name = paste0(comment, "_", data_date)) {
  run.ccw.analysis_with_formulas(formulas_without_group_month, protocol_analysis_fixed, folder_name, mode = "cox")
  invisible(NULL)
}

single_time_noGrp_analysis("240_20_", "240_20_multi_noGrp")

# Here we visualize weights for the 'recommended' protocol
results_all <- readRDS(
  file = file.path("240_20_multi_noGrp_20260424", paste0(gsub("\\.", "", paste0("res_", "cox")), ".rds")))

res_data <- results_all$no_cov$shorter$cases

# Visualize the distribution of the inverse probability weights
ggplot(res_data, aes(x = 1/prob)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Inverse Probability Weights (IPW)",
    subtitle = "Protocol: Two doses",
    x = "Weight (1/Probability of remaining uncensored)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  geom_vline(aes(xintercept = quantile(1/prob, 0.99)), color = "red", linetype = "dashed") +
  annotate("text", x = quantile(1/res_data$prob, 0.99), y = 0,
           label = "99th Percentile", vjust = -1, color = "red", angle = 90)

res_data_2 <- results_all$no_cov$`1 dose`$cases

# Visualize the distribution of the inverse probability weights
ggplot(res_data_2, aes(x = 1/prob)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Inverse Probability Weights (IPW)",
    subtitle = "Protocol: Dose 1 only",
    x = "Weight (1/Probability of remaining uncensored)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  geom_vline(aes(xintercept = quantile(1/prob, 0.99)), color = "red", linetype = "dashed") +
  annotate("text", x = quantile(1/res_data$prob, 0.99), y = 0,
           label = "99th Percentile", vjust = -1, color = "red", angle = 90)


# helper function from `2-0_tte_func.R` ------------------
res_list_formula_long <- function(res_list_formula, res_col = "res") {
  res_long <- purrr::map_dfr(names(res_list_formula), function(nm) {
    res_list_formula[[nm]][[res_col]] %>%
      dplyr::mutate(group = nm)
  })
  
  res_long
}

res_list_formula_wide <- function(res_list_formula,
                                  res_col = "res",
                                  value_col = "cumrisk") {
  res_long <- res_list_formula_long(res_list_formula, res_col = res_col)
  
  res_wide <- res_long %>%
    dplyr::select(time, group, all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = group,
                       values_from = !!rlang::sym(value_col))
  
  res_wide
}

legend_labels_range_str <- function(p, rmin, rmax) {
  range_str <- if (rmin > end.eval) {
    ""
  } else if (is.infinite(rmax)) {
    paste0(" (", rmin, "+ days)")
  } else if (rmin == rmax) {
    paste0(" (", rmin, " days)")
  } else {
    paste0(" (", rmin, "-", rmax, " days)")
  }
  
  paste0(p, range_str)
}

protocol_legend_labels <- function(protocol_mat, res_map) {
  sapply(seq_along(protocol_mat$protocol), function(i) {
    legend_labels_range_str(protocol_mat$label_fig[i],
                            res_map[[i]]$range_min,
                            res_map[[i]]$range_max)
  })
}

ggplot_cumrisk_core <- function(res_1,
                                res_2,
                                protocol_mat,
                                dataset_labels = c("CCW", "Counterfactual"),
                                y_var = c("cumrisk", "rr")) {
  y_var <- match.arg(y_var)
  
  if (y_var == "rr") {
    add_rr_to_res <- function(res) {
      if (is.null(res)) {
        return(NULL)
      }
      risk_table <- res_list_formula_wide(res)
      for (j in seq_len(length(res))) {
        res[[j]]$res <- res[[j]]$res %>%
          dplyr::left_join(risk_table, by = dplyr::join_by(time))
        res[[j]]$res$rr <- res[[j]]$res$cumrisk / res[[j]]$res$shorter
      }
      res
    }
    res_1 <- add_rr_to_res(res_1)
    res_2 <- add_rr_to_res(res_2)
  }
  
  make_df <- function(res_map, dataset) {
    bind_rows(lapply(protocol_mat$label, function(lbl) {
      df <- res_map[[lbl]]$res
      df$protocol <- lbl
      df$dataset <- dataset
      df
    }))
  }
  df1 <- make_df(res_1, dataset_labels[1])
  if (!is.null(res_2)) {
    df2 <- make_df(res_2, dataset_labels[2])
    df <- bind_rows(df1, df2)
  } else {
    df <- df1
  }
  
  col_vals <- setNames(protocol_mat$col, protocol_mat$label)
  proto_labels <- setNames(protocol_legend_labels(protocol_mat, res_1),
                           protocol_mat$label)
  
  p <- ggplot(df, aes(x = time, y = .data[[y_var]], color = protocol))
  
  if (y_var == "cumrisk" &
      all(c("cumrisk_lwr", "cumrisk_upr") %in% names(df1))) {
    df1_ci <- df1 %>% filter(is.finite(cumrisk_lwr), is.finite(cumrisk_upr))
    
    p <- p +
      geom_ribbon(
        data = df1_ci,
        aes(
          x = time,
          ymin = cumrisk_lwr,
          ymax = cumrisk_upr,
          fill = protocol,
          group = protocol
        ),
        inherit.aes = FALSE,
        alpha = 0.18,
        color = NA
      ) +
      scale_fill_manual(values = col_vals) + guides(fill = "none")
  }
  
  p <- p +
    geom_line(aes(linetype = dataset, linewidth = dataset)) +
    scale_linewidth_manual(values = c(CCW = 0.2, Counterfactual = 0.35),
                           guide = "none") +
    scale_x_continuous_days() +
    scale_color_manual(
      values = col_vals,
      breaks = protocol_mat$label,
      labels = proto_labels,
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed"),
      breaks = dataset_labels,
      labels = dataset_labels,
      name = NULL
    ) +
    labs(x = xlab_from_vax, y = ylab_cumrisk) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      plot.tag = element_text(
        size = 16,
        face = "bold",
        hjust = 0,
        vjust = 0
      ),
      plot.tag.position = c(0.05, 1.05),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(6, "mm"),
      legend.spacing.y = unit(0.5, "pt"),
      legend.spacing.x = unit(1, "pt"),
      legend.box.spacing = unit(4, "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      plot.margin = margin(
        t = 18,
        r = 3,
        b = 1,
        l = 1,
        unit = "mm"
      )
    ) + guides()
  
  if (!is.null(res_2)) {
    p <- p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
  } else {
    p <- p + guides(color = guide_legend(nrow = 1), linetype = "none")
  }
  
  if (y_var == "rr") {
    p <- p + scale_y_continuous_rr_log()
  } else if (y_var == "cumrisk") {
    p <- p + scale_y_continuous_percent(ylim_vec = ylim_vec_cumrisk)
  }
  
  p
}

# `source("3-0_res_func.R")` -----------------------------

library(patchwork)

add_fig_tag <- function(fig_list,
                        fig_tag_text_list = rep("", length(fig_list)),
                        panel_labels = LETTERS) {
  lapply(seq_along(fig_list), function(i) {
    fig_list[[i]] +
      labs(tag = paste0(panel_labels[i], ") ", fig_tag_text_list[i]))
  })
}

ggplot_cumrisk_overlay <- function(fm_name,
                                   mode,
                                   prefix_1,
                                   prefix_2,
                                   protocol_mat,
                                   folder_name_1 = paste0(prefix_1, "_", data_date),
                                   folder_name_2 = paste0(prefix_2, protocol_mat$protocol[1], "_", data_date),
                                   y_var = "cumrisk") {
  res_1 <- readRDS(file.path(folder_name_1, paste0("res_", mode, ".rds")))[[fm_name]]
  # res_2 <- readRDS(file.path(folder_name_2, "res_counterfactual.rds"))
  ggplot_cumrisk_core(res_1, res_2 = NULL, protocol_mat, y_var = y_var)
}

# Filter protocols to keep only "late", "recommended", and "unvax" ----------
# Note: "unvax" will be added later, for now keeping "late" and "recommended"
protocol_analysis_filtered <- lapply(protocol_analysis_fixed, function(x) {
  if (is.null(names(x))) {
    keep_idx <- protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")
    x[keep_idx]
  } else {
    x
  }
})
protocol_analysis_filtered$protocol <- protocol_analysis_fixed$protocol[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]
protocol_analysis_filtered$label <- protocol_analysis_fixed$label[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]
protocol_analysis_filtered$min <- protocol_analysis_fixed$min[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]
protocol_analysis_filtered$max <- protocol_analysis_fixed$max[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]
protocol_analysis_filtered$col <- protocol_analysis_fixed$col[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]
protocol_analysis_filtered$label_fig <- protocol_analysis_fixed$label_fig[protocol_analysis_fixed$protocol %in% c("late", "recommended", "unvax")]

# plot cum risk ----------

width_image_cm <- 15
height_image_cm <- 8

p_vac <- seq(20, 60, 20)
p_unvac <- 100 - p_vac

fig_tag_text_list_time1 <-
  paste0(p_unvac, "% unvaccinated; ", p_vac, "% received ≥1 dose")

ylim_vec_cumrisk <- c(0, 0.4)

list_plot_cumrisk_overlay <- list()
list_plot_cumrisk_overlay[[1]] <- ggplot_cumrisk_overlay("no_cov",
                                                         "cox",
                                                         "240_20_multi_noGrp",
                                                         "240_20_",
                                                         protocol_analysis_filtered)
list_plot_cumrisk_overlay <- add_fig_tag(list_plot_cumrisk_overlay, fig_tag_text_list = fig_tag_text_list_time1)


plot_cumrisk_overlay_combined <- wrap_plots(list_plot_cumrisk_overlay,
                                            ncol = 1,
                                            guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
plot_cumrisk_overlay_combined

ggsave(
  file.path(
    paste0("240_20_multi_noGrp", "_", data_date),
    "plot_cumrisk_overlay_no_cov-km.tiff"
  ),
  plot_cumrisk_overlay_combined,
  width = width_image_cm,
  height = height_image_cm * 1,
  units = "cm",
  dpi = 1800,
  compression = "lzw"
)

# plot RR ----------

ylim_vec_rr <- c(0.3, 26)
list_plot_rr_overlay <- list()
list_plot_rr_overlay[[1]] <- ggplot_cumrisk_overlay(
  "no_cov",
  "cox",
  "240_20_multi_noGrp",
  "240_20_",
  protocol_analysis_filtered,
  y_var = "rr"
)
list_plot_rr_overlay <- add_fig_tag(list_plot_rr_overlay, fig_tag_text_list = fig_tag_text_list_time1)

plot_rr_overlay_combined <- wrap_plots(list_plot_rr_overlay, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
plot_rr_overlay_combined

ggsave(
  file.path(
    paste0("240_20_multi_noGrp", "_", data_date),
    "plot_rr_overlay_no_cov-km.tiff"
  ),
  plot_rr_overlay_combined,
  width = width_image_cm,
  height = height_image_cm * 1,
  units = "cm",
  dpi = 1800,
  compression = "lzw"
)


