#!/usr/bin/env Rscript

# ---
# R Script: Fastp Variant Comparison vs Trimmomatic
# Plots: Wall Time, RAM Usage, and Verified Accuracy
# ---

# Check/Install packages
packages <- c("ggplot2", "dplyr", "readr", "stringr", "patchwork", "scales")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)
library(scales)

# --- Configuration ---
DIR_V2 <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling"
DIR_F11_NO <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_fastp-v1.1.0-no-adapters"
DIR_F13_WITH <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_fastp-v1.3.3"
DIR_F13_NO <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_fastp-v1.3.3-no-adapters"

OUTPUT_BASE <- file.path(dirname(DIR_V2), "Figure_Fastp_Comparison")

OLD_RESIDUALS_PLANT <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/residuals_detailed_plant/verified_residuals.txt"
OLD_RESIDUALS_HUMAN <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/residuals_detailed_human/verified_residuals.txt"
NEW_RESIDUALS_PLANT <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling_fastp-v1.1.0-no-adapters/residuals_detailed_plant/verified_residuals.txt"
NEW_RESIDUALS_HUMAN <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling_fastp-v1.1.0-no-adapters/residuals_detailed_human/verified_residuals.txt"

# --- Parser 1: Time & Performance Logs ---
parse_scaling_logs <- function(dir) {
  files <- list.files(dir, pattern = "_[0-9]+t\\.time\\.log$", full.names = TRUE)
  
  results <- data.frame()
  if (length(files) == 0) {
    warning("No time logs found in ", dir)
    return(results)
  }
  
  for (fpath in files) {
    fname <- basename(fpath)
    # Extract Metadata
    parts <- str_match(fname, "^([A-Za-z]+)_([A-Za-z0-9.-]+)_([0-9]+)t\\.time\\.log$")
    
    if (!any(is.na(parts[1, ]))) {
      dataset <- parts[1, 2]
      tool <- parts[1, 3]
      threads <- as.numeric(parts[1, 4])
      
      lines <- readLines(fpath, warn = FALSE)
      
      # 1. Parse Time
      time_line <- grep("Elapsed \\(wall clock\\) time", lines, value = TRUE, ignore.case = TRUE)
      t_min <- NA
      if (length(time_line) > 0) {
        time_str <- str_extract(tail(time_line, 1), "[0-9:.]+$")
        time_parts <- as.numeric(unlist(str_split(time_str, ":")))
        if (length(time_parts) == 3) { 
          t_min <- time_parts[1]*60 + time_parts[2] + time_parts[3]/60 
        } else if (length(time_parts) == 2) { 
          t_min <- time_parts[1] + time_parts[2]/60 
        }
      }
      
      # 2. Parse Memory
      mem_line <- grep("Maximum resident set size", lines, value = TRUE, ignore.case = TRUE)
      mem_gb <- NA
      if (length(mem_line) > 0) {
        mem_kb <- as.numeric(str_extract(tail(mem_line, 1), "\\d+"))
        mem_gb <- mem_kb / 1024 / 1024
      }
      
      results <- rbind(results, data.frame(
        Tool = tool, Dataset = dataset, Threads = threads, 
        Time_Min = t_min, Memory_GB = mem_gb, Label = "",
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}

# --- Parser 2: Verified Residuals Accuracy ---
parse_residuals <- function(filepath, dataset_name) {
  if (!file.exists(filepath)) return(data.frame())
  
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[grepl("\\|", lines)]
  lines <- lines[!grepl("Tool.*Read", lines, ignore.case = TRUE) & !grepl("---", lines)]
  
  if(length(lines) == 0) return(data.frame())
  
  df <- data.frame()
  for (l in lines) {
    parts <- str_trim(unlist(str_split(l, "\\|")))
    if (length(parts) >= 4) {
      tool <- parts[1]
      true_res <- as.numeric(gsub("[^0-9]", "", parts[4]))
      if (!is.na(true_res)) {
        df <- rbind(df, data.frame(Tool = tool, True_Residuals = true_res, stringsAsFactors = FALSE))
      }
    }
  }
  
  if(nrow(df) > 0) {
      df <- df %>% 
        group_by(Tool) %>% 
        summarise(True_Residuals = sum(True_Residuals), .groups="drop")
      df$Dataset <- dataset_name
  }
  return(df)
}

# --- Main Logic ---
cat("Parsing scaling logs...\n")

# Load and explicitly assign names depending on directory origin
df_v2 <- parse_scaling_logs(DIR_V2)
df_trim <- df_v2 %>% filter(grepl("trimmomatic", tolower(Tool))) %>% mutate(Tool = "Trimmomatic v0.40")
df_f11_with <- df_v2 %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.1.0 (with-adapters)")

df_f11_no <- parse_scaling_logs(DIR_F11_NO) %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.1.0 (no-adapters)")
df_f13_with <- parse_scaling_logs(DIR_F13_WITH) %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.3.3 (with-adapters)")
df_f13_no <- parse_scaling_logs(DIR_F13_NO) %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.3.3 (no-adapters)")

# Combine and subset allowed threads
df_perf <- bind_rows(df_trim, df_f11_with, df_f11_no, df_f13_with, df_f13_no) %>%
  filter(Threads %in% c(1, 2, 4, 8, 16, 40, 80, 160, 240))

tool_levels <- c("Trimmomatic v0.40", 
                 "fastp v1.1.0 (with-adapters)", "fastp v1.1.0 (no-adapters)",
                 "fastp v1.3.3 (with-adapters)", "fastp v1.3.3 (no-adapters)")

df_perf$Dataset <- factor(df_perf$Dataset, levels = c("Plant", "Human"))
df_perf$Tool <- factor(df_perf$Tool, levels = tool_levels)

cat("Parsing accuracy metrics...\n")
df_acc_old_plant <- parse_residuals(OLD_RESIDUALS_PLANT, "Plant")
df_acc_old_human <- parse_residuals(OLD_RESIDUALS_HUMAN, "Human")
df_acc_new_plant <- parse_residuals(NEW_RESIDUALS_PLANT, "Plant")
df_acc_new_human <- parse_residuals(NEW_RESIDUALS_HUMAN, "Human")

df_acc_old <- rbind(df_acc_old_plant, df_acc_old_human)
df_acc_new <- rbind(df_acc_new_plant, df_acc_new_human)

# Duplicate/Rename correctly for accuracy files
df_acc_trim <- df_acc_old %>% filter(grepl("trimmomatic", tolower(Tool))) %>% mutate(Tool = "Trimmomatic v0.40")
df_acc_f11_with <- df_acc_old %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.1.0 (with-adapters)")
df_acc_f13_with <- df_acc_old %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.3.3 (with-adapters)")

df_acc_f11_no <- df_acc_new %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.1.0 (no-adapters)")
df_acc_f13_no <- df_acc_new %>% filter(grepl("fastp", tolower(Tool))) %>% mutate(Tool = "fastp v1.3.3 (no-adapters)")

df_acc <- bind_rows(df_acc_trim, df_acc_f11_with, df_acc_f13_with, df_acc_f11_no, df_acc_f13_no)

if(nrow(df_acc) > 0) {
  # Re-summarize to combine R1/R2 rows that were collapsed into the same tool name
  df_acc <- df_acc %>% 
    group_by(Tool, Dataset) %>% 
    summarise(True_Residuals = sum(True_Residuals), .groups="drop")

  df_acc$Tool <- factor(df_acc$Tool, levels = tool_levels)
  df_acc$Dataset <- factor(df_acc$Dataset, levels = c("Plant", "Human"))
}

# --- Plotting ---

tool_colors <- c(
  "Trimmomatic v0.40" = "#E41A1C", 
  "fastp v1.1.0 (with-adapters)" = "#377EB8", 
  "fastp v1.1.0 (no-adapters)" = "#984EA3",
  "fastp v1.3.3 (with-adapters)" = "#4DAF4A",
  "fastp v1.3.3 (no-adapters)" = "#FF7F00"
)

my_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

x_breaks <- c(1, 2, 4, 8, 16, 40, 80, 160, 240)

# Plot A: Wall Time
p1 <- ggplot(df_perf, aes(x = Threads, y = Time_Min, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_color_manual(name = "Tool", values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "A. Wall Clock Time", y = "Time (min)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot B: RAM Usage
p2 <- ggplot(df_perf, aes(x = Threads, y = Memory_GB, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_color_manual(name = "Tool", values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "B. Peak Memory Usage", y = "RAM (GB)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot C: Verified Accuracy
if(nrow(df_acc) > 0) {
  p3 <- ggplot(df_acc, aes(x = Tool, y = True_Residuals, fill = Tool)) +
    geom_col(color="black", size=0.2) +
    geom_text(aes(label = scales::comma(True_Residuals)), vjust = -0.5, size=3.5) +
    facet_wrap(~Dataset, scales = "free_y") +
    scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                       breaks = c(0, 10, 100, 1000, 10000, 100000), 
                       labels = scales::comma, 
                       expand = expansion(mult = c(0, 0.2))) +
    scale_fill_manual(values = tool_colors) +
    guides(fill = "none") + 
    labs(title = "C. Verified Accuracy", 
         y = "Total Residual Adapter Reads", x = NULL) +
    my_theme + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    
  layout <- (p1 + p2) / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
} else {
  warning("No valid residual files found. Plot C will not be generated.")
  layout <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

# Save
ggsave(paste0(OUTPUT_BASE, ".pdf"), layout, width = 14, height = 12)
ggsave(paste0(OUTPUT_BASE, ".png"), layout, width = 14, height = 12, dpi = 300)

cat("\nPlots saved to:", paste0(OUTPUT_BASE, ".png"), "\n")