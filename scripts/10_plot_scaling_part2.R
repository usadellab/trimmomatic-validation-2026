#!/usr/bin/env Rscript

# ---
# R Script: Large Scale Stress Test (1 - 240 Threads) + Accuracy
# Plots: Wall Time, RAM Usage, Speedup Factor, and Verified Accuracy
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
LOG_DIR <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_extended_v2"
OUTPUT_BASE <- file.path(LOG_DIR, "260320_Figure_Scaling_Benchmark")

RESIDUALS_PLANT <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/residuals_detailed_plant/verified_residuals.txt"
RESIDUALS_HUMAN <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/residuals_detailed_human/verified_residuals.txt"

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
    # Extract Metadata: [Dataset]_[Tool]_[Threads]t.time.log
    parts <- str_match(fname, "^([A-Za-z]+)_([A-Za-z0-9]+)_([0-9]+)t\\.time\\.log$")
    
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
      
      # 3. Check RabbitTrim compression tools
      label <- ""
      if (tolower(tool) == "rabbittrim") {
        pragzip_zero <- any(grepl("pragzip[^0-9]*\\b0\\b|rapidgzip[^0-9]*\\b0\\b", lines, ignore.case = TRUE))
        pigz_zero <- any(grepl("pigz[^0-9]*\\b0\\b", lines, ignore.case = TRUE))
        
        if (pragzip_zero && pigz_zero) {
          label <- "*"   # pragzip and pigz are not used
        }
      }
      
      results <- rbind(results, data.frame(
        Tool = tool, Dataset = dataset, Threads = threads, 
        Time_Min = t_min, Memory_GB = mem_gb, Label = label,
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
cat("Parsing scaling logs from:", LOG_DIR, "\n")
df_perf <- parse_scaling_logs(LOG_DIR)

standardize_tool <- function(t) {
  t_low <- tolower(t)
  case_when(
    t_low == "trimmomatic" ~ "Trimmomatic",
    t_low == "rabbittrim" ~ "RabbitTrim",
    t_low == "fastp" ~ "fastp",
    t_low == "bbduk" ~ "BBDuk",
    t_low == "skewer" ~ "Skewer",
    t_low == "cutadapt" ~ "Cutadapt",
    TRUE ~ t
  )
}

df_perf$Tool <- sapply(df_perf$Tool, standardize_tool)
df_perf$Dataset <- factor(df_perf$Dataset, levels = c("Plant", "Human"))

tool_levels <- c("Trimmomatic", "RabbitTrim", "fastp", "BBDuk", "Skewer", "Cutadapt")
df_perf$Tool <- factor(df_perf$Tool, levels = tool_levels)

# Calculate Speedup (Baseline = 1 Thread)
df_speedup <- df_perf %>%
  group_by(Tool, Dataset) %>%
  mutate(
    Base_Time = if(any(Threads == 1)) Time_Min[Threads == 1] else NA,
    Speedup = Base_Time / Time_Min
  ) %>%
  ungroup() %>%
  filter(!is.na(Speedup))

cat("Parsing accuracy metrics...\n")
df_acc_plant <- parse_residuals(RESIDUALS_PLANT, "Plant")
df_acc_human <- parse_residuals(RESIDUALS_HUMAN, "Human")
df_acc <- rbind(df_acc_plant, df_acc_human)

if(nrow(df_acc) > 0) {
  df_acc$Tool <- sapply(df_acc$Tool, standardize_tool)
  df_acc$Tool <- factor(df_acc$Tool, levels = tool_levels)
  df_acc$Dataset <- factor(df_acc$Dataset, levels = c("Plant", "Human"))
}

# --- Plotting ---

tool_colors <- c(
  "Trimmomatic" = "#E41A1C", "RabbitTrim" = "#377EB8", 
  "fastp" = "#4DAF4A", "BBDuk" = "#984EA3",
  "Skewer" = "#FF7F00", "Cutadapt" = "#E6C200"
)

# Shared Legend Title containing the Explanation
legend_title <- "Tool   (RabbitTrim: [*] pragzip & pigz not used)"

my_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Define X-axis spanning the complete tested range
x_breaks <- c(1, 2, 4, 8, 16, 20, 40, 80, 160, 240)
ideal_scaling_df <- data.frame(Threads = seq(1, 240, length.out = 100))
ideal_scaling_df$Speedup <- ideal_scaling_df$Threads / 1

# Plot A: Wall Time
p1 <- ggplot(df_speedup, aes(x = Threads, y = Time_Min, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 4, color = "black", show.legend = FALSE) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_color_manual(name = legend_title, values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "A. Wall Clock Time", y = "Time (min)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot B: RAM Usage
p2 <- ggplot(df_speedup, aes(x = Threads, y = Memory_GB, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 4, color = "black", show.legend = FALSE) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_color_manual(name = legend_title, values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "B. Peak Memory Usage", y = "RAM (GB)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot C: Speedup Factor
# Determine a reasonable dynamic cap for the Y-axis so the ideal line doesn't compress the plot
max_actual_speedup <- max(df_speedup$Speedup, na.rm = TRUE)
plot_c_y_limit <- ceiling(max_actual_speedup / 5) * 5 + 5 # Round up to nearest 5 and add 5 buffer

p3 <- ggplot(df_speedup, aes(x = Threads, y = Speedup, color = Tool, group = Tool)) +
  geom_line(data = ideal_scaling_df, aes(x = Threads, y = Speedup), 
            color = "gray50", linetype = "dashed", inherit.aes = FALSE) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 4, color = "black", show.legend = FALSE) +
  facet_wrap(~Dataset) +
  scale_color_manual(name = legend_title, values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  coord_cartesian(ylim = c(0, plot_c_y_limit)) + 
  annotate("text", x = 16, y = 18, label = "Ideal Scaling", color = "gray50", angle = 35, size = 3.5) +
  labs(title = "C. Scaling Factor (vs 1 Thread)", y = "Fold Speedup", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot D: Verified Accuracy
if(nrow(df_acc) > 0) {
  p4 <- ggplot(df_acc, aes(x = Tool, y = True_Residuals, fill = Tool)) +
    geom_col(color="black", size=0.2) +
    geom_text(aes(label = scales::comma(True_Residuals)), vjust = -0.5, size=3.5) +
    facet_wrap(~Dataset, scales = "free_y") +
    scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                       breaks = c(0, 10, 100, 1000, 10000, 100000), 
                       labels = scales::comma, 
                       expand = expansion(mult = c(0, 0.2))) +
    scale_fill_manual(values = tool_colors) +
    guides(fill = "none") + 
    labs(title = "D. Verified Accuracy", 
         y = "Total Residual Adapter Reads", x = NULL) +
    my_theme + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    
  layout <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
} else {
  warning("No valid residual files found. Plot D will not be generated.")
  layout <- (p1 + p2) / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

# Save
ggsave(paste0(OUTPUT_BASE, ".pdf"), layout, width = 14, height = 12)
ggsave(paste0(OUTPUT_BASE, ".png"), layout, width = 14, height = 12, dpi = 300)

cat("\nPlots saved to:", paste0(OUTPUT_BASE, ".png"), "\n")