#!/usr/bin/env Rscript

# ---
# R Script: Scaling Stress Test (20 vs 80 Threads)
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

# --- Configuration ---
LOG_DIR <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling"
OUTPUT_BASE <- file.path(LOG_DIR, "Figure_Scaling_StressTest")

# --- Parsing Function ---
parse_scaling_logs <- function(dir) {
  # List all time logs in the directory
  files <- list.files(dir, pattern = "_[0-9]+t\\.time\\.log$", full.names = TRUE)
  
  results <- data.frame(
    Tool = character(),
    Dataset = character(),
    Threads = numeric(),
    Time_Min = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (length(files) == 0) {
    stop("No logs found in ", dir, "! Make sure the path is correct and logs exist.")
  }
  
  for (fpath in files) {
    fname <- basename(fpath)
    
    # Extract Metadata from filename: e.g., "Plant_Trimmomatic_80t.time.log"
    # Pattern: [Dataset]_[Tool]_[Threads]t.time.log
    parts <- str_match(fname, "^([A-Za-z]+)_([A-Za-z]+)_([0-9]+)t\\.time\\.log$")
    
    if (!any(is.na(parts[1, ]))) {
      dataset <- parts[1, 2]
      tool <- parts[1, 3]
      threads <- as.numeric(parts[1, 4])
      
      # Extract Time from file content
      lines <- readLines(fpath)
      time_line <- grep("Elapsed (wall clock) time", lines, value = TRUE, fixed = TRUE)
      
      t_min <- NA
      if (length(time_line) > 0) {
        time_str <- str_extract(tail(time_line, 1), "[0-9:.]+$")
        time_parts <- as.numeric(unlist(str_split(time_str, ":")))
        
        # Convert to Minutes
        if (length(time_parts) == 3) { # h:mm:ss
          t_min <- time_parts[1]*60 + time_parts[2] + time_parts[3]/60
        } else if (length(time_parts) == 2) { # m:ss
          t_min <- time_parts[1] + time_parts[2]/60
        }
      }
      
      results <- rbind(results, data.frame(Tool = tool, Dataset = dataset, Threads = threads, Time_Min = t_min))
    }
  }
  return(results)
}

# --- Main Logic ---

cat("Parsing scaling logs from:", LOG_DIR, "\n")
df <- parse_scaling_logs(LOG_DIR)

# Sort Factor Levels
df$Tool <- factor(df$Tool, levels = c("Trimmomatic", "RabbitTrim", "fastp", "BBDuk"))
df$Dataset <- factor(df$Dataset, levels = c("Plant", "Human"))
df$Threads_Label <- paste0(df$Threads, " Threads") # For plot labels

# Calculate Scaling Factor (Speedup)
# Speedup = Time_20 / Time_80
scaling_df <- df %>%
  group_by(Tool, Dataset) %>%
  arrange(Threads) %>%
  summarise(
    Time_20 = Time_Min[Threads == 20],
    Time_80 = Time_Min[Threads == 80],
    Speedup = Time_20 / Time_80,
    .groups = "drop"
  ) %>%
  filter(!is.na(Speedup))

print(scaling_df) # Show table in console

# --- Plotting ---

# Colors for Thread Counts
thread_colors <- c("20 Threads" = "#999999", "80 Threads" = "#E69F00")

# Theme
my_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Plot A: Raw Time Comparison
p1 <- ggplot(df, aes(x = Tool, y = Time_Min, fill = factor(paste(Threads, "Threads")))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color="black", size=0.2) +
  geom_text(aes(label = sprintf("%.1f m", Time_Min)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_fill_manual(values = thread_colors, name = "") +
  labs(title = "A. Raw Processing Time (20 vs 80 Threads)", y = "Wall Time (minutes)", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  my_theme

# Plot B: Scaling Efficiency (Speedup Factor)
# Ideal scaling (4x threads) would be 4.0x speedup
p2 <- ggplot(scaling_df, aes(x = Tool, y = Speedup, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color="black", size=0.2) +
  geom_hline(yintercept = 4.0, linetype = "dashed", color = "red") + # Theoretical max
  geom_text(aes(label = sprintf("%.1fx", Speedup)), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  annotate("text", x = 0.6, y = 4.1, label = "Ideal Scaling (4x)", color = "red", size = 3, hjust=0) +
  scale_fill_manual(values = c("Plant" = "#009E73", "Human" = "#0072B2")) +
  labs(title = "B. Observed Speedup (Going from 20 -> 80 Threads)", 
       subtitle = "Calculated as: Time(20t) / Time(80t)",
       y = "Fold Speedup", x = NULL) +
  scale_y_continuous(limits = c(0, 4.5), expand = expansion(mult = c(0, 0.05))) +
  my_theme

# Combine
final_plot <- p1 / p2 + plot_layout(heights = c(1, 0.8))

# Save
ggsave(paste0(OUTPUT_BASE, ".pdf"), final_plot, width = 10, height = 10)
ggsave(paste0(OUTPUT_BASE, ".png"), final_plot, width = 10, height = 10, dpi = 300)

cat("\nPlot saved to:", paste0(OUTPUT_BASE, ".png"), "\n")