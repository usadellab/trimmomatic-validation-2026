#!/usr/bin/env Rscript

# ---
# R Script: Benchmark on RabbitTrim Paper Dataset (SRR7890824)
# Plots: Wall Time and RAM Usage for Trimmomatic vs RabbitTrim
# ---

# Check/Install packages
packages <- c("ggplot2", "dplyr", "readr", "stringr", "patchwork", "scales")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(scales)

# --- Configuration ---
LOG_DIR <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/rabbit_scaling"
FIG_DIR <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/figures"
OUTPUT_BASE <- file.path(FIG_DIR, "Figure_RabbitTrim_Dataset_Benchmark")

# --- Parser Function ---
parse_scaling_logs <- function(dir) {
  # Matches the filename structure: SRR7890824_Tool_8t_c6.time.log
  files <- list.files(dir, pattern = "_[0-9]+t_c6\\.time\\.log$", full.names = TRUE)
  
  results <- data.frame()
  if (length(files) == 0) {
    warning("No time logs found in ", dir)
    return(results)
  }
  
  for (fpath in files) {
    fname <- basename(fpath)
    # Extract Metadata
    parts <- str_match(fname, "^SRR7890824_([A-Za-z0-9]+)_([0-9]+)t_c6\\.time\\.log$")
    
    if (!any(is.na(parts[1, ]))) {
      tool <- parts[1, 2]
      threads <- as.numeric(parts[1, 3])
      
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
      
      # 3. Check RabbitTrim compression tools (pragzip / pigz usage)
      label <- ""
      if (tolower(tool) == "rabbittrim") {
        pragzip_zero <- any(grepl("pragzip[^0-9]*\\b0\\b|rapidgzip[^0-9]*\\b0\\b", lines, ignore.case = TRUE))
        pigz_zero <- any(grepl("pigz[^0-9]*\\b0\\b", lines, ignore.case = TRUE))
        
        if (pragzip_zero && pigz_zero) {
          label <- "**"   # Neither pragzip nor pigz used
        } else if (pragzip_zero && !pigz_zero) {
          label <- "*"    # Only pragzip not used
        }
      }
      
      results <- rbind(results, data.frame(
        Tool = tool, Threads = threads, 
        Time_Min = t_min, Memory_GB = mem_gb, Label = label,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}

# --- Main Logic ---
cat("Parsing logs from:", LOG_DIR, "\n")
df_perf <- parse_scaling_logs(LOG_DIR)

# Standardize tool capitalization safely
standardize_tool <- function(t) {
  t_low <- tolower(t)
  case_when(
    t_low == "trimmomatic" ~ "Trimmomatic",
    t_low == "rabbittrim" ~ "RabbitTrim",
    TRUE ~ t
  )
}

if(nrow(df_perf) > 0) {
    df_perf$Tool <- standardize_tool(df_perf$Tool)
    df_perf$Tool <- factor(df_perf$Tool, levels = c("Trimmomatic", "RabbitTrim"))
} else {
    stop("No data parsed! Check log file path and format.")
}

# --- Plotting ---

tool_colors <- c("Trimmomatic" = "#E41A1C", "RabbitTrim" = "#377EB8")

# Shared Legend Title containing the Explanation
legend_title <- "Tool  (RabbitTrim marks: [**] without pragzip & pigz  |  [*] without pragzip)"

my_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

x_breaks <- c(8, 16, 32, 64, 128)

# Plot A: Wall Time
p1 <- ggplot(df_perf, aes(x = Threads, y = Time_Min, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 5, color = "black", show.legend = FALSE) +
  scale_color_manual(name = legend_title, values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  labs(title = "A. Wall Clock Time (SRR7890824 Dataset)", y = "Time (min)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Plot B: RAM Usage
p2 <- ggplot(df_perf, aes(x = Threads, y = Memory_GB, color = Tool, group = Tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 5, color = "black", show.legend = FALSE) +
  scale_color_manual(name = legend_title, values = tool_colors) +
  scale_x_continuous(trans = "log2", breaks = x_breaks) +
  labs(title = "B. Peak Memory Usage (SRR7890824 Dataset)", y = "RAM (GB)", x = "Threads") +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Combine Side by Side
layout <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save
ggsave(paste0(OUTPUT_BASE, ".pdf"), layout, width = 12, height = 6)
ggsave(paste0(OUTPUT_BASE, ".png"), layout, width = 12, height = 6, dpi = 300)

cat("\nPlots saved to:", paste0(OUTPUT_BASE, ".png"), "\n")