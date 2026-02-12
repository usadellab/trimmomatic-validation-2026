#!/usr/bin/env Rscript

# ---
# R Script: Trimming Tools Benchmark (Performance & Accuracy)
# Publication-Ready Figure Generation (Grouped: Plant vs Human)
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
LOG_DIR <- "/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs"
OUTPUT_BASE <- file.path(LOG_DIR, "Figure_Benchmark_Grouped")

# --- Parsing Functions ---

parse_performance_logs <- function(dir, suffix, dataset_label) {
  # Define expected filenames pattern
  tools <- c("Trimmomatic", "fastp", "BBDuk", "Skewer", "Cutadapt", "RabbitTrim")
  results <- data.frame(Tool = character(), Time_Min = numeric(), Memory_GB = numeric(), Dataset = character(), stringsAsFactors = FALSE)
  
  for (t in tools) {
    # Handle RabbitTrim naming variations if necessary
    fname_base <- tolower(t)
    if(t == "RabbitTrim" && !file.exists(file.path(dir, paste0("rabbittrim", suffix)))) {
       if(file.exists(file.path(dir, paste0("rabbit", suffix)))) fname_base <- "rabbit"
    }
    
    fpath <- file.path(dir, paste0(fname_base, suffix))
    
    if (file.exists(fpath)) {
      lines <- readLines(fpath)
      
      # Time Parsing
      time_line <- grep("Elapsed (wall clock) time", lines, value = TRUE, fixed = TRUE)
      if (length(time_line) > 0) {
        time_str <- str_extract(tail(time_line, 1), "[0-9:.]+$")
        parts <- as.numeric(unlist(str_split(time_str, ":")))
        t_min <- if(length(parts)==3) parts[1]*60 + parts[2] + parts[3]/60 else parts[1] + parts[2]/60
      } else { t_min <- NA }
      
      # Memory Parsing
      mem_line <- grep("Maximum resident set size", lines, value = TRUE, fixed = TRUE)
      if (length(mem_line) > 0) {
        mem_kb <- as.numeric(str_extract(tail(mem_line, 1), "\\d+"))
        mem_gb <- mem_kb / 1024 / 1024
      } else { mem_gb <- NA }
      
      results <- rbind(results, data.frame(Tool = t, Time_Min = t_min, Memory_GB = mem_gb, Dataset = dataset_label))
    }
  }
  return(results)
}

parse_accuracy_log <- function(filepath, dataset_label) {
  if (!file.exists(filepath)) return(NULL)
  lines <- readLines(filepath)
  data_lines <- lines[grep("\\|", lines)]
  data_lines <- data_lines[!grepl("matches|---|Matches|Tool", data_lines, ignore.case = TRUE)]
  
  acc_data <- data.frame(Tool = character(), Residuals = numeric(), Dataset = character(), stringsAsFactors = FALSE)
  
  for (line in data_lines) {
    parts <- str_trim(unlist(str_split(line, "\\|")))
    if (length(parts) >= 3) {
      tool <- parts[1]
      r1 <- as.numeric(gsub("[^0-9]", "", parts[2]))
      r2 <- as.numeric(gsub("[^0-9]", "", parts[3]))
      if (!is.na(r1) && !is.na(r2)) {
        acc_data <- rbind(acc_data, data.frame(Tool = tool, Residuals = r1 + r2, Dataset = dataset_label))
      }
    }
  }
  return(acc_data)
}

# --- Main Logic ---

cat("Parsing Manihot logs...\n")
perf_plant <- parse_performance_logs(LOG_DIR, ".time.log", "Plant")
acc_plant  <- parse_accuracy_log(file.path(LOG_DIR, "adapter_assessment.log"), "Plant")

cat("Parsing Human logs...\n")
perf_human <- parse_performance_logs(LOG_DIR, ".human.time.log", "Human")
acc_human  <- parse_accuracy_log(file.path(LOG_DIR, "human_adapter_assessment.log"), "Human")

# Merge Performance Data
perf_total <- rbind(perf_plant, perf_human)

# Merge Accuracy Data
if (!is.null(acc_plant) && !is.null(acc_human)) {
  acc_total <- rbind(acc_plant, acc_human)
  final_df <- merge(perf_total, acc_total, by = c("Tool", "Dataset"), all = TRUE)
} else {
  final_df <- perf_total
  final_df$Residuals <- NA
}

# Order Tools
# Here we order alphabetically or explicitly
tool_order <- c("Trimmomatic", "RabbitTrim", "fastp", "BBDuk", "Skewer", "Cutadapt")
final_df$Tool <- factor(final_df$Tool, levels = rev(tool_order)) # rev for coord_flip

# Order Dataset (Plant first, then Human)
final_df$Dataset <- factor(final_df$Dataset, levels = c("Plant", "Human"))

write_csv(final_df, paste0(OUTPUT_BASE, ".csv"))

# --- Plotting ---

# Colors
dataset_colors <- c("Plant" = "#009E73", "Human" = "#0072B2")

# Theme
pub_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.text.y = element_text(color = "black", size = 11),
    axis.text.x = element_text(color = "black", size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot A: Speed
p1 <- ggplot(final_df, aes(x = Tool, y = Time_Min, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Time_Min)), 
            position = position_dodge(width = 0.8), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = dataset_colors, breaks = c("Human", "Plant")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "A", y = "Wall time (min)", x = NULL) +
  pub_theme +
  theme(legend.position = "none")

# Plot B: Memory
p2 <- ggplot(final_df, aes(x = Tool, y = Memory_GB, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Memory_GB)), 
            position = position_dodge(width = 0.8), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = dataset_colors, breaks = c("Human", "Plant")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "B", y = "Peak RAM (GB)", x = NULL) +
  pub_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")

# Plot C: Accuracy
p3 <- ggplot(final_df, aes(x = Tool, y = Residuals, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = comma(Residuals)), 
            position = position_dodge(width = 0.8), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_fill_manual(values = dataset_colors, breaks = c("Human", "Plant")) + 
  scale_y_log10(labels = label_comma(), expand = expansion(mult = c(0, 0.3))) + 
  labs(title = "C", y = "Residual Adapter Reads (Log Scale)", x = NULL) +
  pub_theme +
  theme(legend.position = "none")

# Combine and Save
layout <- (p1 + p2) / p3 + plot_layout(heights = c(1, 1.2))
ggsave(paste0(OUTPUT_BASE, ".pdf"), layout, width = 10, height = 8)
ggsave(paste0(OUTPUT_BASE, ".png"), layout, width = 10, height = 8, dpi = 300)

cat("\nPublication figures saved to:", paste0(OUTPUT_BASE, ".png"), "\n")