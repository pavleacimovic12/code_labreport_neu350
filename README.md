# code_labreport_neu350

# Load necessary libraries
library(tidyverse)
library(EMD)  # Required for Empirical Mode Decomposition

# Define file path
file_path <- "~/Downloads/PrimeNumbers11trials.txt"

# Read raw data as lines
raw_data <- readLines(file_path, encoding = "latin1", warn = FALSE)

# Identify the first row that contains numeric values (skip metadata)
data_start <- which(grepl("^[0-9]+\\s+[-+]?[0-9]*\\.?[0-9]+$", raw_data))[1]

# Extract only valid data lines
data_lines <- raw_data[data_start:length(raw_data)]

# Remove non-printable characters
data_lines <- iconv(data_lines, from = "latin1", to = "ASCII", sub = "")

# Keep only lines that contain exactly two numeric values
data_lines <- data_lines[grepl("^[0-9]+\\s+[-+]?[0-9]*\\.?[0-9]+$", data_lines)]

# Convert cleaned data into a dataframe
eeg_data <- read.table(text = paste(data_lines, collapse="\n"),
                       header = FALSE, col.names = c("Time", "Amplitude"),
                       stringsAsFactors = FALSE)

# Ensure numeric values and remove NA rows
eeg_data <- eeg_data %>% mutate(across(everything(), as.numeric)) %>% na.omit()

# Extract time and amplitude
time <- eeg_data$Time
signal <- eeg_data$Amplitude

emd_result <- emd(signal, max.imf = 9, tol = 0.05)  


# Extract IMFs
imfs <- emd_result$imf
num_imfs <- ncol(imfs)

par(mar = c(4, 4, 2, 1))  # (bottom, left, top, right)
plot(time, signal, type = "l", col = "black", main = "Original EEG Signal",
     xlab = "Time (s)", ylab = "Amplitude")


# Plot each IMF
for (i in 1:num_imfs) {
  plot(time, imfs[, i], type = "l", col = rainbow(num_imfs)[i],
       main = paste("IMF", i), xlab = "Time (s)", ylab = "Amplitude")
}

# Reset plotting layout
par(mfrow = c(1, 1))

# Save IMFs to CSV for further analysis
write.csv(imfs, "~/Downloads/eeg_imfs.csv", row.names = FALSE)

print(paste("Extracted", num_imfs, "IMFs from EEG signal."))




library(pracma)  # Needed for Hilbert transform

# Apply Hilbert Transform to each IMF
instantaneous_freqs <- apply(imfs, 2, function(imf) {
  phase <- Arg(hilbert(imf))  # Compute phase
  freq <- diff(phase) / (2 * pi * dt)  # Compute instantaneous frequency
  return(c(NA, freq))  # Keep same length
})

# Convert to DataFrame for analysis
instantaneous_freqs_df <- as.data.frame(instantaneous_freqs)
colnames(instantaneous_freqs_df) <- paste("IMF", 1:ncol(imfs))

# Save frequency data
write.csv(instantaneous_freqs_df, "~/Downloads/instantaneous_frequencies.csv", row.names = FALSE)

print("Instantaneous frequencies computed and saved.")


# Assign IMFs to EEG bands
  # Delta (slowest)
library(signal)
delta_band <- butter(2, 4/100, type="low") %>% filtfilt(delta_band)
theta_band <- imfs[, 4]  # Theta
alpha_band <- filter(alpha_band, rep(1/3, 3), sides=2)
beta_band  <- imfs[, 2]  # Beta
gamma_band <- imfs[, 1]  # Gamma (fastest)

# Create a dataframe with EEG bands
eeg_bands <- data.frame(
  Time = time,
  Delta = delta_band,
  Theta = theta_band,
  Alpha = alpha_band,
  Beta = beta_band,
  Gamma = gamma_band
)

# Save reconstructed bands
write.csv(eeg_bands, "~/Downloads/eeg_reconstructed_bands.csv", row.names = FALSE)





library(ggplot2)


eeg_long <- eeg_bands %>%
  pivot_longer(cols = -Time, names_to = "Band", values_to = "Amplitude")

ggplot(eeg_long, aes(x = Time, y = Amplitude, color = Band)) +
  geom_line(size = 1) +
  labs(title = "Reconstructed EEG Frequency Bands",
       x = "Time (s)", y = "Amplitude") +
  scale_color_manual(values = c("Delta" = "blue", "Theta" = "purple", 
                                "Alpha" = "green", "Beta" = "red", "Gamma" = "orange")) +
  theme_minimal()


