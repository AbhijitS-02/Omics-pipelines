offset_values <- c(10, 50, 100)

par(mfrow = c(1, length(offset_values)))  # Layout for plots

for (offset in offset_values) {
  corrected <- backgroundCorrect(raw.GPL6480, method = "normexp", offset = offset)
  hist(log2(corrected$E[,1]), main = paste("Offset =", offset), xlab = "log2(Intensity)", col = "skyblue")
}


offset_values <- c(10, 50, 100)

par(mfrow = c(1, length(offset_values)))  # Layout for plots

for (offset in offset_values) {
  corrected <- backgroundCorrect(raw.GPL14550, method = "normexp", offset = offset)
  hist(log2(corrected$E[,1]), main = paste("Offset =", offset), xlab = "log2(Intensity)", col = "skyblue")
}