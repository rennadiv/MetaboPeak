barPeaks_ggplot <- function(x, mass, n){

  library(ggplot2)
  library(tidyr)

  number.of.samples <- n
  n3 <- number.of.samples + 3

  a <- x[substr(x$m.z, 1, nchar(mass)) == mass, c(2:n3)]

  if (nrow(a) == 0) {
    stop("This mass is not in the data frame.")
  }

  aa <- a[, -c(number.of.samples+1, number.of.samples+2)]

  Mean <- rowMeans(aa, na.rm = TRUE)
  SD <- apply(aa, 1, sd, na.rm = TRUE)
  CV <- SD / Mean

  ord <- order(CV)
  top_idx <- ord[1:min(3, length(ord))]

  a_top <- a[top_idx, ]
  aa_top <- aa[top_idx, ]
  CV_top <- CV[top_idx]

  plots <- list()

  for (i in seq_len(nrow(aa_top))) {

    df_plot <- data.frame(
      sample = colnames(aa_top),
      intensity = as.numeric(aa_top[i, ])
    )

    p <- ggplot(df_plot, aes(x = sample, y = intensity)) +
      geom_col(fill = "navy") +
      labs(
        title = paste(
          "m/z =", round(a_top[i, number.of.samples+1], 3),
          ", RT =", round(a_top[i, number.of.samples+2], 2),
          ", CV =", round(CV_top[i], 3)
        )
      )

    plots[[i]] <- p
  }

  return(plots)
}
