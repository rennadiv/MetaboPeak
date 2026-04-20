#' @title Barplots of 3 lowest variable peaks based on mass (m/z).
#' @description
#' The function calculates CV coefficient (SD/Mean) of each peak with the given mass (m/z)
#' and plots abundances of all samples for the max 3 lowest CV coefficient.
#'
#' @usage barPeaks(x, mass, n)
#' @param x data frame
#' @param mass peak mass (m/z)
#' @param n number of samples
#'
#' @details
#' The format of the data frame should be: Peak names|abundances of samples|m.z|RT|...
#'
#' The mass should be put as a string.
#'
#' @returns Bar plots of max 3 peaks with similar m/z and their RT and CV coefficient.
#'
#' @examples
#' barPeaks(neg, '355.10', 48)
#'
#' @importFrom graphics barplot
#' @importFrom graphics par
#' @importFrom stats complete.cases
#' @importFrom stats sd
#'
#' @export


barPeaks <- function(x, mass, n){

  # ---- checks ----
  if (!is.data.frame(x)) {
    stop("Input x must be a data.frame")
  }

  if (!is.numeric(n) || length(n) != 1) {
    stop("n must be a single number (number of samples)")
  }

  # ---- data prep ----
  number.of.samples <- n
  data <- x
  n3 <- number.of.samples + 3

  a <- data[substr(data$m.z, 1, nchar(mass)) == mass, c(2:n3)]

  if (nrow(a) == 0) {
    stop("This mass is not in the data frame.")
  }

  aa <- a[, -c(number.of.samples+1, number.of.samples+2)]

  Mean <- rowMeans(aa, na.rm = TRUE)
  SD <- apply(aa, 1, sd, na.rm = TRUE)
  CV <- SD / Mean

  # ---- select top peaks ----
  ord <- order(CV)
  top_idx <- ord[1:min(3, length(ord))]

  a_top <- a[top_idx, ]
  aa_top <- aa[top_idx, ]
  CV_top <- CV[top_idx]

  peak_names <- toy_df[ord,1]

  # ---- plotting ----
  par(mfrow = c(nrow(a_top), 1))

  mids_list <- list()

  for (i in seq_len(nrow(a_top))) {

    mids <- barplot(
      as.matrix(aa_top[i, ]),
      main = paste(
        "m/z =", round(a_top[i, number.of.samples+1], 3),
        ", RT =", round(a_top[i, number.of.samples+2], 2),
        ", CV =", round(CV_top[i], 3)
      ),
      names.arg = rep("", n),
      col = "navy"
    )

    mids_list[[i]] <- mids
  }

  return(list(
    peak_names = peak_names,
    mids = mids_list,
    CV = CV_top,
    data = aa_top
  ))
}
