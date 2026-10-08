test_that("differential-expression direction uses red-up and blue-down", {
  colors <- .de_direction_colors()
  expect_identical(unname(colors[["up"]]), "red2")
  expect_identical(unname(colors[["down"]]), "royalblue")

  labels <- .de_direction_label(
    c(2, -2, 0, NA_real_),
    c(TRUE, TRUE, TRUE, TRUE)
  )
  expect_identical(
    as.character(labels),
    c("Upregulated", "Downregulated", "Not significant", "Not significant")
  )
  expect_identical(
    unname(.de_direction_scale()),
    c("grey70", "royalblue", "red2")
  )
})

test_that("directional MA plot has shared directional scale", {
  skip_if_not_installed("ggplot2")
  res <- data.frame(
    baseMean = c(10, 20, 30),
    log2FoldChange = c(2, -2, 0.1),
    padj = c(0.001, 0.001, 0.8)
  )

  plot <- .plot_directional_ma(res, "test", 0.05, 1)
  expect_identical(
    ggplot2::ggplot_build(plot)$data[[1]]$colour,
    c("red2", "royalblue", "grey70")
  )
  expect_identical(
    as.character(plot$data$direction),
    c("Upregulated", "Downregulated", "Not significant")
  )
})
