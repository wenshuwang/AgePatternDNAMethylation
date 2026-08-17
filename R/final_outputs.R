# Final tables and filtering figures ------------------------------------------

require_final_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for the final figures. Install it ",
      "and source run_final_outputs.R again.",
      call. = FALSE
    )
  }
}

load_final_aacpg_objects <- function(project_root = ".") {
  paths <- c(
    all_list = "relevant_rds/all_list.rds",
    male_list = "relevant_rds/male_list.rds",
    female_list = "relevant_rds/female_list.rds",
    all_results = "relevant_rds/all_results.rds",
    male_results = "relevant_rds/male_results.rds",
    female_results = "relevant_rds/female_results.rds"
  )
  full_paths <- file.path(project_root, paths)
  missing <- paths[!file.exists(full_paths)]
  if (length(missing)) {
    stop(
      "Final-output inputs are missing:\n- ",
      paste(missing, collapse = "\n- "),
      "\nComplete filtering before generating the final outputs.",
      call. = FALSE
    )
  }

  objects <- lapply(full_paths, readRDS)
  names(objects) <- names(paths)
  objects
}

build_aacpg_table <- function(all_list, male_list, female_list) {
  windows <- names(all_list)
  if (!identical(names(male_list), windows) ||
      !identical(names(female_list), windows)) {
    stop("All, male, and female lists do not contain the same age windows.",
         call. = FALSE)
  }

  same_window_overlaps <- Map(
    intersect,
    male_list[windows],
    female_list[windows]
  )

  table_data <- data.frame(
    `Age Window` = paste0("[", sub("-", ", ", windows), "]"),
    All = as.integer(lengths(all_list[windows])),
    Male = as.integer(lengths(male_list[windows])),
    Female = as.integer(lengths(female_list[windows])),
    `Overlap (M/F)` = as.integer(lengths(same_window_overlaps)),
    check.names = FALSE
  )

  # The manuscript's total overlap is the number of unique CpGs appearing in
  # a male/female overlap within the same age window. It is not the broader
  # intersection of the two cohort-wide unions, where windows may differ.
  total_row <- data.frame(
    `Age Window` = "Total",
    All = length(unique(unlist(all_list, use.names = FALSE))),
    Male = length(unique(unlist(male_list, use.names = FALSE))),
    Female = length(unique(unlist(female_list, use.names = FALSE))),
    `Overlap (M/F)` = length(unique(unlist(
      same_window_overlaps,
      use.names = FALSE
    ))),
    check.names = FALSE
  )

  rbind(table_data, total_row)
}

draw_aacpg_table <- function(table_data, include_caption = TRUE) {
  require_final_package("grid")

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = 0.5, y = 0.5,
    width = grid::unit(0.94, "npc"),
    height = grid::unit(0.92, "npc")
  ))

  if (include_caption) {
    grid::grid.text(
      paste0(
        "Table 1. Number of aaCpGs selected for the whole cohort, male-only,\n",
        "and female-only groups in different age windows."
      ),
      x = 0.03, y = 0.98, just = c("left", "top"),
      gp = grid::gpar(
        fontfamily = "sans",
        fontface = "bold",
        fontsize = 14,
        lineheight = 1.15
      )
    )
  }

  top_y <- if (include_caption) 0.84 else 0.94
  header_y <- top_y - 0.055
  subheader_y <- top_y - 0.11
  rule_y <- top_y - 0.145
  row_step <- 0.044
  row_start <- rule_y - 0.042

  column_x <- c(0.10, 0.37, 0.52, 0.67, 0.88)
  column_just <- c("left", "right", "right", "right", "right")

  grid::grid.lines(
    x = c(0.03, 0.97), y = c(top_y, top_y),
    gp = grid::gpar(lwd = 1.1)
  )
  grid::grid.text(
    "# aaCpGs", x = mean(column_x[2:5]), y = header_y,
    gp = grid::gpar(fontfamily = "sans", fontsize = 15)
  )
  grid::grid.text(
    "Age Window", x = column_x[1], y = subheader_y,
    just = "left", gp = grid::gpar(fontfamily = "sans", fontsize = 14)
  )

  headers <- colnames(table_data)[2:5]
  for (j in seq_along(headers)) {
    grid::grid.text(
      headers[j], x = column_x[j + 1], y = subheader_y,
      just = "right", gp = grid::gpar(fontfamily = "sans", fontsize = 14)
    )
  }

  grid::grid.lines(
    x = c(0.03, 0.97), y = c(rule_y, rule_y),
    gp = grid::gpar(lwd = 1.6)
  )
  grid::grid.lines(
    x = c(0.03, 0.97), y = c(rule_y - 0.008, rule_y - 0.008),
    gp = grid::gpar(lwd = 0.65)
  )

  formatted <- table_data
  for (j in 2:5) {
    formatted[[j]] <- format(
      table_data[[j]],
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  }

  for (i in seq_len(nrow(formatted))) {
    y <- row_start - (i - 1) * row_step
    is_total <- i == nrow(formatted)
    row_gp <- grid::gpar(
      fontfamily = "sans",
      fontsize = 13.5,
      fontface = if (is_total) "bold" else "plain"
    )
    for (j in seq_len(ncol(formatted))) {
      grid::grid.text(
        as.character(formatted[i, j]),
        x = column_x[j], y = y,
        just = column_just[j], gp = row_gp
      )
    }
  }

  bottom_y <- row_start - (nrow(formatted) - 0.35) * row_step
  grid::grid.lines(
    x = c(0.03, 0.97), y = c(bottom_y, bottom_y),
    gp = grid::gpar(lwd = 1.1)
  )
  grid::popViewport()
  invisible(table_data)
}

save_aacpg_table <- function(table_data, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  csv_path <- file.path(output_dir, "table1_aacpg_counts.csv")
  pdf_path <- file.path(output_dir, "table1_aacpg_counts.pdf")
  png_path <- file.path(output_dir, "table1_aacpg_counts.png")

  utils::write.csv(table_data, csv_path, row.names = FALSE)

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  pdf_device(pdf_path, width = 11, height = 8.5)
  draw_aacpg_table(table_data)
  grDevices::dev.off()

  grDevices::png(png_path, width = 11, height = 8.5, units = "in", res = 300)
  draw_aacpg_table(table_data)
  grDevices::dev.off()

  invisible(c(csv = csv_path, pdf = pdf_path, png = png_path))
}

save_figure3 <- function(objects, clock_cpgs, output_dir) {
  packages <- c("ComplexHeatmap", "cowplot", "ggpubr", "ggplot2", "grid")
  invisible(lapply(packages, require_final_package))

  all_list <- objects$all_list
  male_list <- objects$male_list
  female_list <- objects$female_list
  all_results <- objects$all_results
  male_results <- objects$male_results
  female_results <- objects$female_results

  all_cpg <- unique(unlist(all_list, use.names = FALSE))
  male_cpg <- unique(unlist(male_list, use.names = FALSE))
  female_cpg <- unique(unlist(female_list, use.names = FALSE))

  window_upset <- ComplexHeatmap::make_comb_mat(all_list, mode = "intersect")
  panel_a <- grid::grid.grabExpr({
    grid::grid.newpage()
    ComplexHeatmap::draw(ComplexHeatmap::UpSet(
      window_upset,
      set_order = rownames(window_upset),
      pt_size = grid::unit(1.5, "mm"),
      lwd = 1
    ))
  })

  cohort_sets <- list(
    "All aaCpGs" = all_cpg,
    "Male aaCpGs" = male_cpg,
    "Female aaCpGs" = female_cpg,
    "Clock CpGs" = clock_cpgs
  )
  cohort_upset <- ComplexHeatmap::make_comb_mat(cohort_sets, mode = "intersect")
  panel_b <- grid::grid.grabExpr({
    grid::grid.newpage()
    ComplexHeatmap::draw(ComplexHeatmap::UpSet(
      cohort_upset,
      set_order = rownames(cohort_upset)
    ))
  })

  upset_row <- cowplot::plot_grid(
    panel_a, panel_b,
    ncol = 2,
    rel_widths = c(1.2, 0.8),
    labels = c("A", "B")
  )

  window_columns <- grep("-corr_res$", colnames(all_results), value = TRUE)
  scatter_plots <- lapply(window_columns, function(column) {
    window <- sub("-corr_res$", "", column)
    male_window <- male_list[[window]]
    female_window <- female_list[[window]]
    plotted_cpgs <- union(male_window, female_window)

    plot_data <- data.frame(
      male = male_results[plotted_cpgs, column],
      female = female_results[plotted_cpgs, column],
      row.names = plotted_cpgs
    )
    overlap <- intersect(male_window, female_window)
    plot_data$Group <- ifelse(
      rownames(plot_data) %in% overlap,
      "Overlap",
      ifelse(
        rownames(plot_data) %in% male_window,
        "Male Unique",
        "Female Unique"
      )
    )
    plot_data$Group <- factor(
      plot_data$Group,
      levels = c("Male Unique", "Overlap", "Female Unique")
    )

    ggplot2::ggplot(plot_data, ggplot2::aes(
      x = male, y = female, color = Group
    )) +
      ggplot2::geom_point(size = 0.25, alpha = 0.5, na.rm = TRUE) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
      ggplot2::scale_color_manual(values = c(
        "Male Unique" = "#ADD8E6",
        "Overlap" = "#777777",
        "Female Unique" = "#FFB6C1"
      )) +
      ggplot2::coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
      ggplot2::labs(title = paste0("[", window, "]"), x = "Male", y = "Female") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.border = ggplot2::element_rect(fill = NA, color = "black")
      )
  })

  scatter_grid <- ggpubr::ggarrange(
    plotlist = scatter_plots,
    ncol = 5,
    nrow = 3,
    common.legend = TRUE,
    legend = "top"
  )
  scatter_grid <- cowplot::plot_grid(scatter_grid, labels = "C")
  figure3 <- cowplot::plot_grid(
    upset_row,
    scatter_grid,
    ncol = 1,
    rel_heights = c(1.9, 4)
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(output_dir, "figure3.pdf"),
    plot = figure3,
    device = grDevices::cairo_pdf,
    width = 14,
    height = 12,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "figure3.png"),
    plot = figure3,
    width = 14,
    height = 12,
    dpi = 200,
    bg = "white"
  )

  male_upset <- ComplexHeatmap::make_comb_mat(male_list, mode = "intersect")
  male_grob <- grid::grid.grabExpr({
    ComplexHeatmap::draw(ComplexHeatmap::UpSet(
      male_upset,
      set_order = rownames(male_upset)
    ), newpage = FALSE)
  })
  female_upset <- ComplexHeatmap::make_comb_mat(female_list, mode = "intersect")
  female_grob <- grid::grid.grabExpr({
    ComplexHeatmap::draw(ComplexHeatmap::UpSet(
      female_upset,
      set_order = rownames(female_upset)
    ), newpage = FALSE)
  })
  supplementary_overlap <- cowplot::plot_grid(
    male_grob,
    female_grob,
    ncol = 1,
    labels = c("A", "B")
  )
  ggplot2::ggsave(
    file.path(output_dir, "supplementary_figure4.pdf"),
    plot = supplementary_overlap,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  invisible(figure3)
}

generate_quick_final_outputs <- function(project_root = ".",
                                         output_dir = "final_outputs") {
  objects <- load_final_aacpg_objects(project_root)
  table_data <- build_aacpg_table(
    objects$all_list,
    objects$male_list,
    objects$female_list
  )

  expected <- c(All = 19423L, Male = 32210L, Female = 11422L,
                `Overlap (M/F)` = 8105L)
  observed <- unlist(table_data[nrow(table_data), names(expected)], use.names = TRUE)
  observed <- as.integer(observed)
  names(observed) <- names(expected)
  if (!identical(observed, expected)) {
    stop(
      "The saved aaCpG lists do not match the validated final checkpoint.\n",
      "Observed: ", paste(names(observed), observed, collapse = ", "), "\n",
      "Expected: ", paste(names(expected), expected, collapse = ", "),
      call. = FALSE
    )
  }

  output_dir <- file.path(project_root, output_dir)
  table_paths <- save_aacpg_table(table_data, output_dir)

  clock_path <- file.path(project_root, "relevant_rds", "clock_cpgs.rds")
  if (!file.exists(clock_path)) {
    fallback <- file.path(project_root, "intermediates", "POST_clock_sites.rds")
    if (!file.exists(fallback)) {
      stop("No filtered clock CpG list was found.", call. = FALSE)
    }
    clock_cpgs <- readRDS(fallback)
  } else {
    clock_cpgs <- readRDS(clock_path)
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  save_figure3(objects, clock_cpgs, output_dir)

  message("Quick final outputs complete in: ", normalizePath(output_dir))
  invisible(list(table = table_data, paths = table_paths))
}
