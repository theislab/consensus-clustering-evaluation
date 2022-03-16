#!/usr/bin/env Rscript

"
Plot metrics

Usage:
    plot_metrics.R --out-file=<path> [options] <file>

Options:
    -h --help             Show this screen.
    --out-file=<path>     Path to output file.
" -> doc

plot_metrics <- function(metrics) {

    message("Making metrics plots...")
    base_plot <- ggplot2::ggplot(
        metrics,
        ggplot2::aes(
            y      = .data$Score,
            colour = .data$Method,
            shape  = .data$Dataset
        )
    ) +
    ggplot2::geom_jitter(size = 2, alpha = 0.8, width = 0.2) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = 90)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = NA),
            strip.background = ggplot2::element_rect(fill = "black"),
            strip.text       = ggplot2::element_text(colour = "white")
        )

    list(
        base_plot +
            ggplot2::aes(x = .data$Method) +
            ggplot2::facet_grid(.data$Metric ~ .),
        base_plot +
            ggplot2::aes(x = .data$Dataset) +
            ggplot2::facet_grid(.data$Metric ~ .)
    )
}

plot_k <- function(metrics) {

    metrics <- dplyr::filter(metrics, .data$Metric == "k")

    ref_lines <- dplyr::filter(metrics, .data$Method == "random")

    plot <- ggplot2::ggplot(
        metrics,
        ggplot2::aes(
            x    = .data$Method,
            y    = .data$Score,
            fill = .data$Method
        )
    ) +
        ggplot2::geom_col() +
        ggplot2::geom_hline(
            data = ref_lines,
            ggplot2::aes(yintercept = .data$Score),
            colour = "red"
        ) +
        ggplot2::guides(x = ggplot2::guide_axis(angle = 90)) +
        ggplot2::labs(y = "Number of clusters") +
        ggplot2::facet_grid(.data$Dataset ~ ., scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position  = "none",
            panel.background = ggplot2::element_rect(fill = NA),
            strip.background = ggplot2::element_rect(fill = "black"),
            strip.text       = ggplot2::element_text(colour = "white")
        )

    list(plot)

}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file     <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    message("Reading metrics from '", file, "'...")
    all_metrics <- readr::read_tsv(file)
    metrics <- dplyr::filter(all_metrics, !(Metric %in% c("k", "k_prop")))
    metrics_plots <- plot_metrics(metrics)
    k_metrics <- dplyr::filter(all_metrics, Metric %in% c("k", "k_prop"))
    k_plots <- plot_k(k_metrics)
    plots <- c(metrics_plots, k_plots)
    message("Writing plots to '", out_file, "'...")
    pdf(out_file, width = 8.3, height = 11.7)
    for (plot in plots) {
        print(plot)
    }
    dev.off()
    message("Done!")
}
