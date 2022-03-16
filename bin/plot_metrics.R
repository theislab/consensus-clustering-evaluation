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

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file     <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    message("Reading metrics from '", file, "'...")
    metrics <- readr::read_tsv(file)
    plots <- plot_metrics(metrics)
    message("Writing plots to '", out_file, "'...")
    pdf(out_file, width = 8.3, height = 11.7)
    for (plot in plots) {
        print(plot)
    }
    dev.off()
    message("Done!")
}
