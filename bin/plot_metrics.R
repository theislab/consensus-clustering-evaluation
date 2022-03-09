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
    ggplot2::ggplot(
        metrics,
        ggplot2::aes(
            x = .data$Method, y = .data$Score,
            colour = .data$Method,
            shape  = .data$Dataset
        )
    ) +
        ggplot2::geom_jitter(size = 2, alpha = 0.8, width = 0.2) +
        ggplot2::facet_grid(.data$Metric ~ .) +
        ggplot2::guides(x = ggplot2::guide_axis(angle = 90)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = NA),
            strip.background = ggplot2::element_rect(fill = "black"),
            strip.text       = ggplot2::element_text(colour = "white")
        )

}

if (sys.nframe() == 0) {
    args <- docopt::docopt(doc)

    file     <- args[["<file>"]]
    out_file <- args[["--out-file"]]

    message("Reading metrics from '", file, "'...")
    metrics <- readr::read_tsv(file)
    plot <- plot_metrics(metrics)
    message("Writing plots to '", out_file, "'...")
    ggplot2::ggsave(out_file, plot = plot, width = 21, height = 29.7,
                    units = "cm")
    message("Done!")
}
