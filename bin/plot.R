#!/usr/bin/Rscript --vanilla

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Cairo))

load.table <- function(infile, header=FALSE, col.names=c()) {
    result <- tryCatch(
        if (header) {
            read.table(infile, header=TRUE)
        } else {
            read.table(infile, col.names=col.names)
        },
        error=function(e) NULL,
        warning=function(w) {})
    if (is.null(result)) {
        cat(paste("Skipping missing file '", infile, "'.\n", sep=""))
    }
    return(result)
}

plot.pressure <- function(infile="pressure.ssv") {
    p.data <- load.table(infile, col.names=c("psi", "radius", "pressure"))
    if (is.null(p.data)) { return(NULL) }
    # p.data <- read.table(infile, col.names=c("psi", "radius", "pressure"))
    p.data$psi <- factor(p.data$psi, levels=rev(sort(unique(p.data$psi))),
                         ordered=TRUE)

    p <- ggplot(p.data, aes(x=radius, y=pressure, colour=psi, group=psi)) +
        geom_line() +
        scale_colour_hue(expression(Psi)) +
        scale_x_continuous(breaks=seq(0, 1.8, by=0.2)) +
        scale_y_continuous(breaks=seq(0, 30, by=5)) +
        coord_cartesian(xlim=c(0, 1.8), ylim=c(0, 30)) +
        xlab(expression(r[rho])) +
        ylab("P (kPa)")

    return(p)
}

plot.stress.strain <- function(infile="stress_strain.ssv") {
    p.data <- load.table(infile, col.names=c("psi", "eps", "str", "str_a"))
    if (is.null(p.data)) { return(NULL) }
    p.data$psi <- factor(p.data$psi, levels=rev(sort(unique(p.data$psi))),
                         ordered=TRUE)

    p <- ggplot(p.data, aes(x=eps, y=str, colour=psi, group=psi)) +
        geom_line() +
        geom_line(aes(y=str_a), linetype="dashed") +
        scale_colour_hue(expression(Psi)) +
        scale_x_continuous(breaks=seq(-1, 0.8, by=0.2)) +
        scale_y_continuous(breaks=seq(0, 150, by=25)) +
        coord_cartesian(xlim=c(-1, 0.8), ylim=c(0, 150)) +
        xlab(expression(epsilon)) +
        ylab(expression(paste(sigma, " (kPa)")))

    return(p)
}

plot.radius <- function(infile="radius.ssv") {
    p.data <- load.table(infile, col.names=c("p", "r", "psi"))
    if (is.null(p.data)) { return(NULL) }
    p.data$psi <- factor(p.data$psi, levels=rev(sort(unique(p.data$psi))),
                         ordered=TRUE)
    p <- ggplot(p.data, aes(x=r, y=p, colour=psi, group=psi)) +
        geom_path() +
        scale_colour_hue(expression(Psi[MYO](P))) +
        scale_x_continuous(breaks=seq(0, 1.8, by=0.2)) +
        scale_y_continuous(breaks=seq(0, 30, by=5)) +
        coord_cartesian(xlim=c(0, 1.8), ylim=c(0, 30)) +
        xlab(expression(r[rho])) +
        ylab("P (kPa)")

    return(p)
}

plot.aa.state <- function(infile="profiles.ssv") {
    p.data <- load.table(infile, col.names=c("x", "norm_rad", "abs_rad", "res",
                                     "pa", "psi", "psimyo", "psitgf", "shape",
                                     "p_art"))
    if (is.null(p.data)) { return(NULL) }
    p.data$shape <- factor(p.data$shape)
    p.data$p_art <- factor(round(p.data$p_art * 7.50061683))
    p.data$pa <- p.data$pa * 7.50061683

    p1 <- ggplot(p.data, aes(x=x, y=abs_rad, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("r (um)") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p2 <- ggplot(p.data, aes(x=x, y=norm_rad, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("r (norm)") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p3 <- ggplot(p.data, aes(x=x, y=res, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("res") +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p4 <- ggplot(p.data, aes(x=x, y=pa, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("P") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p5 <- ggplot(p.data, aes(x=x, y=psi, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi)) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p6 <- ggplot(p.data, aes(x=x, y=psimyo, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi[MYO])) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p7 <- ggplot(p.data, aes(x=x, y=psitgf, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi[TGF])) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p8 <- ggplot(p.data, aes(x=pa, y=psimyo, colour=p_art, linetype=shape)) +
        geom_path() +
        #geom_point(data=p.data[seq(1,nrow(df), 10), ]) +
        geom_point() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("P (kPa)") +
        ylab(expression(Psi[MYO]))

    plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)

    return(plots)
}

plot.aa.state.glom <- function(infile="profiles_glom.ssv") {
    p.data <- load.table(infile, col.names=c("x", "norm_rad", "abs_rad", "res",
                                     "pa", "psi", "psimyo", "psitgf", "shape",
                                     "p_art", "sngfr"))
    if (is.null(p.data)) { return(NULL) }
    p.data$shape <- factor(p.data$shape)
    p.data$p_eff <- round(p.data$p_art * 7.50061683 * 1.5)
    p.data$p_art <- factor(round(p.data$p_art * 7.50061683))
    p.data$pa <- p.data$pa * 7.50061683

    p1 <- ggplot(p.data, aes(x=x, y=abs_rad, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("r (um)") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p2 <- ggplot(p.data, aes(x=x, y=norm_rad, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("r (norm)") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p3 <- ggplot(p.data, aes(x=x, y=res, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("res") +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p3a <- ggplot(p.data, aes(x=p_eff, y=sngfr)) +
        geom_line() +
        geom_point() +
        xlab(expression(P[ART] ~ " (mmHg)")) +
        ylab("SNGFR")

    p4 <- ggplot(p.data, aes(x=x, y=pa, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab("P") +
        expand_limits(y=0) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p5 <- ggplot(p.data, aes(x=x, y=psi, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi)) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p6 <- ggplot(p.data, aes(x=x, y=psimyo, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi[MYO])) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p7 <- ggplot(p.data, aes(x=x, y=psitgf, colour=p_art, linetype=shape)) +
        geom_line() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("x") +
        ylab(expression(Psi[TGF])) +
        scale_x_continuous(breaks=seq(from=0, to=1, by=0.2))

    p8 <- ggplot(p.data, aes(x=pa, y=psimyo, colour=p_art, linetype=shape)) +
        geom_path() +
        #geom_point(data=p.data[seq(1,nrow(df), 10), ]) +
        geom_point() +
        scale_colour_hue(expression(P[ART])) +
        scale_linetype_discrete("shape") +
        xlab("P (kPa)") +
        ylab(expression(Psi[MYO]))

    plots <- list(p1, p2, p3, p3a, p4, p5, p6, p7, p8)

    return(plots)
}

plot.aa.bf <- function(infile="aa_bf.ssv") {
    p.data <- load.table(infile, header=TRUE)
    if (is.null(p.data)) { return(NULL) }
    p.data$pa <- p.data$pa * 7.5 # convert kPa to mmHg
    p.data$Q <- p.data$Q * 60 # convert nl/s to nl/min
    p.data$n <- factor(p.data$n)
    p.data$slope <- factor(p.data$slope)
    p <- ggplot(p.data, aes(x=pa, y=Q, colour=n, linetype=slope)) +
        geom_line() +
        scale_colour_hue("Segments") +
        scale_linetype_discrete("slope") +
        #scale_x_continuous(breaks=seq(0, 25, by=5)) +
        ## scale_y_continuous(breaks=seq(0, 15, by=5)) +
        #coord_cartesian(xlim=c(0, 25), ylim=c(0, 15)) +
        xlab(expression(paste(P[a], " (mmHg)"))) +
        ylab("Q (nl/min)")

    return(p)
}

plot.sngfr.lookup <- function(infile="sngfr_lookup_base.ssv") {
    p.data <- load.table(infile, header=TRUE)
    if (is.null(p.data)) { return(NULL) }
    p.data$ptgf <- factor(p.data$ptgf)
    p.data$rpp <- factor(p.data$rpp)
    # kPa . s / um^3 -> mmHg . min / nL
    p.data$res <- p.data$res * 1.2501028e5

    s.data <- p.data[p.data$slope == 2,]

    p1 <- ggplot(s.data, aes(x=rpp, y=res, colour=ptgf, group=ptgf)) +
        geom_line() +
        scale_colour_hue(expression(Psi[TGF])) +
        xlab("RPP (mmHg)") +
        ## ylab("Resistance (kPa . s / um^3)") +
        ylab("Resistance (mmHg . min / nL)") +
        scale_x_discrete(breaks=seq(40, 180, by=20)) +
        guides(col = guide_legend(ncol=3))

    f.data <- data.frame(
        rpp  = as.numeric(s.data$rpp[s.data$ptgf == "0.01"]) + 39,
        rmin = s.data$res[s.data$ptgf == "0.01"],
        rmax = s.data$res[s.data$ptgf == "0.11"])

    p2 <- ggplot(f.data, aes(x=rpp)) +
        geom_line(aes(y=rmin)) +
        geom_line(aes(y=rmax)) +
        geom_ribbon(aes(ymin=rmin, ymax=rmax)) +
        xlab("RPP (mmHg)") +
        ## ylab("Resistance (kPa . s / um^3)") +
        ylab("Resistance (mmHg . min / nL)") +
        scale_x_continuous(breaks=seq(40, 180, by=20)) +
        scale_y_continuous(breaks=seq(0, 0.3, by=0.05)) +
        coord_cartesian(xlim=c(40, 180), ylim=c(0, 0.3)) +
        theme(legend.position="none")

    plots <- list(p1, p2)

    return(plots)
}

plot.to.file <- function(filename, plots, width=5, height=5) {
    if (! is.null(plots)) {
        CairoPDF(file=filename, bg="transparent", width=width, height=height)
        if (is.null(names(plots))) {
            # Since the list has no named attributes, it isn't a plot.
            # Therefore, we assume that it is a list of plots.
            for (i in seq_along(plots)) {
                suppressWarnings(print(plots[[i]]))
            }
        } else {
            suppressWarnings(print(plots))
        }
        ignore.this <- dev.off()
    }
}

plot.figures <- function(set.theme=TRUE) {
    CairoFonts(regular="Charis SIL:style=Regular", bold="Charis SIL:style=Bold")

    if (set.theme) {
        base_size <- 16
        plot.theme <- theme_grey(base_size)
        plot.theme$panel.background <- element_rect(fill = NA, colour = NA)
        plot.theme$panel.grid.major <- element_line(colour = "grey25",
                                                    size = 0.50)
        plot.theme$panel.grid.minor <- element_line(colour = "grey50",
                                                    size = 0.25)
        plot.theme$axis.text.x <- element_text(colour = "black",
                                               size = base_size * 0.8,
                                               vjust = 1, lineheight = 0.9)
        plot.theme$axis.text.y <- element_text(colour = "black",
                                               size = base_size * 0.8,
                                               hjust = 1, lineheight = 0.9)
        plot.theme$panel.border <- element_rect(fill = NA, colour = "black",
                                                size = 1.5)
        plot.theme$strip.background <- element_rect(colour = NA, fill = NA)
        plot.theme$strip.text.x <- element_text(colour = "black",
                                                size = base_size)
        plot.theme$legend.background <- element_rect(colour = "black",
                                                     fill="#fdf6e3")
        plot.theme$plot.background <- element_rect(fill = NA, colour = NA)
        theme_set(plot.theme)
    }

    plot.to.file("stress_strain.pdf", plot.stress.strain())
    plot.to.file("pressure.pdf", plot.pressure())
    plot.to.file("radius.pdf", plot.radius())
    plot.to.file("profiles.pdf", plot.aa.state())
    plot.to.file("profiles_glom.pdf", plot.aa.state.glom())
    plot.to.file("aa_bf.pdf", plot.aa.bf())
    plot.to.file("sngfr_lookup_base.pdf", plot.sngfr.lookup(), width=10)
    plot.to.file("sngfr_lookup_glom.pdf",
                 plot.sngfr.lookup("sngfr_lookup_glom.ssv"), width=10)

}

plot.figures(set.theme=FALSE)
