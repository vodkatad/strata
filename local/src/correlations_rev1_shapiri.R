library(cowplot)
library(ggpubr)
library(ggplot2)

d <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/figures/IHC_m29_merged.tsv', sep="\t", header=T, row.names=1)
reci <- d$recist
d$recist <- NULL
shapi <- apply(d, 2, function(x) { shapiro.test(x)$p.value})
round(shapi, digits=3)


d$recist <- reci
sp <- ggscatter(d, x = "H2AX_induction", y = "dvw3",
                color = "recist", palette = "jco",
                size = 3, alpha = 0.6)+border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(d, "H2AX_induction", fill="darkolivegreen3")
yplot <- ggdensity(d, "dvw3",  fill="darkolivegreen3")+ rotate()
# Cleaning the plots
sp <- sp# + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot

plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

qqnorm(d$POLD1)
qqline(d$POLD1)

qqnorm(d$H2AX_induction)
qqline(d$H2AX_induction)

qqnorm(scale(d$H2AX_induction))
qqline(scale(d$H2AX_induction))


sp <- ggscatter(d, x = "Illumina", y = "POLD1",
                color = "recist", palette = "jco",
                size = 3, alpha = 0.6)+border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(d, "Illumina", fill="darkolivegreen3")
yplot <- ggdensity(d, "POLD1",  fill="darkolivegreen3")+ rotate()
# Cleaning the plots
sp <- sp# + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot

plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

pdf('~/qqplot.pdf')
qqnorm(log2(d$H2AX_induction))
qqline(log2(d$H2AX_induction))
graphics.off()

d$loggamma <- log2(d$H2AX_induction)
sp <- ggscatter(d, x = "loggamma", y = "dvw3",
                color = "recist", palette = "jco",
                size = 3, alpha = 0.6)+border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(d, "loggamma", fill="darkolivegreen3")
yplot <- ggdensity(d, "dvw3",  fill="darkolivegreen3")+ rotate()
# Cleaning the plots
sp <- sp# + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot

plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))


pcordens <- function(myd, xlab, ylab, collab) {
  sp <- ggscatter(myd, x = xlab, y = ylab,
                  color = collab, palette = "jco",
                  size = 3, alpha = 0.6)+border()                                         
  # Marginal density plot of x (top panel) and y (right panel)
  xplot <- ggdensity(myd, xlab, fill="darkolivegreen3")
  yplot <- ggdensity(myd, ylab,  fill="darkolivegreen3")+ rotate()
  # Cleaning the plots
  sp <- sp# + rremove("legend")
  yplot <- yplot + clean_theme() + rremove("legend")
  xplot <- xplot + clean_theme() + rremove("legend")
  # Arranging the plot using cowplot
  
  plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
            rel_widths = c(2, 1), rel_heights = c(1, 2))
}
