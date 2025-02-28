library(ggplot2)
library(RColorBrewer)

m_f <- snakemake@input[['m']]
outplot_f <- snakemake@output[['outplot']]
log_f <- snakemake@log[['log']]

size <- 8


axis_names <- list(
  H2AX_induction= 'γ-2AX induction post-treatment',
  dvw3= 'Relative tumor growth (%t baseline)'
)

x_lim <- list(
  H2AX_induction= c(0, 60)
)

y_lim <- list(
  dvw3 = c(-100, 300)
)


m <- read.table(m_f, sep="\t", stringsAsFactors = F)

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
#showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
#showtext_auto(enable = TRUE)

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- (size) * death_conversion_dpi96

# arial family gave:
#Error in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  : 
#                              invalid font type
                            

unmute_theme <- theme(
  text = element_text(size = textSize), #, family='Arial'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size=0.508/0.564), # origin of this ratio is honestly not known, empirical
  axis.ticks = element_line(color = "black", size=0.508/0.564),
  axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),
  panel.background = element_blank()
)

guess_ticks <- function(values, nticks=5, fixed_max=NULL, fixed_min=0) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- ceiling(vmax)
  } else {
    round_max <- fixed_max
  }
  my_breaks <- seq(fixed_min, round_max, length.out=nticks)
  return(my_breaks)
}


plot_lmsmooth_info <- function(xvar, yvar, md) {
  pear <- cor.test(md[,xvar], md[,yvar])
  title <- paste0("pearson=", signif(pear$estimate, 3), " P=", signif(pear$p.value, 3))
  xl <- x_lim[[xvar]]
  yl <- y_lim[[yvar]]
  y_breaks <- guess_ticks(md[,yvar], fixed_min=yl[1], fixed_max=yl[2])
  x_breaks <- guess_ticks(md[,xvar], fixed_min=xl[1], fixed_max=xl[2])
  p <- ggplot(data=md, aes_string(x=xvar, y=yvar))+geom_point(size=1)+theme_bw()+geom_smooth(method="lm", se=FALSE, size=1)+
    unmute_theme+theme(legend.position="none") + xlab(xvar)+ylab(yvar)+
    scale_y_continuous(breaks=y_breaks, limits=yl, expand=c(0,0))+
    scale_x_continuous(breaks=x_breaks, limits=xl, expand=c(0,0))+
    xlab(axis_names[[xvar]])+
    ylab(axis_names[[yvar]])
  return(list(p, title))
}

p <- plot_lmsmooth_info(snakemake@wildcards[['axisX']], snakemake@wildcards[['axisY']], m)

sink(log_f)
print(p[[2]])
sink()

#
#ggsave(p, file=plot_cor_png, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

ggsave(p[[1]], file=outplot_f, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

save.image('pippo.Rdata')