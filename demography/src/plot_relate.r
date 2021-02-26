## Plot Relate demography.
## JMIB [jmi.barth@gmail.com], 2020-06-10

##Install packages if necessary and load them.
install.packages("remotes",repos = "http://cran.us.r-project.org")
remotes::install_github("leospeidel/relater")
library(relater)
library(ggplot2)

##Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
coal_file_name <- args[1]
generation_time <- as.numeric(args[2])
pdf_file_name <- args[3]

##Read coalescence rate input file.
coal <- read.coal(coal_file_name)
coal$popsize <- 0.5/coal$haploid.coalescence.rate #diploid effective population size is the 0.5* inverse coalescence rate
coal$epoch.start <- generation_time * coal$epoch.start #multiply epochs times to scale to years
coal$within <- c("across groups","within group")[(coal$group1 == coal$group2)+1] #add a column on whether coalescence rate is within or across groups

##Split within and across rates.
coal_split <- split(coal, coal$within) #split data frame in separate tables for within and across rates

##Plot the inverse of all within-population coalescence rates as the population size over time.
pdf(file=pdf_file_name, height=5, width=10)
ggplot(coal_split$`within group`) +
  geom_step(aes(x = epoch.start, y = popsize, color = group2, linetype = within, alpha=0.8)) +
  scale_x_continuous(trans = "log10", limit = c(1e2,1e7)) +
  scale_y_continuous(trans = "log10", limit = c(1e2,1e8)) +
  annotation_logticks(sides = "lb") +
  xlab("years ago") +
  ylab("diploid effective population size")
dev.off()
