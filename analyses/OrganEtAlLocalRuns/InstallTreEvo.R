local({r <- getOption("repos")
       r["CRAN"] <- "http://cloud.r-project.org/" 
       options(repos=r)
})
#install.packages("devtools")
install.packages(c("pbkrtest", "geiger", "pls", "stats", "corpcor", "coda", "foreach", "doMC", "rgl", "partitions", "car", "mvtnorm", "methods", "grDevices", "graphics", "phytools", "ape"), type="source")
devtools::install_github("bomeara/treevo", type="source")

