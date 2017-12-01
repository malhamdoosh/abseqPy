# https://gist.github.com/stevenworthington/3178163 with minor modifications
depsInstaller <- function(pkgs) {
    new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
    if (length(new.pkgs) > 0) {
        install.packages(new.pkgs, repos="http://cran.rstudio.com/", dependencies = TRUE)
    }
}

deps <- c("ggplot2", "RColorBrewer", "circlize", "reshape2", "VennDiagram", "plyr")
depsInstaller(deps)