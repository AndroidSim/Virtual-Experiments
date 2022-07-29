#! /bin/bash setR

if (length(argv) <= 0) {
  print("Usage: ei-plots.r <data files *.csv>")
  print("  Data should be in the format:")
  print("  Indep Name, Dep1 Name, Dep2 Name, Dep3 Name")
  print("  indep val, dep1 val, dep2 val, dep3 val")
  print("               ...")
  quit()
}

for (f in argv) {
  # get filename root
  nb <- unlist(strsplit(f, "\\.csv"))[1]
  # read the data
  d <- read.csv(file = f)

  ind <- colnames(d)[1]
  for (dep in colnames(d)[-1]) {
    # Start the graphics device
    #png(paste("graphics/",nb, "-", dep, ".png",sep=""), width=864, height=432)
    svg(paste("graphics/",nb, "-", dep, ".svg",sep=""), width=8, height=4)
    plot(d[[ind]], d[[dep]], main=nb, type="l", xlab=ind, ylab=dep)
    dev.off() # otherwise "too many open devices"
  }

}
q()
