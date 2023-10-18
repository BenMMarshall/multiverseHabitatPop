# Download package tarball from CRAN archive

url <- "https://cran.r-project.org/src/contrib/Archive/rgeos/rgeos_0.6-4.tar.gz"
pkgFile <- "rgeos_0.6-4.tar.gz"
download.file(url = url, destfile = pkgFile)

# Expand the zip file using whatever system functions are preferred

# look at the DESCRIPTION file in the expanded package directory

# Install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
unlink(pkgFile)
