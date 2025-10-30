savepkgs <- function(saveloc)
{
  instpkg <- installed.packages()
  cranpkg <- available.packages()
  iscran <- instpkg[,'Package'] %in% cranpkg[,'Package']
  pkgdata <- data.frame(package=instpkg[,'Package'], iscran=iscran, instvers=instpkg[,'Version'])
  write.csv(pkgdata, saveloc)
}

restore_pkgs <- function(pkglistfile)
{
  pkgtbl <- read.csv(pkglistfile, stringsAsFactors = FALSE)
  pkglist <- pkgtbl$package[pkgtbl$iscran]
  ## install the packages we can get from CRAN
  install.packages(pkglist, dependencies=TRUE, type='source', Ncpus=8)
  ## Message about the ones we can't get from CRAN
  nocran <- pkgtbl$package[!pkgtbl$iscran]
  nocranstr <- paste(nocran, collapse=', ')
  message('The following packages were not available from CRAN: ', nocranstr)
  ## Return the list of not-installed packages for the user's consideration
  nocran
}