library(devtools)

# creating documentation (i.e. the Rd files in man/)
devtools::document()

# checking documentation
devtools::check_man()

# running tests
devtools::test()

# building tarball (e.g. oski_0.1.tar.gz)
devtools::build()

# checking install
devtools::install()

devtools::use_vignette("introduction")





#########
library(devtools)
# creating documentation (i.e. the Rd files in man/)
devtools::document()
# checking documentation
devtools::check_man() # run tests
devtools::test()
# checking documentation
devtools::build_vignettes()
# building tarball (e.g. oski_0.1.tar.gz)
devtools::build()
# checking install
devtools::install()
