# bcgp
bcgp

The function rstantools::rstan_package_skeleton() isn't working right. It calls usethis::create_package(), which in turn tells me my project is nested in another package. Continuing to create the package from there creates the package fine, but then there are thousands of irrelevant documents in the Git pane. Instead of messing with it more, I've just created this package and will pull it when I want to write a package involving compiled Stan code.
