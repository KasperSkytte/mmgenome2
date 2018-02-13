.onAttach <- function(lib, pkg)  {
  local_version <- utils::packageVersion("mmgenome2")
  packageStartupMessage("This is ", pkg, " version ", local_version, ". Great documentation is available at the mmgenome2 website: https://madsalbertsen.github.io/mmgenome2/\n", appendLF = TRUE)
  options(scipen = 6)
}