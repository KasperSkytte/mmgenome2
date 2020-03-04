.onAttach <- function(lib, pkg) {
  options(scipen = 6)
  # Check for new github version (master branch).
  if (!interactive()) {
    return()
  } else {
    installed_version <- as.character(utils::packageVersion(pkg))
    packageStartupMessage(
      paste0(
        "This is ",
        pkg,
        " version ",
        installed_version,
        ". Great documentation is available at the mmgenome2 website: https://kasperskytte.github.io/mmgenome2/"),
      appendLF = TRUE)
    gitHubUser <- "kasperskytte"
    errMsg <- "Couldn't reach GitHub to check for new version just now."
    tryCatch(
      {
        DESCRIPTION <- readLines(
          paste0(
            "https://raw.githubusercontent.com/",
            gitHubUser,
            "/",
            pkg,
            "/master/DESCRIPTION"
          )
        )
        remote_version <- gsub("Version:\\s*", "", DESCRIPTION[grep("Version:", DESCRIPTION)])
        if (installed_version < remote_version) {
          packageStartupMessage(
            paste0(
              "\nNewer version of ", 
              pkg, 
              " (", 
              remote_version, 
              ") is available! Install the latest version with the following command (copy/paste): \nremotes::install_github(\"kasperskytte/mmgenome2\")"),
            appendLF = TRUE
          )
        }
      },
      error = function(e) {
        packageStartupMessage(errMsg, appendLF = TRUE)
      },
      warning = function(e) {
        packageStartupMessage(errMsg, appendLF = TRUE)
      }
    )
  }
}
