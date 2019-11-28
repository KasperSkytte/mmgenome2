.onAttach <- function(lib, pkg) {
  options(scipen = 6)
  # Check for new github release version. (Not master branch version, release version!)
  if (!interactive()) {
    return()
  } else {
    installed_version <- as.character(utils::packageVersion(pkg))
    gitHubUser <- "kasperskytte"
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
        startupMsg <- ""
      },
      error = function(e) {
        startupMsg <- "\nCan't reach GitHub to check for new version just now. Trying again next time."
        remote_version <- "0"
      }
    )

    if (installed_version < remote_version) {
      startupMsg <- paste0("\nNew version of ", pkg, " (", remote_version, ") is available! Install the latest version with the following command (copy/paste): \nremotes::install_github(\"kasperskytte/mmgenome2\")")
    }
    startupMsg <- paste0(
      "This is ",
      pkg,
      " version ",
      installed_version,
      if (installed_version >= remote_version) {
        " (up to date)"
      },
      ". Great documentation is available at the mmgenome2 website: https://kasperskytte.github.io/mmgenome2/",
      startupMsg
    )
    packageStartupMessage(startupMsg, appendLF = TRUE)
  }
}
