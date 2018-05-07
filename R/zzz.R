.onAttach <- function(lib, pkg)  {
  options(scipen = 6)
  #Check for new github release version. (Not master branch version, release version!)
  if (!interactive()) {
    return()
  } else {
    local_version <- utils::packageVersion("mmgenome2")
    packageStartupMessage("This is ", pkg, " version ", local_version, ". Great documentation is available at the mmgenome2 website: https://kasperskytte.github.io/mmgenome2/", appendLF = TRUE)
    if(requireNamespace("remotes", quietly = TRUE)) {
      tryCatch({
        github_ref <- remotes:::github_resolve_ref(
          remotes::github_release(), 
          remotes:::parse_git_repo("kasperskytte/mmgenome2@*release"))$ref
        github_version <- package_version(gsub("v", "", github_ref))
        if(local_version < github_version) {
          packageStartupMessage(
            "\nNew release of ", pkg, " (", github_version, ") is available! Install the latest release with (copy/paste): \nremotes::install_github(\"kasperskytte/mmgenome2@*release\")\n\nRead the release notes at: https://github.com/kasperskytte/mmgenome2/releases/tag/", github_version)
        }
      }, error=function(e) {
        packageStartupMessage("\nCan't reach GitHub to check for new releases just now. Trying again next time.")
      })
    }
  }
}