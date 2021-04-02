#' @title Merge mm object with more data
#'
#' @description Internal function used in mmload. \code{y} can be a named vector, dataframe, or a list containing any of these types of data to be merged with \code{x}.
#'
#' @param x mm object
#' @param y Object to merge with \code{x}
#' @param type Character string defining the type of data being merged
#'
#' @importFrom dplyr left_join intersect filter
#' @importFrom tibble enframe
#'
#' @return A tibble.
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
mmmerge <- function(x, y, type) {
  # must be a data frame, named atomic vector, or a list of data frames and/or named vectors
  if (any(class(y) %in% c("list", "data.frame", "tbl", "tbl_df")) | is.atomic(y) | is.factor(y)) {
    # wrap non-lists in a list to work with the for-loop
    if (any(!class(y) %in% "list")) {
      y <- list(y)
    }
    # merge mm with each element in the provided list
    for (i in 1:length(y)) {
      string <- if (length(y) == 1) paste0("'", type, "'") else paste0("'", type, "'", " element ", i)
      if (is.factor(y[[i]])) {
        y[[i]] <- as.character(y[[i]])
      }
      # first column must be the sequence names
      if (any(class(y[[i]]) %in% c("data.frame", "tbl", "tbl_df")) & length(y[[i]]) < 2) {
        stop(paste0(string, " not accepted: Data frames must contain at least 2 columns where the first column contains the sequence names exactly matching those of the assembly."), call. = FALSE)
      }

      # vectors must be named to be able to merge with mm
      if (is.atomic(y[[i]]) & is.null(names(y[[i]]))) {
        stop(paste0(string, " not accepted: The vector is not a named vector. The vector elements must be named by sequence names exactly matching those of the assembly."), call. = FALSE)
      }

      # column names are preserved from data frames, but not from vectors. Use names of the provided list, or else create a dummy name
      if (is.atomic(y[[i]]) & !is.null(names(y[[i]]))) {
        y[[i]] <- tibble::enframe(y[[i]], name = "scaffold", value = if ((is.null(names(y)[[i]]) || names(y)[[i]] == "")) paste0(type, i) else names(y)[[i]])
      }

      # if y[[i]] is a 2 column data frame, name the resulting column in x by the name of the dataframe as provided in the list, if any.
      # Otherwise if the y[[i]] has no name in the list, keep the current column name
      if (any(class(y[[i]]) %in% c("data.frame", "tbl", "tbl_df")) & length(y[[i]]) == 2 & !(is.null(names(y)[[i]]) || names(y)[[i]] == "")) {
        colnames(y[[i]])[2] <- names(y)[[i]]
      } else {
        colnames(y[[i]])[2] <- colnames(y[[i]])[2]
      }

      # merge x and y[[i]] by scaffold
      colnames(y[[i]])[1] <- "scaffold" # first columns must have same name
      y[[i]][1] <- lapply(y[[i]][1], as.character) # and must be character
      if (any(duplicated(y[[i]][1]))) {
        stop(paste0("All scaffold names must be unique. Duplicate scaffold names found in ", string, "."), call. = FALSE)
      }
      sharedScaffolds <- dplyr::intersect(x$scaffold, y[[i]][["scaffold"]]) # which scaffolds are shared between x and y[[i]]

      # print missing or excess scaffolds between x and y[[i]]
      if (!all(x$scaffold %in% y[[i]][["scaffold"]])) {
        # missingScaffolds <- dplyr::filter(x, !scaffold %in% sharedScaffolds)[[1]]
        warning(paste0(
          "Only ",
          length(sharedScaffolds),
          " of all ",
          length(x$scaffold),
          " scaffolds in the assembly match in ",
          string,
          "."
          # , " The following "
          # , length(missingScaffolds),
          # ," scaffolds are missing:\n\""
          # ,paste(missingScaffolds,
          #       collapse = "\", \""),
          # "\""
        ),
        call. = FALSE
        )
      } else if (!all(y[[i]][["scaffold"]] %in% x$scaffold)) {
        excessScaffolds <- filter(y[[i]], !scaffold %in% sharedScaffolds)[[1]]
        warning(paste0(
          string,
          " contains ",
          length(excessScaffolds),
          " more scaffolds than the assembly. "
          # ,"The following ",
          # ," scaffolds have not been loaded:\n\"",
          # paste(excessScaffolds,
          #       collapse = "\", \""),
          # "\""
        ),
        call. = FALSE
        )
      } else if (!any(x$scaffold %in% y[[i]][["scaffold"]])) {
        # no match sucks
        stop("No scaffold names match between the assembly and ", string, ". ", call. = FALSE)
      }
      x <- dplyr::left_join(x,
        y[[i]],
        by = "scaffold"
      )
    }
  } else {
    stop("Data must be provided as a data frame, named vector, or a list of multiple data frames and/or named vectors.", call. = FALSE)
  }
  return(x)
}

#' @title Locate points in a ggplot2 plot
#'
#' @description Internal function used within mmplot and mmplot_network
#'
#' @param plot ggplot2 object.
#' @param x_scale x axis scale, \code{"sqrt"} or \code{"log10"}.
#' @param y_scale y axis scale, \code{"sqrt"} or \code{"log10"}.
#'
#' @return A data frame with the x/y coordinates of the mousepositions clicked in the ggplot2 plot.
#'
#' @import ggplot2
#' @importFrom shiny actionButton div fillPage icon observeEvent p plotOutput reactiveValues renderPlot runApp shinyApp stopApp
#' @importFrom clipr write_clip
#' @importFrom rstudioapi viewer
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Rasmus Hansen Kirkegaard \email{rhk@@bio.aau.dk}
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
mmlocator <- function(plot, x_scale = NULL, y_scale = NULL) {
  app <- shinyApp(
    ui = fillPage(
      padding = c(5, 5, 50),
      div(
        style = "display: inline-block;",
        actionButton(
          inputId = "undo",
          label = "Undo",
          icon = icon("undo")
        ),
        actionButton(
          inputId = "stop",
          label = "Finish",
          icon = icon("check")
        )
      ),
      p(),
      plotOutput("plot",
        click = "plot_click",
        height = "100%",
        width = "100%"
      )
    ),
    server = function(input, output, session) {
      clickData <- reactiveValues(
        x = numeric(),
        y = numeric()
      )

      observeEvent(input$plot_click, {
        new_x <- input$plot_click$x
        new_y <- input$plot_click$y

        # adjust points if square-root scale have been applied in the plot
        if (!is.null(x_scale)) {
          if (x_scale == "sqrt") {
            new_x <- new_x^2
          }
        }
        if (!is.null(y_scale)) {
          if (y_scale == "sqrt") {
            new_y <- new_y^2
          }
        }

        # add the click coordinates to clickData
        clickData[["x"]] <- append(clickData[["x"]], new_x)
        clickData[["y"]] <- append(clickData[["y"]], new_y)

        # Update selection object in global environment
        df <- get(".current_selection", envir = globalenv())
        assign(".current_selection",
          rbind(
            df,
            data.frame(
              x = new_x,
              y = new_y
            )
          ),
          envir = globalenv()
        )
        # message("point clicked:\n",
        #        "\t", colnames(df)[1], " = ", round(df[nrow(df), 1], 3),
        #        "\t", colnames(df)[2], " = ", round(df[nrow(df), 2], 3))
      })

      observeEvent(input$undo, {
        # remove the last point clicked
        lastPoint <- length(clickData[["x"]])
        clickData[["x"]] <- clickData[["x"]][-lastPoint]
        clickData[["y"]] <- clickData[["y"]][-lastPoint]
        # Update selection object in global environment
        df <- get(".current_selection", envir = globalenv())
        assign(".current_selection", df[-lastPoint, ], envir = globalenv())
      })

      observeEvent(input$stop, {
        stopApp()
      })

      output$plot <- renderPlot({
        p <- plot +
          geom_point(
            data = data.frame(
              x = clickData[["x"]],
              y = clickData[["y"]]
            ),
            aes_string(
              x = "x",
              y = "y"
            ),
            color = "black",
            size = 2,
            inherit.aes = FALSE,
            na.rm = TRUE
          ) +
          geom_polygon(
            data = data.frame(
              x = clickData[["x"]],
              y = clickData[["y"]]
            ),
            aes_string(
              x = "x",
              y = "y"
            ),
            fill = NA,
            size = 0.5,
            lty = 2,
            color = "darkred",
            inherit.aes = FALSE,
            na.rm = TRUE
          )
        return(p)
      })
    },
    onStart = function() {
      df <- data.frame(x = numeric(), y = numeric())
      assign(".current_selection", df, envir = globalenv())
    }
  )
  suppressWarnings(
    runApp(app,
      quiet = TRUE,
      launch.browser = rstudioapi::viewer
    )
  )
  df <- get(".current_selection", df, envir = globalenv())
  colnames(df) <- c(
    gsub("~", "", deparse(plot[["mapping"]][["x"]])),
    gsub("~", "", deparse(plot[["mapping"]][["y"]]))
  )
  assign(".current_selection", df, envir = globalenv())
  if (nrow(df) > 0) {
    selection <- paste0(
      "data.frame(",
      colnames(df[1]),
      " = ",
      paste0(round(df[1], 3)),
      ",\n           ",
      colnames(df[2]),
      " = ",
      paste0(round(df[2], 3)),
      ")"
    )
    message(paste0("The following selection has been copied to clipboard:\n", selection))
    clipr::write_clip(selection)
  }
  return(df)
}

#' @title Check for installed package(s)
#' @description Returns an error if required package(s) are not installed. Mostly used for checking whether packages listed under the Suggests field in the DESCRIPTION file are installed.
#'
#' @param pkgs Character vector with package(s) to check for.
#' @param msg Optionally additional text appended (with \code{paste0}) to the default error message.
#'
#' @return Returns error and message if not installed, otherwise \code{invisible(TRUE)}
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
checkReqPkgs <- function(pkgs, msg = "") {
  stopifnot(is.character(pkgs), length(pkgs) > 0L, nchar(pkgs) > 0, is.character(msg))
  missingPkgs <- character()
  # loop over pkgs appending missing ones to vector, otherwise load
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missingPkgs <- c(missingPkgs, pkg)
    } else {
      require(pkg, quietly = TRUE, character.only = TRUE)
    }
  }
  if (length(missingPkgs) > 0L) {
    stop(
      paste0(
        "The following package",
        if (length(missingPkgs) > 1L) {
          "s are "
        } else {
          " is "
        },
        "required but not installed: ",
        paste0(pkgs, collapse = ", "),
        msg
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# This function is sourced from the DescTools CRAN package
# version 0.99.24 for pretty printing numbers
StrAlign <- function(x, sep = "\\r") {
  id.na <- is.na(x)
  if (length(grep("\\", sep, fixed = TRUE)) == 0) {
    idx <- !grepl(x = x, pattern = sep, fixed = TRUE)
    x[idx] <- paste(x[idx], sep, sep = "")
  }
  if (sep == "\\c") {
    return(StrPad(x,
      width = max(nchar(x), na.rm = TRUE),
      pad = " ", adj = "center"
    ))
  }
  x <- StrPad(x, max(nchar(x), na.rm = TRUE))
  if (sep == "\\l") {
    return(sub("(^ +)(.+)", "\\2\\1", x))
  }
  if (sep == "\\r") {
    return(sub("(.+?)( +$)", "\\2\\1", x))
  }
  bef <- substr(x, 1, StrPos(x, sep, fix = TRUE))
  aft <- substr(x, StrPos(x, sep, fix = TRUE) + 1, nchar(x))
  aft <- substr(aft, 1, max(nchar(StrTrim(aft, method = "right"))))
  res <- paste(replace(StrPad(bef, max(nchar(bef), na.rm = TRUE),
    " ",
    adj = "right"
  ), is.na(bef), ""), replace(
    StrPad(aft,
      max(nchar(aft), na.rm = TRUE), " ",
      adj = "left"
    ), is.na(aft),
    ""
  ), sep = "")
  res[id.na] <- NA
  if (length(grep("\\", sep, fixed = TRUE)) == 0) {
    res[idx] <- gsub(sep, " ", res[idx], fixed = TRUE)
  }
  return(res)
}

# This function is sourced from the DescTools CRAN package
# version 0.99.24 for pretty printing numbers
StrPad <- function(x, width = NULL, pad = " ", adj = "left") {
  .pad <- function(x, width, pad = " ", adj = "left") {
    if (is.na(x)) {
      return(NA)
    }
    mto <- match.arg(adj, c("left", "right", "center"))
    free <- max(0, width - nchar(x))
    fill <- substring(paste(rep(pad, ceiling(free / nchar(pad))),
      collapse = ""
    ), 1, free)
    if (free <= 0) {
      x
    } else if (mto == "left") {
      paste(x, fill, sep = "")
    } else if (mto == "right") {
      paste(fill, x, sep = "")
    } else {
      paste(substring(fill, 1, free %/% 2), x, substring(
        fill,
        1 + free %/% 2, free
      ), sep = "")
    }
  }
  if (is.null(width)) {
    width <- max(nchar(x), na.rm = TRUE)
  }
  lgp <- Recycle(
    x = x, width = width, pad = pad,
    adj = adj
  )
  sapply(1:attr(lgp, "maxdim"), function(i) {
    .pad(
      lgp$x[i],
      lgp$width[i], lgp$pad[i], lgp$adj[i]
    )
  })
}

# This function is sourced from the DescTools CRAN package
# version 0.99.24 for pretty printing numbers
StrPos <- function(x, pattern, pos = 1, ...) {
  pos <- rep(pos, length.out = length(x))
  x <- substr(x, start = pos, stop = nchar(x))
  i <- as.vector(regexpr(pattern = pattern, text = x, ...))
  i[i < 0] <- NA
  return(i)
}

# This function is sourced from the DescTools CRAN package
# version 0.99.24 for pretty printing numbers
StrTrim <- function(x, pattern = " \t\n", method = "both") {
  switch(match.arg(arg = method, choices = c(
    "both", "left",
    "right"
  )),
  both = {
    gsub(
      pattern = gettextf("^[%s]+|[%s]+$", pattern, pattern),
      replacement = "", x = x
    )
  },
  left = {
    gsub(
      pattern = gettextf("^[%s]+", pattern), replacement = "",
      x = x
    )
  },
  right = {
    gsub(
      pattern = gettextf("[%s]+$", pattern), replacement = "",
      x = x
    )
  }
  )
}

# This function is sourced from the DescTools CRAN package
# version 0.99.24 for pretty printing numbers
Recycle <- function(...) {
  lst <- list(...)
  maxdim <- max(unlist(lapply(lst, length)))
  # recycle all params to maxdim
  res <- lapply(lst, rep_len, length.out = maxdim)
  attr(res, "maxdim") <- maxdim

  return(res)
}
