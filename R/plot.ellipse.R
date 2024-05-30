#' Plot method for "ellipse" class
#'
#' @param x An object of class "ellipse"
#' @param ... Other arguments passed to or from other methods
#' @export plot.ellipse
#' @export
plot.ellipse <- function(x, ...) {
  dimPlot <- length(x[[1]])
  varNames <- names(x[[1]])

  colList <- utils::combn(1:dimPlot, 2, simplify = FALSE)

  triangle_line_lower <- function(k, n){1:(k) + n*(k)}
  triangle_line_upper <- function(k, n){(k + 1):n + n*(k - 1)}

  lower_panel_index <- function(dimPlot){
    Reduce(c, Map(function(k, n) triangle_line_lower(k, n), 1:(dimPlot - 1), dimPlot ))
  }

  upper_panel_index <- function(dimPlot){
    Reduce(c, Map(function(k, n) triangle_line_upper(k, n), 1:(dimPlot - 1), dimPlot ))
  }
  lowerPanelIndex <- lower_panel_index(dimPlot)
  upperPanelIndex <- upper_panel_index(dimPlot)

  dataToPlot <- Map(function(cols){ellipse_points(x, cols)}, colList)

  diagSeq <- seq(1, dimPlot^2, dimPlot + 1)
  op <- graphics::par(mfrow = c(dimPlot, dimPlot))
  invisible(

  Map(function(dataP){
    if(dataP %in% diagSeq ){
      m <- match(dataP, diagSeq)
      graphics::plot.new()
      graphics::text(0.5, 0.5, varNames[m],cex = 2)
    } else if(dataP %in% lowerPanelIndex) {
      m <- match(dataP, lowerPanelIndex)

      plot(data = dataToPlot[[m]],
           x ~ y,
           type = "l",
           xlab = varNames[ colList[[m]][1] ],
           ylab = varNames[ colList[[m]][2] ]
      )

    } else if(dataP %in% upperPanelIndex) {
      m <- match(dataP, upperPanelIndex)

      plot(data = dataToPlot[[m]],
           y ~ x,
           type = "l",
           xlab = varNames[ colList[[m]][2] ],
           ylab = varNames[ colList[[m]][1] ]
      )
    }
  }, 1:dimPlot^2 )
  )
  op <- graphics::par(mfrow = c(1, 1))
}
