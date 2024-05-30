plot_test <- function(ellipse){
  dimPlot <- length(ellipse[[1]])
  varNames <- names(ellipse[[1]])

  colList <- combn(1:dimPlot, 2, simplify = FALSE)

  #variables <- Map(function(x){paste0(x, collapse = "_")},  colList)

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

  dataToPlot <- Map(function(cols){ellipse_points(ellipse_par_example, cols)}, colList)

  diagSeq <- seq(1, dimPlot^2, dimPlot + 1)
  op <- par(mfrow = c(dimPlot, dimPlot))

  Map(function(dataP){
    if(dataP %in% diagSeq ){
      m <- match(dataP, diagSeq)
      plot.new()
      text(0.5, 0.5, varNames[m],cex = 2)
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

}

M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
Mshp <- terra::unwrap(readr::read_rds(M_path))


stack_1_12 <- get_example_data("stack_1_12_19")
library(nicher)
pars <- get_ENM_par(rawSpOccPnts, stack_1_12, Mshp, method = "vmmin")

plot_test(pars)
