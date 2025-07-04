#' @name scoreHVT
#' @title Score which cell each point in the test dataset belongs to.
#' @description This function scores each data point in the test dataset based on a trained hierarchical Voronoi tessellations model. 
#' @param dataset Data frame. A data frame which to be scored. Can have categorical columns if `analysis.plots` are required.
#' @param hvt.results.model List. A list obtained from the trainHVT function 
#' @param child.level Numeric. A number indicating the depth for which the heat map is to be plotted. 
#' @param mad.threshold Numeric. A numeric value indicating the permissible Mean Absolute Deviation.
#' @param line.width Vector. A vector indicating the line widths of the tessellation boundaries for each layer. 
#' @param color.vec Vector. A vector indicating the colors of the tessellation boundaries at each layer. 
#' @param normalize Logical. A logical value indicating if the dataset should be normalized. When set to TRUE,
#'  the data (testing dataset) is standardized by ‘mean’ and ‘sd’ of the training dataset referred from the trainHVT(). 
#'  When set to FALSE, the data is used as such without any changes.
#' @param distance_metric Character. The distance metric can be L1_Norm(Manhattan) or L2_Norm(Eucledian). L1_Norm is selected by default.
#' The distance metric is used to calculate the distance between an n dimensional point and centroid.
#'  The distance metric can be different from the one used during training.
#' @param error_metric Character. The error metric can be mean or max. max is selected by default. 
#' max will return the max of m values and mean will take mean of m values where
#' each value is a distance between a point and centroid of the cell.
#' @param yVar Character. A character or a vector representing the name of the dependent variable(s)
#' @param analysis.plots Logical. A logical value indicating that the scored plot should be plotted or not. If TRUE, 
#' the identifier column(character column) name should be supplied in `names.column` argument. The output will
#' be a 2D heatmap plotly which gives info on the cell id and the observations of a cell.
#' @param names.column Character. A character or a vector representing the name of the identifier column/character column.
#' @returns Dataframe containing scored data, plots and summary
#' @author Shubhra Prakash <shubhra.prakash@@mu-sigma.com>, Sangeet Moy Das <sangeet.das@@mu-sigma.com> ,
#' Vishwavani <vishwavani@@mu-sigma.com>
#' @seealso \code{\link{trainHVT}} \cr \code{\link{plotHVT}}
#' @keywords Scoring
#' @importFrom magrittr %>%
#' @examples
#' data("EuStockMarkets")
#' dataset <- data.frame(date = as.numeric(time(EuStockMarkets)),
#'                      DAX = EuStockMarkets[, "DAX"],
#'                      SMI = EuStockMarkets[, "SMI"],
#'                      CAC = EuStockMarkets[, "CAC"],
#'                      FTSE = EuStockMarkets[, "FTSE"])
#' rownames(EuStockMarkets) <- dataset$date
#' # Split in train and test
#' train <- EuStockMarkets[1:1302, ]
#' test <- EuStockMarkets[1303:1860, ]
#' #model training
#' hvt.results<- trainHVT(train,n_cells = 60, depth = 1, quant.err = 0.1,
#'                       distance_metric = "L1_Norm", error_metric = "max",
#'                       normalize = TRUE,quant_method = "kmeans")
#' scoring <- scoreHVT(test, hvt.results)
#' data_scored <- scoring$scoredPredictedData
#' @export scoreHVT


scoreHVT <- function(dataset,
                       hvt.results.model,
                       child.level = 1,
                       mad.threshold = 0.2,
                       line.width = 0.6,
                       color.vec = c("navyblue", "slateblue", "lavender"),
                       normalize = TRUE,
                       distance_metric = "L1_Norm",
                       error_metric = "max",
                       yVar = NULL,
                       analysis.plots = FALSE,
                       names.column = NULL
                       ) {
  
  suppressWarnings({
  set.seed(300)
  sum_n = NULL
  requireNamespace("dplyr")
  requireNamespace("purrr")
  requireNamespace("data.table")
  
  if (any(is.na(dataset))) {
    stop("Input dataframe contains missing (NA) values. Please handle missing data before using this function (e.g., via na.omit(), na.fill(), or imputation).")
  } 
  

  if (!("Cell.ID" %in% colnames(hvt.results.model[[3]]$summary))) {
    hvt.results.model[[3]]$summary <- getCellId(hvt.results = hvt.results.model)
  }
 
  if (analysis.plots == TRUE && is.null(names.column)) {
   stop("names.column is not defined")
  }
  
  
  hvt.results.model[[3]]$summary <- cbind(hvt.results.model[[3]]$summary, centroidRadius = unlist(hvt.results.model[[3]]$max_QE))

  numeric_columns <- sapply(as.data.frame(dataset), is.numeric)
  data <- dataset[, numeric_columns, drop = FALSE]

  summary_list <- hvt.results.model[[3]]
  train_colnames <- names(summary_list[["nodes.clust"]][[1]][[1]])

    if (!all(train_colnames %in% colnames(data))) {
    stop("Not all training columns are part of test dataset")
  }
  
  data_structure_score  <- dim(data)
  
  if (!all(is.na(summary_list$scale_summary)) && normalize == TRUE) {
    scaled_test_data <- scale(
      data[, train_colnames],
      center = summary_list$scale_summary$mean_data[train_colnames],
      scale = summary_list$scale_summary$std_data[train_colnames]
    )
  } else {
     scaled_test_data <- data[, train_colnames]
 
  }

  colnames(scaled_test_data) <- train_colnames
 
  level <- child.level

  
  if (!is.null(yVar)) {
    yVardf <- data[, yVar]
    if (length(yVar) != 1) {
      colnames(yVardf) <- paste0("Scored.", yVar)
    }
  }
   


  find_path <- function(data_vec, centroid_data) {
  
    if (distance_metric == "L1_Norm") {
      centroidDist <- which.min(colSums(abs(centroid_data - data_vec), na.rm = TRUE))
      Quant.Error <- (colSums(abs(centroid_data - data_vec), na.rm = TRUE))[centroidDist]
    } else {
      centroidDist <- which.min(sqrt(colSums((centroid_data - data_vec)^2, na.rm = TRUE)))
      Quant.Error <- sqrt(colSums((centroid_data - data_vec)^2, na.rm = TRUE))[centroidDist]
    }
    return(data.frame("Index" = centroidDist, "Quant.Error" = Quant.Error / length(train_colnames)))
  }

 
  newdfMapping <- summary_list$summary

  innermostCells2 <- newdfMapping %>%
    dplyr::filter((n > 0 & Segment.Level == level) | (Segment.Level < level & (Quant.Error < mad.threshold | n <= 3)))

  transposedCells <- innermostCells2 %>%
    select(all_of(train_colnames)) %>%
    t()

  cent_dist_df2 <- apply(data.frame(scaled_test_data), 1, find_path, transposedCells) %>%
    bind_rows()

  groupCols2 <- c(paste0("Segment.", c("Level", "Parent", "Child")), yVar)
  
  if ("Cell.ID" %in% names(innermostCells2)) { 
       groupCols2 <- c(groupCols2, "Cell.ID", "centroidRadius")
    } else{
      groupCols2 <- groupCols2
    }

  predict_test_data2 <-
    cbind(data.frame(scaled_test_data, "n" = 1), cent_dist_df2) %>%
    dplyr::left_join(
      innermostCells2 %>%
        select(all_of(groupCols2)) %>%
        cbind(Index = as.integer(row.names(.))),
      by = "Index"
    ) %>%
    select(-Index) %>%
    select(names(newdfMapping))

  # Considering margin of error
  predict_test_data3 <- predict_test_data2 %>% mutate(diff = centroidRadius - Quant.Error)
  predict_test_data3 <- predict_test_data3 %>% mutate(anomalyFlag = ifelse(Quant.Error < (mad.threshold), 0, 1))

  # Renaming fitted yVar with prefix Fitted
  if (!is.null(yVar)) {
    if (length(yVar) == 3) {
      for (i in 1:length(yVar)) {
        indexScored <- which(colnames(predict_test_data2) == yVar[i])
        colnames(predict_test_data2)[indexScored] <- paste0("Fitted.", yVar[i])
      }
    } else {
      indexScored <- which(colnames(predict_test_data2) == yVar)
      colnames(predict_test_data2)[indexScored] <- paste0("Fitted.", yVar)
    }
  }


  # Adding dep Variables back to scored data
  if (!is.null(yVar)) {
    predict_test_data2 <- merge(predict_test_data2, yVardf, by = 0) %>% select(-"Row.names")
    if (length(yVar) == 1) {
      indexFitted <- which(colnames(predict_test_data2) == "y")
      colnames(predict_test_data2)[indexFitted] <- paste0("Scored.", yVar)
    }
  }

  groupCols2 <- c(paste0("Segment.", c("Level", "Parent", "Child")), "anomalyFlag")
  # Calculating mean for scored data for each centroid
  groupCols3 <- c(paste0("Segment.", c("Level", "Parent", "Child")))
  # filter anomalous values for test dataset

  if (error_metric == "mean") {
    predictQE2 <- predict_test_data3 %>%
      group_by_at(groupCols2) %>%
      dplyr::summarise(
        n = sum(n),
        Quant.Error = mean(Quant.Error)
      )
    predictQE3 <- predict_test_data3 %>%
      group_by_at(groupCols3) %>%
      dplyr::summarise(
        n = sum(n),
        Quant.Error = mean(Quant.Error)
      )
  } else {
    predictQE2 <- predict_test_data3 %>% # with anamoly flags
      group_by_at(groupCols2) %>%
      dplyr::summarise(
        n = sum(n),
        Quant.Error = max(Quant.Error)
      )
    predictQE3 <- predict_test_data3 %>% # wo anamoly flags
      group_by_at(groupCols3) %>%
      dplyr::summarise(
        n = sum(n),
        Quant.Error = max(Quant.Error)
      )
  }

  # Calculating mean for scored data for each centroid
  groupCols2 <- c(paste0("Segment.", c("Level", "Parent", "Child")))
  newdfMapping <- newdfMapping %>% mutate(sumOriginal = Quant.Error * n)
  df_temp <- inner_join(predictQE2,
    newdfMapping %>% select(c(groupCols2, Quant.Error, sumOriginal, n)),
    by = groupCols2
  )
  df_temp2 <- inner_join(predictQE3,
    newdfMapping %>% select(c(groupCols2, Quant.Error, sumOriginal, n)),
    by = groupCols2
  )


  if (error_metric == "mean") {
    df_temp <- df_temp %>% mutate(Scored.Quant.Error = (sumOriginal + (Quant.Error.x * n.x)) / (n.x + n.y)) # sum original is wrong
    df_temp2 <- df_temp2 %>% mutate(Scored.Quant.Error = (sumOriginal + (Quant.Error.x * n.x)) / (n.x + n.y))
  } else {
    df_temp <- df_temp %>% mutate(Scored.Quant.Error = max(Quant.Error.x, Quant.Error.y))
    df_temp2 <- df_temp2 %>% mutate(Scored.Quant.Error = (sumOriginal + (Quant.Error.x * n.x)) / (n.x + n.y)) # should be maxScored
  }

  QECompareDf2 <- df_temp %>%
    mutate(
      Quant.Error.Diff = abs(Scored.Quant.Error - Quant.Error.y),
      `Quant.Error.Diff (%)` = abs(Scored.Quant.Error - Quant.Error.y) / Quant.Error.y * 100
    ) %>%
    dplyr::rename(Fitted.Quant.Error = Quant.Error.y, n = n.x) %>%
    select(-c("Quant.Error.x", "sumOriginal", "n.y"))

  plotList <- hvt.results.model[[2]] %>%
    unlist(., recursive = FALSE) %>%
    unlist(., recursive = FALSE)

  

  
  
  
  hvt_res1 <- hvt.results.model[[2]][[1]]$`1`
  hvt_res2 <- hvt.results.model[[3]]$summary$Cell.ID
  a <- 1: length(hvt_res1)
  b <- a[hvt_res2]
  b <-  as.vector(b)
  hvt_res2 <- stats::na.omit(b)
  
  coordinates_value1 <- lapply(1:length(hvt_res1), function(x) {
    centroids1 <- hvt_res1[[x]]
    coordinates1 <- centroids1$pt})
  cellID_coordinates <- do.call(rbind.data.frame, coordinates_value1)
  colnames(cellID_coordinates) <- c("x", "y")
  cellID_coordinates$Cell.ID <- hvt_res2
 # cellID_coordinates <- cellID_coordinates %>% arrange(Cell.ID)
  
  ##################
  boundaryCoords2 <-
    lapply(plotList, function(x) {
      data.frame(
        "Segment.Level" = x[["Segment.Level"]],
        "Segment.Parent" = x[["Segment.Parent"]],
        "Segment.Child" = x[["Segment.Child"]],
        "x" = x$pt["x"],
        "y" = x$pt["y"],
        "bp.x" = I(x$x),
        "bp.y" = I(x$y)
      )
    }) %>%
    bind_rows(.) 
  
  
  boundaryCoords2 <- merge(boundaryCoords2, cellID_coordinates, by = c("x" ,"y")) %>%
    right_join(.,
               QECompareDf2 %>% dplyr::filter(anomalyFlag == 1),
               by = paste0("Segment.", c("Level", "Parent", "Child"))
    )
  
  
  anomalyPlot <- plotHVT(
    hvt.results.model,
    child.level = child.level
  ) + ggtitle(paste(
    "Hierarchical Voronoi Tessellation for Level",
    child.level
  )) +
    theme(
      plot.title = element_text(
        size = 18,
        hjust = 0.5,
        margin = margin(0, 0, 20, 0)
      ),
      # legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  
  colour_scheme <- c(
    "#6E40AA", "#6B44B2", "#6849BA", "#644FC1", "#6054C8", "#5C5ACE", "#5761D3", "#5268D8", "#4C6EDB", "#4776DE", "#417DE0", "#3C84E1", "#368CE1",
    "#3194E0", "#2C9CDF", "#27A3DC", "#23ABD8", "#20B2D4", "#1DBACE", "#1BC1C9", "#1AC7C2", "#19CEBB", "#1AD4B3", "#1BD9AB", "#1DDFA3", "#21E39B",
    "#25E892", "#2AEB8A", "#30EF82", "#38F17B", "#40F373", "#49F56D", "#52F667", "#5DF662", "#67F75E", "#73F65A", "#7FF658", "#8BF457", "#97F357", "#A3F258"
  )
  
  if (nrow(boundaryCoords2) != 0) {
    hoverText <- paste(
      " Cell ID:",
      boundaryCoords2$Cell.ID,
      "<br>",
      "Segment.Level:",
      boundaryCoords2$Segment.Level,
      "<br>",
      "Segment.Parent:",
      boundaryCoords2$Segment.Parent,
      "<br>",
      "Segment.Child:",
      boundaryCoords2$Segment.Child,
      "<br>",
      "Number of observations:",
      boundaryCoords2$n,
      "<br>"
    )
  } else {
    hoverText <- NULL
  }
  # browser()
  anomalyPlot <- anomalyPlot + geom_polygon(
    data = boundaryCoords2,
    aes(
      x = bp.x,
      y = bp.y,
      group = interaction(Segment.Level, Segment.Parent, Segment.Child),
      fill = n,
      text = hoverText
    ),
    color = "red",
    size = 1
  ) +
    geom_point(data = boundaryCoords2 %>% distinct(x, y), aes(x = x, y = y), size = 1.5) +
    scale_fill_gradientn(colours = colour_scheme) +
    guides(colour = "none")
  
  plotlyPredict <- plotly::ggplotly(anomalyPlot, tooltip = "text")
  
  hoverText <- lapply(plotlyPredict$x$data, function(x) {
    if (!is.null(x$text)) {
      return(x$text)
    }
  }) %>% unlist()

  checkCell <- substr(hoverText, 1, 5) %in% " Cell"
  trace_vec <- seq_along(checkCell)[!checkCell]

  plotlyPredict <- plotlyPredict %>%
    plotly::layout(
      hoverlabel = list(bgcolor = "rgba(255,255,0,0.2)"),
      legend = list(
        title = list(text = "Level"),
        itemdoubleclick = FALSE,
        itemclick = "toggleothers",
        traceorder = "reversed"
      )
    ) %>%
    plotly::style(plotlyPredict, hoverinfo = "none", traces = trace_vec) %>%
    plotly::config(displayModeBar = FALSE)
  
  
  predict_test_data3 <- predict_test_data3 %>% mutate_if(is.numeric, round, digits = 4) 
  predict_test_data3 <- predict_test_data3 %>% dplyr::select(c("Segment.Level",  "Segment.Parent" ,"Segment.Child" , "n" ,            
                                                               "Cell.ID", "Quant.Error",  "centroidRadius" ,"diff" , "anomalyFlag", everything()))
  predict_test_dataRaw <- predict_test_data3
  predict_test_dataRaw[, train_colnames] <- data[, train_colnames]
################################################# Changes #################
  predicted_result <- hvt.results.model[[3]]$summary 
  current_predicted <- colnames(predicted_result)
  new_names <- paste0("pred_", current_predicted)
  colnames(predicted_result) <- new_names
  Cell.ID <- data.frame(predicted_result$pred_Cell.ID)
  predicted_result <- predicted_result %>% select(-c("pred_Segment.Level", "pred_Segment.Parent", "pred_Segment.Child", "pred_n", "pred_Cell.ID", "pred_Quant.Error", "pred_centroidRadius"))
  predicted_result <- cbind(Cell.ID, predicted_result)
  predicted_result <- data.table::setnames(predicted_result, "predicted_result.pred_Cell.ID", "Cell.ID")



  actuals <- predict_test_data3
  current_actual <- colnames(actuals)
  new_names <- paste0("act_", current_actual)
  colnames(actuals) <- new_names

  actuals$Row.No <- row.names(data)
  data_with_row <- data.frame(actuals$Row.No)
  data_with_cell <- data.frame(actuals$act_Cell.ID)
  actuals <- actuals %>% select(-c("act_Segment.Level", "act_Segment.Parent", "act_Segment.Child", "act_n", "act_Cell.ID", "act_Quant.Error", "act_centroidRadius", "act_diff", "act_anomalyFlag", "Row.No"))
  actuals_data <- cbind(data_with_row, actuals, data_with_cell)
  actuals_data <- data.table::setnames(actuals_data, "actuals.act_Cell.ID", "Cell.ID")
  actuals_data <- data.table::setnames(actuals_data, "actuals.Row.No", "Row.No")



  merged_df <- merge(actuals_data, predicted_result, by = "Cell.ID")
  merged_result <- merged_df %>% arrange(as.numeric(merged_df$Row.No))

  subtract_predicted_actual <- function(data, actual_prefix = "act_", predicted_prefix = "pred_") {
    actual_cols <- grep(paste0("^", actual_prefix), names(data), value = TRUE)
    df_new <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
    temp0 <<- data.frame(matrix(nrow = nrow(data)))
    for (col in actual_cols) {
      predicted_col <- gsub(actual_prefix, predicted_prefix, col)

      if (predicted_col %in% names(data)) {
        temp0[[predicted_col]] <<- abs(data[[col]] - data[[predicted_col]])
      }
    }
    temp0 <- temp0 %>% purrr::discard(~ all(is.na(.) | . == ""))
    df_new[, 1] <- rowMeans(temp0)
    return(df_new)
  }



  diff <- subtract_predicted_actual(merged_result)
  merged_result$diff <- diff$matrix.ncol...1..nrow...nrow.data..
  merged_result <- rename(merged_result, c("diff" = "diff"))


  desired_order <- c("Row.No", grep("^act_", colnames(merged_result), value = TRUE), "Cell.ID", grep("^pred_", colnames(merged_result), value = TRUE), "diff")
  df_reordered <- merged_result[, desired_order]
#################################################

  
if(analysis.plots) { 
  
  
    scored_data <- data.frame(Cell.ID = predict_test_data3$Cell.ID, names.column)
    scored_data <- scored_data[order(scored_data$Cell.ID), c("Cell.ID", "names.column")]
    reformed_data <- stats::aggregate(names.column ~ Cell.ID, scored_data, FUN = function(x) paste(x, collapse = ", "))
    
    boundaryCoords2_1 <-
      lapply(plotList, function(x) {
        data.frame(
          "Segment.Level" = x[["Segment.Level"]],
          "Segment.Parent" = x[["Segment.Parent"]],
          "Segment.Child" = x[["Segment.Child"]],
          "x" = x$pt["x"],
          "y" = x$pt["y"],
          "bp.x" = I(x$x),
          "bp.y" = I(x$y)
        )
      }) %>%
      bind_rows(.) 
    
    boundaryCoords2_1 <- merge(boundaryCoords2_1, cellID_coordinates, by = c("x" ,"y")) %>%
      right_join(.,
                 QECompareDf2,
                 by = paste0("Segment.", c("Level", "Parent", "Child"))
      )
    
    boundaryCoords2_1 <- merge(boundaryCoords2_1, reformed_data, by = "Cell.ID")
    
    
    boundaryCoords2_1<- boundaryCoords2_1 %>%
      group_by(Cell.ID) %>%
      mutate(
        sum_n = ifelse(any(anomalyFlag == 0) & any(anomalyFlag == 1), sum(n[anomalyFlag == 0][1], n[anomalyFlag == 1][1]), NA),
        n = ifelse(!is.na(sum_n), sum_n, n)
      ) %>%
      select(-sum_n)

    
    scoredPlot <- plotHVT(
      hvt.results.model,
      child.level = child.level
    ) +
      theme(
        plot.title = element_text(
          size = 18,
          hjust = 0.5,
          margin = margin(0, 0, 20, 0)
        ),
        # legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom"
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    
    
    colour_scheme <- c(
      "#6E40AA", "#6B44B2", "#6849BA", "#644FC1", "#6054C8", "#5C5ACE", "#5761D3", "#5268D8", "#4C6EDB", "#4776DE", "#417DE0", "#3C84E1", "#368CE1",
      "#3194E0", "#2C9CDF", "#27A3DC", "#23ABD8", "#20B2D4", "#1DBACE", "#1BC1C9", "#1AC7C2", "#19CEBB", "#1AD4B3", "#1BD9AB", "#1DDFA3", "#21E39B",
      "#25E892", "#2AEB8A", "#30EF82", "#38F17B", "#40F373", "#49F56D", "#52F667", "#5DF662", "#67F75E", "#73F65A", "#7FF658", "#8BF457", "#97F357", "#A3F258"
    )
    

    wrap_and_limit_text <- function(text, line_length = 100, max_chars = 500) {
      # First, limit the total text to max_chars
      if (nchar(text) > max_chars) {
        text <- substr(text, 1, max_chars - 3)
        text <- paste0(text, "...")
      }
      
      # Then wrap the text
      wrapped <- sapply(seq(1, nchar(text), line_length), function(i) {
        substr(text, i, min(i + line_length - 1, nchar(text)))
      })
      
      # Join the wrapped lines with <br> for HTML line breaks
      paste(wrapped, collapse = "<br>")
    }
    
    if (nrow(boundaryCoords2_1) != 0) {
      boundaryCoords2_1$hoverText <- apply(boundaryCoords2_1, 1, function(row) {
        full_names <- row["names.column"]
        formatted_names <- wrap_and_limit_text(full_names, line_length = 100, max_chars = 500)
        
        paste(
          "Cell ID:", row["Cell.ID"],
          "<br>",
          "Number of observations:", row["n"],
          "<br>",
          "Name of observations:<br>", formatted_names
        )
      })
    } else {
      boundaryCoords2_1$hoverText <- NULL
    }
    
    # Create the plot
    scoredPlot <- scoredPlot +
      geom_polygon(
        data = boundaryCoords2_1,
        aes(
          x = bp.x,
          y = bp.y,
          group = interaction(Segment.Level, Segment.Parent, Segment.Child),
          fill = n,
          text = hoverText
        ),
        color = "black",
        size = 0.5
      ) +
      scale_fill_gradientn(colours = colour_scheme) +
      guides(colour = "none") +
      geom_point(
        data = boundaryCoords2_1 %>% distinct(x, .keep_all = TRUE),
        aes(x = x, y = y, text = hoverText),
        size = 1
      )
    
   
    plotlyscored <- plotly::ggplotly(scoredPlot, tooltip = "text")
    plotlyscored <- plotlyscored %>%
      plotly::layout(
        title = list(
          text = "Scored Heatmap with Cell-Level Observations",
          x = 0, 
          y = 0.99,
          xanchor = "left", font = list(size = 20)),
        hovermode = 'closest',
        hoverdistance = 100,
        hoverlabel = list(
          bgcolor = "rgba(255,255,0,0.2)",
          font = list(size = 10),
          align = "left",
          namelength = -1),
        legend = list(
          title = list(text = "Level"),
          itemdoubleclick = FALSE,
          itemclick = "toggleothers",
          traceorder = "reversed"
        )
      ) %>%
      plotly::config(displayModeBar = TRUE)

#############################  
    names_data <- boundaryCoords2_1 %>% dplyr::select("Cell.ID","names.column")    
    names_data <- names_data %>% 
      distinct(Cell.ID, .keep_all = TRUE) 
    a <- cellID_coordinates %>% arrange(Cell.ID)
    centroid_data = merge(a, names_data, by = 'Cell.ID') %>% as.data.frame()
    

  
  states_data <- boundaryCoords2_1 %>% dplyr::select("Cell.ID","names.column")
  states_data <- states_data %>% 
    distinct(Cell.ID, .keep_all = TRUE)
}
 #################################################
  #MODEL INFO Rewriting
  input_dataset <- hvt.results.model[["model_info"]][["input_parameters"]][["input_dataset"]]
  n_cells <-hvt.results.model[["model_info"]][["input_parameters"]][["n_cells"]]
  compression_percentage <- compression_percentage <- paste0(hvt.results.model[[3]][["compression_summary"]][["percentOfCellsBelowQuantizationErrorThreshold"]] * 100, "%")
  quantization_error <-  hvt.results.model[["model_info"]][["input_parameters"]][["quant.err"]]
  
  trained_model <- list(
    input_dataset = input_dataset, 
    no_of_cells = n_cells, 
    compression_percentage = compression_percentage,
    quantization_error = quantization_error)
  
  score_dataset <- paste0(data_structure_score[1] ," Rows & ", data_structure_score[2], " Columns")
  qe_range_min <-min(predict_test_data3$Quant.Error)
  qe_range_max <- max(predict_test_data3$Quant.Error)
  qe_range <- paste0( qe_range_min, " to " ,qe_range_max)
  no_of_anomaly_cells <- sum(QECompareDf2$anomalyFlag == 1)
  no_of_anomaly_data <- sum(QECompareDf2$n[QECompareDf2$anomalyFlag == 1])
  
  scored_model <- list(
    input_dataset = score_dataset,
    scored_qe_range = qe_range,
    mad.threshold = mad.threshold,
    no_of_anomaly_datapoints = no_of_anomaly_data,
    no_of_anomaly_cells = no_of_anomaly_cells
  )
  
  
  ####################################
  if (analysis.plots){
  prediction_list <- list(
    scoredPredictedData = predict_test_data3,
    actual_predictedTable = df_reordered,
    QECompareDf = QECompareDf2,
    anomalyPlot = plotlyPredict,
    scoredPlotly = plotlyscored,
  scoredPlot_test= scoredPlot,
  cellID_coordinates =cellID_coordinates,
    centroidData = centroid_data,
    states_data = states_data,
    predictInput = c("depth" = child.level, "quant.err" = mad.threshold),
    model_mad_plots = list(),
    model_info = list(type = "hvt_prediction",
                      trained_model_summary = trained_model,
                      scored_model_summary =scored_model)
  
  )}else{
  
  prediction_list <- list(
    scoredPredictedData = predict_test_data3,
    actual_predictedTable = df_reordered,
    QECompareDf = QECompareDf2,
    anomalyPlot = plotlyPredict,
    cellID_coordinates =cellID_coordinates,
    predictInput = c("depth" = child.level, "quant.err" = mad.threshold),
    model_mad_plots = list(),
    model_info = list(type = "hvt_prediction", 
                      trained_model_summary = trained_model,
                      scored_model_summary =scored_model)
  )}
  
  model_mad_plots <- NA
 
  if (!all(is.na(hvt.results.model[[4]]))) {
    mtrain <- hvt.results.model[[4]]$mad_plot_train + ggtitle("Mean Absolute Deviation Plot: Calibration on Train Data")
  }
  if (!all(is.na(hvt.results.model[[5]]))) {
    mtest <- hvt.results.model[[5]][["mad_plot"]] + ggtitle("Mean Absolute Deviation Plot:Validation")
  }

  if (!all(is.na(hvt.results.model[[4]])) & !all(is.na(hvt.results.model[[5]]))) {
    model_mad_plots <- list(mtrain = mtrain, mtest = mtest)
  }


  prediction_list$model_mad_plots <- model_mad_plots
  
  class(prediction_list) <- "hvt.object"
  return(prediction_list)
})
  
}

