#' @name plotHVT 
#' @title Plot the hierarchical tessellations.
#' @description This is the main plotting function to construct hierarchical voronoi tessellations in 1D,2D or
#' Interactive surface plot.
#' @param hvt.results (1D/2DProj/2Dhvt/2Dheatmap/surface_plot) List. A list containing the output of \code{trainHVT} function
#' which has the details of the tessellations to be plotted.
#' @param plot.type Character. An option to indicate which type of plot should be generated. Accepted entries are 
#' '1D','2Dproj','2Dhvt','2Dheatmap'and 'surface_plot'. Default value is '2Dhvt'.
#' @param line.width (2Dhvt/2Dheatmap) Numeric Vector. A vector indicating the line widths of the
#' tessellation boundaries for each level.
#' @param color.vec (2Dhvt/2Dheatmap) Vector. A vector indicating the colors of the boundaries of
#' the tessellations at each level.
#' @param centroid.size (2Dhvt/2Dheatmap) Numeric Vector. A vector indicating the size of centroids
#' for each level.
#' @param centroid.color (2Dhvt/2Dheatmap) Numeric Vector. A vector indicating the color of centroids
#' for each level.
#' @param child.level (2Dheatmap/surface_plot) Numeric. Indicating the level for which the plot should
#' be displayed
#' @param hmap.cols (2Dheatmap/surface_plot) Numeric or Character. The column number or column name from
#' the dataset indicating the variables for which the heat map is to be plotted.
#' @param quant.error.hmap (2Dheatmap) Numeric. A number representing the quantization error threshold to be highlighted in the heatmap. 
#' When a value is provided, it will emphasize cells with quantization errors equal or less than the specified threshold,
#' indicating that these cells cannot be further subdivided in the next depth layer. The default value is NULL,
#' meaning all cells will be colored in the heatmap across various depths.
#' @param separation_width (surface_plot) Numeric. An integer indicating the width between hierarchical levels in surface plot
#' @param layer_opacity (surface_plot) Numeric. A vector indicating the opacity of each hierarchical levels in surface plot
#' @param dim_size  (surface_plot) Numeric. An integer controls the resolution or granularity of the 3D surface grid
#' @param cell_id (2Dhvt/2Dheatmap) Logical. A logical indicating whether the cell IDs should be displayed
#' @param cell_id_position (2Dhvt/2Dheatmap) Character. A character indicating the position of the cell IDs. Accepted entries are 'top' , 
#' 'bottom', 'left' and 'right'.
#' @param cell_id_size (2Dhvt/2Dheatmap) Numeric. A numeric vector indicating the size of the cell IDs for all levels.
#' @returns plot object containing the visualizations of reduced dimension(1D/2D) for the given dataset.
#' @author Shubhra Prakash <shubhra.prakash@@mu-sigma.com>, Sangeet Moy Das <sangeet.das@@mu-sigma.com>, Vishwavani <vishwavani@@mu-sigma.com>
#' @seealso \code{\link{trainHVT}} 
#' @keywords Tessellation_and_Heatmap
#' @importFrom magrittr %>%
#' @import ggplot2
#' @examples
#' data("EuStockMarkets")
#' hvt.results <- trainHVT(EuStockMarkets, n_cells = 60, depth = 1, quant.err = 0.1, 
#'                        distance_metric = "L1_Norm", error_metric = "max",
#'                        normalize = TRUE,quant_method="kmeans")
#'                        
#' #change the 'plot.type' argument to '2Dproj' or '2DHVT' to visualize respective plots.                      
#' plotHVT(hvt.results, plot.type='1D')
#' 
#' #change the 'plot.type' argument to 'surface_plot' to visualize the Interactive surface plot                   
#' plotHVT(hvt.results,child.level = 1, 
#' hmap.cols = "DAX", plot.type = '2Dheatmap')
#' @export plotHVT


plotHVT <- function(hvt.results, 
                    line.width = 0.5, 
                    color.vec =  'black', 
                    centroid.size = 0.6, 
                    centroid.color = "black", 
                    child.level = 1, 
                    hmap.cols,
                    separation_width = 7, 
                    layer_opacity = c(0.5, 0.75, 0.99), 
                    dim_size = 1000, 
                    plot.type = '2Dhvt',
                    quant.error.hmap = NULL,
                    cell_id = FALSE,
                    cell_id_position="bottom", 
                    cell_id_size = 2.6) {
  
  lev <- NULL
  n_cells <- max(stats::na.omit(hvt.results[[3]][["summary"]][["Cell.ID"]]))
  n_cells.hmap <- hvt.results[["model_info"]][["input_parameters"]][["n_cells"]]
  pch1 = 21
    
  if (is.null(plot.type)) {
    plot.type <- '2Dhvt'
  }
  
  if (plot.type == '1D') {
    generic_col=c("Segment.Level","Segment.Parent","Segment.Child","n","Quant.Error")
    hvq_k <- hvt.results[[6]]  
    temp_summary=hvq_k[["summary"]] %>% dplyr::select(!generic_col) %>% dplyr::mutate(id=row_number())
    cent_val= temp_summary %>% subset(.,stats::complete.cases(.)) 
    set.seed(123)
    sammon_1d_cord <- MASS::sammon(
      d = stats::dist(cent_val %>% dplyr::select(!id),method = "manhattan"),
      niter = 10 ^ 5,
      trace = FALSE,
      k=1
    )$points
    temp_df=data.frame(sammon_1d_cord,id=cent_val$id)%>%dplyr::arrange(sammon_1d_cord) %>% dplyr::mutate(Cell.ID=row_number()) %>% dplyr::select(!sammon_1d_cord)
    temp_summary = dplyr::left_join(temp_summary,temp_df,by="id") %>% select(!"id")
    hvq_k[["summary"]]$Cell.ID=temp_summary$Cell.ID
    
    x <-sammon_1d_cord
    y <- hvq_k[["summary"]][["Cell.ID"]]
    data_plot <- data.frame(x,y)
    
      if(length(y) <= 100) {dot_size <- 5} else if(length(y) <= 500) {dot_size <- 2.5} 
      else if(length(y) <= 1000) {dot_size <- 1.5} else {dot_size <- 1}
    
    plotly_obj <-  plotly::plot_ly(data_plot, x = ~y, y = ~x, type = 'scatter', mode = 'markers',
            marker = list(size = dot_size, color = 'blue', symbol = 'circle'),
            text = ~paste('Cell ID:', y, '<br>1D point:', round(x,4)), 
            hoverinfo = 'text') %>% 
            plotly::layout(title = 'Sammons 1D x Cell ID',
             xaxis = list(title = 'Cell ID', zeroline = FALSE),
             yaxis = list(title = '1D points',zeroline = FALSE))
    
     return(suppressMessages(plotly_obj))  
      
  } else if (plot.type == '2Dproj'){
   
     hvt_centroids_list <- hvt.results
      
    hvt_coordinates<- hvt_centroids_list[[2]][[1]][["1"]]
    centroids <- list()
    coordinates_value <- lapply(1:length(hvt_coordinates), function(x){
      centroids <-hvt_coordinates[[x]]
      coordinates <- centroids$pt
    })
    centroid_coordinates<<- do.call(rbind.data.frame, coordinates_value)  
    colnames(centroid_coordinates) <- c("x_coord","y_coord")
    centroid_coordinates$Row.No <- as.numeric(row.names(centroid_coordinates)) 
    centroid_coordinates <- centroid_coordinates %>% dplyr::select(Row.No,x_coord,y_coord)
    centroid_coordinates1 <- centroid_coordinates %>% data.frame() %>% round(4)

    gg_proj <- ggplot(centroid_coordinates1, aes(x_coord, y_coord)) +
      geom_point(color = "blue") +
      labs(title = "2D Projection plot of centroids",
           x = "X", y = "Y")
    
    
    return(suppressMessages(gg_proj))
    
     } else if (plot.type == '2Dhvt') {
      
    hvt_list <- hvt.results
    
    child.level <- min(child.level, max(hvt_list[[3]][["summary"]] %>% stats::na.omit() %>% dplyr::select("Segment.Level")))
    
    
    min_x <- 1e9
    min_y <- 1e9
    max_x <- -1e9
    max_y <- -1e9
    depthVal <- c()
    clusterVal <- c()
    childVal <- c()
    value <- c()
    x_pos <- c()
    y_pos <- c()
    x_cor <- c()
    y_cor <- c()
    depthPos <- c()
    clusterPos <- c()
    childPos <- c()
    levelCluster <- c()
    for (clusterNo in 1:length(hvt_list[[2]][[1]][[1]])) {
      bp_x <- hvt_list[[2]][[1]][[1]][[clusterNo]][["x"]]
      bp_y <- hvt_list[[2]][[1]][[1]][[clusterNo]][["y"]]
      
      
      if (min(bp_x) < min_x) {
        min_x <- min(bp_x)
      }
      if (max(bp_x) > max_x) {
        max_x <- max(bp_x)
      }
      if (min(bp_y) < min_y) {
        min_y <- min(bp_y)
      }
      if (max(bp_y) > max_y) {
        max_y <- max(bp_y)
      }
    }
    
    
    for (depth in 1:child.level) {
      for (clusterNo in 1:length(hvt_list[[2]][[depth]])) {
        for (childNo in 1:length(hvt_list[[2]][[depth]][[clusterNo]])) {
          current_cluster <- hvt_list[[2]][[depth]][[clusterNo]][[childNo]]
          
          x <- as.numeric(current_cluster[["x"]])
          y <- as.numeric(current_cluster[["y"]])
          x_cor <- c(x_cor, as.numeric(current_cluster[["pt"]][["x"]]))
          y_cor <- c(y_cor, as.numeric(current_cluster[["pt"]][["y"]]))
          depthVal <- c(depthVal, depth)
          clusterVal <- c(clusterVal, clusterNo)
          childVal <- c(childVal, childNo)
          depthPos <- c(depthPos, rep(depth, length(x)))
          clusterPos <- c(clusterPos, rep(clusterNo, length(x)))
          childPos <- c(childPos, rep(childNo, length(x)))
          x_pos <- c(x_pos, x)
          y_pos <- c(y_pos, y)
          levelCluster <- c(levelCluster, depth)
        }
      }
    }
    
    valuesDataframe <- data.frame(
      depth = depthVal,
      cluster = clusterVal,
      child = childVal
    )
    
    
    positionsDataframe <- data.frame(
      depth = depthPos,
      cluster = clusterPos,
      child = childPos,
      x = x_pos,
      y = y_pos
    )
    
    
    centroidDataframe <-
      data.frame(x = x_cor, y = y_cor, lev = levelCluster)
    
    datapoly <-
      merge(valuesDataframe,
            positionsDataframe,
            by = c("depth", "cluster", "child")
      )

      hvt_res1 <- hvt_list[[2]][[1]]$`1`
      hvt_res2 <- hvt_list[[3]]$summary$Cell.ID
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
      centroidDataframe_2 <- merge(cellID_coordinates, centroidDataframe, by = c("x" ,"y"))
      
    p <- ggplot2::ggplot()  
    for (i in child.level:1) {
      p <-
        p + ggplot2::geom_polygon(
          data = datapoly[which(datapoly$depth == i), ],
          ggplot2::aes(
            x = x,
            y = y,
            color = factor(depth),
            size = factor(depth),
            group = interaction(depth, cluster, child),
          ),
          fill = NA
        ) +
        ggplot2::scale_colour_manual(values = color.vec) +
        ggplot2::scale_size_manual(values = line.width, guide = "none") +
        ggplot2::labs(color = "Level")
    }
    
    if (cell_id == TRUE) {

    for (depth in 1:child.level) {
      depth_size <- centroid.size[child.level - depth + 1]
      centroid_color <- centroid.color[child.level - depth + 1]
      
      p <- p + ggplot2::geom_point(
        data = centroidDataframe[centroidDataframe["lev"] == depth, ],
        ggplot2::aes(x = x, y = y),
        size = depth_size,
        pch = pch1,
        fill = centroid_color,
        color = centroid_color
      ) 
    }

    subset_data <- subset(centroidDataframe_2, lev == 1)
    
    
    alignments <- list(
      right = list(hjust = -0.5, vjust = 0.5),
      left = list(hjust = 1.5, vjust = 0.5),
      bottom = list(hjust = 0.5, vjust = 1.5),
      top = list(hjust = 0.5, vjust = -0.5) # Default case
    )
    
    # Get alignment values for the current position
    alignment <- rlang::`%||%`(alignments[[cell_id_position]], alignments$bottom)
    
    # Add geom_text to the plot
    p <- p + ggplot2::geom_text(
      data = subset_data,
      ggplot2::aes(x = x, y = y, label = Cell.ID),
      size = cell_id_size,
      color = "black",
      hjust = alignment$hjust,
      vjust = alignment$vjust,
      nudge_x = -0.02 * nchar(as.character(subset_data$Cell.ID)))
     
    }else {
       for (depth in 1:child.level) {
         depth_size <- centroid.size[child.level - depth + 1]
         centroid_color <- centroid.color[child.level - depth + 1]
         
         p <- p + ggplot2::geom_point(
           data = centroidDataframe[centroidDataframe["lev"] == depth, ],
           ggplot2::aes(x = x, y = y),
           size = depth_size,
           pch = pch1,
           fill = centroid_color,
           color = centroid_color
         ) 
       }
       
     }
  
 

    p <- p +
      ggplot2::scale_color_manual(
        name = "Level",
        values = color.vec
      ) +
      ggplot2::theme_bw() + ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        plot.title = element_text(
          size = 20,
          hjust = 0.5,
          margin = margin(0, 0, 20, 0)
        ),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank()
      ) + ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::geom_label(
        label = centroidDataframe$outlier_cell,
        nudge_x = 0.45, nudge_y = 0.1,
        check_overlap = TRUE,
        label.padding = unit(0.55, "lines"),
        label.size = 0.4,
        color = "white",
        fill = "#038225"
      ) +
      ggplot2::ggtitle(paste0("Tesellation Plot of ", n_cells, " cells"))+ theme(
        plot.title = element_text(hjust = 0, size = 10 ))
   

    return(suppressMessages(p))
    
    
  } else if (plot.type == '2Dheatmap') {
    hvt_list <- hvt.results
    
    child.level <- min(child.level, max(hvt_list[[3]][["summary"]] %>% stats::na.omit() %>% dplyr::select("Segment.Level")))
    summaryDF <- hvt_list[[3]][["summary"]]
    valuesDataframe <- data.frame(
      depth = 0,
      cluster = 0,
      child = 0,
      qe = 0,
      n_cluster = 0,
      value = 0,
      cellID = 0
    )
    positionsDataframe <- data.frame(
      depth = 0,
      cluster = 0,
      child = 0,
      x = 0,
      y = 0
    )
    
    centroidDataframe <- data.frame(x = 0, y = 0, lev = 0)
    
    for (depth in 1:child.level) {
      if (depth < 3) {
        for (clusterNo in names(hvt_list[[2]][[depth]])) {
          for (childNo in 1:length(hvt_list[[2]][[depth]][[clusterNo]])) {
            if (!is.null(hvt_list[[2]][[depth]][[clusterNo]])) {
              summaryFilteredDF <-
                summaryDF %>% dplyr::filter(
                  Segment.Level == depth,
                  Segment.Parent == clusterNo,
                  Segment.Child == childNo
                )
              
              current_cluster <- hvt_list[[2]][[depth]][[clusterNo]][[childNo]]
              
              
              val <- summaryFilteredDF[, hmap.cols]
              qe <- summaryFilteredDF[, "Quant.Error"]
              ncluster <- summaryFilteredDF[, "n"]
              x <- as.numeric(current_cluster[["x"]])
              y <- as.numeric(current_cluster[["y"]])
              if ("Cell.ID" %in% colnames(summaryFilteredDF)) {
                cell_ID <- summaryFilteredDF[, "Cell.ID"]
              } else {
                cell_ID <- 0
              }
              valuesDataframe <-
                rbind(
                  valuesDataframe,
                  data.frame(
                    depth = depth,
                    cluster = clusterNo,
                    child = childNo,
                    qe = qe,
                    n_cluster = ncluster,
                    value = val,
                    cellID = cell_ID
                  )
                )
              positionsDataframe <-
                rbind(
                  positionsDataframe,
                  data.frame(
                    depth = rep(depth, length(x)),
                    cluster = rep(clusterNo, length(x)),
                    child = rep(childNo, length(x)),
                    x = x,
                    y = y
                  )
                )
              
              centroidDataframe <-
                rbind(
                  centroidDataframe,
                  data.frame(
                    x = as.numeric(current_cluster[["pt"]][["x"]]),
                    y = as.numeric(current_cluster[["pt"]][["y"]]),
                    lev = depth
                  )
                )
            }
          }
        }
      } else {
        for (clusterNo in 1:n_cells.hmap^(child.level - 1)) {
          for (childNo in 1:length(hvt_list[[2]][[depth]][[as.character(clusterNo)]])) {
            if (!is.null(hvt_list[[2]][[depth]][[as.character(clusterNo)]])) {
              summaryFilteredDF <-
                summaryDF %>% dplyr::filter(
                  Segment.Level == depth,
                  Segment.Parent == clusterNo,
                  Segment.Child == childNo
                )
              
              current_cluster <- hvt_list[[2]][[depth]][[as.character(clusterNo)]][[childNo]]
              
              
              val <- summaryFilteredDF[, hmap.cols]
              qe <- summaryFilteredDF[, "Quant.Error"]
              
              x <- as.numeric(current_cluster[["x"]])
              y <- as.numeric(current_cluster[["y"]])
              if ("Cell.ID" %in% colnames(summaryFilteredDF)) {
                cell_ID <- summaryFilteredDF[, "Cell.ID"]
              } else {
                cell_ID <- 0
              }
              valuesDataframe <-
                rbind(
                  valuesDataframe,
                  data.frame(
                    depth = depth,
                    cluster = clusterNo,
                    child = childNo,
                    qe = qe,
                    n_cluster = ncluster,
                    value = val,
                    cellID = cell_ID
                  )
                )
              positionsDataframe <-
                rbind(
                  positionsDataframe,
                  data.frame(
                    depth = rep(depth, length(x)),
                    cluster = rep(clusterNo, length(x)),
                    child = rep(childNo, length(x)),
                    x = x,
                    y = y
                  )
                )
              
              centroidDataframe <-
                rbind(
                  centroidDataframe,
                  data.frame(
                    x = as.numeric(current_cluster[["pt"]][["x"]]),
                    y = as.numeric(current_cluster[["pt"]][["y"]]),
                    lev = depth
                  )
                )
            }
          }
        }
      }
    }
    valuesDataframe <- valuesDataframe[2:nrow(valuesDataframe), ]
    
    positionsDataframe <-
      positionsDataframe[2:nrow(positionsDataframe), ]
    
    centroidDataframe <-
      centroidDataframe[2:nrow(centroidDataframe), ]
    
    
    datapoly <-
      merge(valuesDataframe,
            positionsDataframe,
            by = c("depth", "cluster", "child")
      )

    p <- ggplot2::ggplot()
    colour_scheme <- c(
      "#6E40AA", "#6B44B2", "#6849BA", "#644FC1", "#6054C8", "#5C5ACE", "#5761D3", "#5268D8", "#4C6EDB", "#4776DE", "#417DE0", "#3C84E1", "#368CE1",
      "#3194E0", "#2C9CDF", "#27A3DC", "#23ABD8", "#20B2D4", "#1DBACE", "#1BC1C9", "#1AC7C2", "#19CEBB", "#1AD4B3", "#1BD9AB", "#1DDFA3", "#21E39B",
      "#25E892", "#2AEB8A", "#30EF82", "#38F17B", "#40F373", "#49F56D", "#52F667", "#5DF662", "#67F75E", "#73F65A", "#7FF658", "#8BF457", "#97F357", "#A3F258"
    )
    data <- datapoly
    if (child.level > 1) {
      for (i in 1:(child.level - 1)) {
        #index_tess <- which(data$depth == i & data$qe > quant.error.hmap & data$n_cluster > 3)
         index_tess <- which( data$qe > quant.error.hmap )
        
        if (length(index_tess) >0){
          data <- data[-index_tess, ] 
        } else {
            data = data
          }
         
        #rm(index_tess)
      }
    }
    
    hvt_res1 <- hvt_list[[2]][[1]]$`1`
    hvt_res2 <- hvt_list[[3]]$summary$Cell.ID
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
    centroidDataframe_2 <- merge(cellID_coordinates, centroidDataframe, by = c("x" ,"y"))
    
    p <- ggplot2::ggplot()
 
   
    
    # changing hoverText for torus demo
    if ("Cell.ID" %in% colnames(summaryFilteredDF)) {
      p <-
        p + ggplot2::geom_polygon(
          data = data,
          ggplot2::aes(
            x = x,
            y = y,
            group = interaction(depth, cluster, child),
            fill = value,
            text = paste(
              " Cell.ID:", cellID,
              "<br>", "Segment.Level:", depth,
              "<br>", "Segment.Parent:", cluster,
              "<br>", "Segment.Child:", child
            )
          )
        ) +
        ggplot2::scale_fill_gradientn(colours = colour_scheme, guide = guide_colorbar(order = 2)) +
        ggplot2::labs(fill = hmap.cols)
    } else {
      p <-
        p + ggplot2::geom_polygon(
          data = data,
          ggplot2::aes(
            x = x,
            y = y,
            group = interaction(depth, cluster, child),
            fill = value
          )
        ) +
        ggplot2::scale_fill_gradientn(colours = colour_scheme,guide = guide_colorbar(order = 2)) +
        ggplot2::labs(fill = hmap.cols)
    }
    
    
    for (i in child.level:1) {
      p <-
        p + ggplot2::geom_polygon(
          data = datapoly[which(datapoly$depth == i), ],
          ggplot2::aes(
            x = x,
            y = y,
            color = factor(depth),
            size = factor(depth),
            group = interaction(depth, cluster, child)
          ),
          fill = NA
        ) +
        ggplot2::scale_colour_manual(values = color.vec,guide = guide_legend(order = 1)) +
        ggplot2::scale_size_manual(values = line.width, guide = FALSE) +
        ggplot2::labs(color = "Level")
    }
    
    if (cell_id == TRUE) {

      for (depth in 1:child.level) {
        depth_size <- centroid.size[child.level - depth + 1]
        centroid_color <- centroid.color[child.level - depth + 1]
        
        p <- p + ggplot2::geom_point(
          data = centroidDataframe[centroidDataframe["lev"] == depth, ],
          ggplot2::aes(x = x, y = y),
          size = depth_size,
          fill = centroid_color,
          color = centroid_color
        ) 
      }
      
      subset_data <- subset(centroidDataframe_2, lev == 1)  
      
      alignments <- list(
        right = list(hjust = -0.5, vjust = 0.5),
        left = list(hjust = 1.5, vjust = 0.5),
        bottom = list(hjust = 0.5, vjust = 1.5),
        top = list(hjust = 0.5, vjust = -0.5) )
      
      alignment <- rlang::`%||%`(alignments[[cell_id_position]], alignments$bottom)
      
      p <- p + ggplot2::geom_text(
        data = subset_data,
        ggplot2::aes(x = x, y = y, label = Cell.ID),
        size = cell_id_size,
        color = "black",
        hjust = alignment$hjust,
        vjust = alignment$vjust,
        nudge_x = -0.02 * nchar(as.character(subset_data$Cell.ID)))
    }else {
      for (depth in 1:child.level) {
        depth_size <- centroid.size[child.level - depth + 1]
        centroid_color <- centroid.color[child.level - depth + 1]
        
        p <- p + ggplot2::geom_point(
          data = centroidDataframe[centroidDataframe["lev"] == depth, ],
          ggplot2::aes(x = x, y = y),
          size = depth_size,
          pch = pch1,
          fill = centroid_color,
          color = centroid_color
        ) 
      }
      
    }
    
      

    
    # Continue with the rest of the plot modifications
    p <- p + ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      plot.title = element_text(
        size = 20,
        hjust = 0.5,
        margin = margin(0, 0, 20, 0)
      ),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_blank()
    ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::ggtitle(paste0("Heatmap of ", n_cells, " cells"))+ theme(
        plot.title = element_text(hjust = 0, size = 10 ))
   
    
   
    
    return(suppressMessages(p))
    
  } else if (plot.type == 'surface_plot') {
    
    child.level <- min(child.level, max(hvt.results[[3]][["summary"]]
                                     %>% stats::na.omit()
                                     %>% dplyr::select("Segment.Level")))
    summaryDF <- hvt.results[[3]][["summary"]]
    valuesDataframe <- data.frame(
      depth = 0,
      cluster = 0,
      child = 0,
      qe = 0,
      n_cluster = 0,
      value = 0,
      cellID = 0
    )
    positionsDataframe <- data.frame(
      depth = 0,
      cluster = 0,
      child = 0,
      x = 0,
      y = 0
    )
    
    centroidDataframe <- data.frame(x = 0, y = 0, lev = 0)
    
    for (depth in 1:child.level) {
      if (depth < 3) {
        for (clusterNo in names(hvt.results[[2]][[depth]])) {
          for (childNo in 1:length(hvt.results[[2]][[depth]][[clusterNo]])) {
            if (!is.null(hvt.results[[2]][[depth]][[clusterNo]])) {
              summaryFilteredDF <-
                summaryDF %>% dplyr::filter(
                  Segment.Level == depth,
                  Segment.Parent == clusterNo,
                  Segment.Child == childNo
                )
              
              current_cluster <- hvt.results[[2]][[depth]][[clusterNo]][[childNo]]
              
              
              val <- summaryFilteredDF[, hmap.cols]
              qe <- summaryFilteredDF[, "Quant.Error"]
              ncluster <- summaryFilteredDF[, "n"]
              x <- as.numeric(current_cluster[["x"]])
              y <- as.numeric(current_cluster[["y"]])
              if ("Cell.ID" %in% colnames(summaryFilteredDF)) {
                cell_ID <- summaryFilteredDF[, "Cell.ID"]
              } else {
                cell_ID <- 0
              }
              valuesDataframe <-
                rbind(
                  valuesDataframe,
                  data.frame(
                    depth = depth,
                    cluster = clusterNo,
                    child = childNo,
                    qe = qe,
                    n_cluster = ncluster,
                    value = val,
                    cellID = cell_ID
                  )
                )
              positionsDataframe <-
                rbind(
                  positionsDataframe,
                  data.frame(
                    depth = rep(depth, length(x)),
                    cluster = rep(clusterNo, length(x)),
                    child = rep(childNo, length(x)),
                    x = x,
                    y = y
                  )
                )
              
              centroidDataframe <-
                rbind(
                  centroidDataframe,
                  data.frame(
                    x = as.numeric(current_cluster[["pt"]][["x"]]),
                    y = as.numeric(current_cluster[["pt"]][["y"]]),
                    lev = depth
                  )
                )
            }
          }
        }
      } else {
        for (clusterNo in 1:n_cells.hmap^(child.level - 1)) {
          for (childNo in 1:length(hvt.results[[2]][[depth]][[as.character(clusterNo)]])) {
            if (!is.null(hvt.results[[2]][[depth]][[as.character(clusterNo)]])) {
              summaryFilteredDF <-
                summaryDF %>% dplyr::filter(
                  Segment.Level == depth,
                  Segment.Parent == clusterNo,
                  Segment.Child == childNo
                )
              
              current_cluster <- hvt.results[[2]][[depth]][[as.character(clusterNo)]][[childNo]]
              
              
              val <- summaryFilteredDF[, hmap.cols]
              qe <- summaryFilteredDF[, "Quant.Error"]
              
              x <- as.numeric(current_cluster[["x"]])
              y <- as.numeric(current_cluster[["y"]])
              if ("Cell.ID" %in% colnames(summaryFilteredDF)) {
                cell_ID <- summaryFilteredDF[, "Cell.ID"]
              } else {
                cell_ID <- 0
              }
              valuesDataframe <-
                rbind(
                  valuesDataframe,
                  data.frame(
                    depth = depth,
                    cluster = clusterNo,
                    child = childNo,
                    qe = qe,
                    n_cluster = ncluster,
                    value = val,
                    cellID = cell_ID
                  )
                )
              positionsDataframe <-
                rbind(
                  positionsDataframe,
                  data.frame(
                    depth = rep(depth, length(x)),
                    cluster = rep(clusterNo, length(x)),
                    child = rep(childNo, length(x)),
                    x = x,
                    y = y
                  )
                )
              
              centroidDataframe <-
                rbind(
                  centroidDataframe,
                  data.frame(
                    x = as.numeric(current_cluster[["pt"]][["x"]]),
                    y = as.numeric(current_cluster[["pt"]][["y"]]),
                    lev = depth
                  )
                )
            }
          }
        }
      }
    }
    valuesDataframe <- valuesDataframe[2:nrow(valuesDataframe), ]
    
    positionsDataframe <-
      positionsDataframe[2:nrow(positionsDataframe), ]
    
    centroidDataframe <-
      centroidDataframe[2:nrow(centroidDataframe), ]
    
    
    datapoly <-
      merge(valuesDataframe,
            positionsDataframe,
            by = c("depth", "cluster", "child")
      )
    
    
    
    Level_list <- lapply(1:child.level, function(x) {
      temp_df <- datapoly[datapoly$depth == x, ]
      min_x <- min(datapoly$x)
      max_x <- max(datapoly$x)
      min_y <- min(datapoly$y)
      max_y <- max(datapoly$y)
      temp_df <- temp_df %>% mutate(name = paste0("Segment", cluster, "-", child))
      # segments_number=unique()
      returnList <- split(temp_df, f = temp_df$name)
      return(list(
        returnList,
        min_x,
        max_x,
        min_y,
        max_y
      ))
    })
    
    
    
    #####################################################################################
    depth_wise_surface <- function(finalList
                                  
    ) {
      column <- hmap.cols
      min_x <- unlist(finalList[[2]])
      max_x <- unlist(finalList[[3]])
      min_y <- unlist(finalList[[4]])
      max_y <- unlist(finalList[[5]])
      
      y <- rep(1:dim_size, times = dim_size)
      x <- rep(1:dim_size, each = dim_size)
      hvtVolcanoMatrix <- rep(0, each = dim_size^2)
      for (k in names(finalList[[1]])) {
        eval(parse(
          text = paste0(
            "polygon_x <- round(scales::rescale(finalList[[1]]$`",
            k,
            "`$x, to=c(0, ", dim_size, "), from  = c(min_x,max_x)))"
          )
        ))
        eval(parse(
          text = paste0(
            "polygon_y <- round(scales::rescale(finalList[[1]]$`",
            k,
            "`$y, to=c(0, ", dim_size, "), from  = c(min_y,max_y)))"
          )
        ))
        present <- sp::point.in.polygon(x, y, polygon_x, polygon_y)
        eval(parse(
          text = paste0(
            "quantError <- finalList[[1]]$`",
            k,
            "`$value[1]"
          )
        ))
        hvtVolcanoMatrix[which(present != 0)] <- quantError
      }
      
      hvtVolcanoMatrix <- matrix(hvtVolcanoMatrix, nrow = dim_size)
      return(hvtVolcanoMatrix)
    }
    
    #####################################################################################
    
    temp_hvtVolcanoMatrix <- lapply(Level_list, depth_wise_surface)
    number_of_layers <- length(temp_hvtVolcanoMatrix)
    p <- plotly::plot_ly(showscale = FALSE)
    
    temp <- lapply(1:number_of_layers, function(i) {
      hvtVolcanoMatrix <- temp_hvtVolcanoMatrix[[i]] - ((i - 1) * separation_width)
      p <<- p %>%
        plotly::add_surface(z = ~hvtVolcanoMatrix, opacity = layer_opacity[i], name = paste("Layer_", i), showlegend = TRUE)
    })
    
    p <- p %>%
      plotly::config(displaylogo = FALSE) %>%
      plotly::layout(
        zaxis = list(title = hmap.cols),
        title = paste(
          "Hierarchical Voronoi Tessellation Interactive surface plot with",
          hmap.cols,
          "Heatmap Overlay"
        ),
        scene = list(zaxis = list(title = hmap.cols))
      )
    
    
    
    return(suppressMessages(p))
    
  } else {
    stop("Invalid value for 'plot.type'. Expected '1D','2Dproj', '2Dhvt','2Dheatmap', or 'surface_plot'.")
  }
}
