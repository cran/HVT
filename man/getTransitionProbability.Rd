% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTransitionProbability.R
\name{getTransitionProbability}
\alias{getTransitionProbability}
\title{Creating Transition Probabilities list}
\usage{
getTransitionProbability(df, cellid_column, time_column)
}
\arguments{
\item{df}{Data frame. The input data frame should contain two columns, 
cell ID from scoreHVT function and time stamp of that dataset.}

\item{cellid_column}{Character. Name of the column containing cell IDs.}

\item{time_column}{Character. Name of the column containing time stamps.}
}
\value{
Prints and stores a nested list of data frames with transition probabilities.
}
\description{
This is the main function to create transition probabilities list.
The transition probability table quantifies the likelihood of transitioning from one state to another. 
States: The table includes the current states and the possible next states.
Probabilities: For each current state, it lists the probability of transitioning to each of the next possible states.
}
\examples{
dataset <- data.frame(date = as.numeric(time(EuStockMarkets)),
                      DAX = EuStockMarkets[, "DAX"],
                      SMI = EuStockMarkets[, "SMI"],
                      CAC = EuStockMarkets[, "CAC"],
                      FTSE = EuStockMarkets[, "FTSE"])
hvt.results<- trainHVT(dataset,n_cells = 60, depth = 1, quant.err = 0.1,
                       distance_metric = "L1_Norm", error_metric = "max",
                       normalize = TRUE,quant_method = "kmeans")
scoring <- scoreHVT(dataset, hvt.results)
cell_id <- scoring$scoredPredictedData$Cell.ID
time_stamp <- dataset$date
dataset <- data.frame(cell_id, time_stamp)
table <- getTransitionProbability(dataset, cellid_column = "cell_id",time_column = "time_stamp")
}
\seealso{
\code{\link{trainHVT}} \cr \code{\link{scoreHVT}}
}
\author{
PonAnuReka Seenivasan <ponanureka.s@mu-sigma.com>
}
\keyword{Transition_or_Prediction}
