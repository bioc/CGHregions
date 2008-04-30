\name{CGHregions-package}
\alias{CGHregions-package}
\docType{package}
\title{
	Dimension Reduction for Array CGH Data with Minimal Information Loss.
}
\description{
	CGHregions takes as input array CGH data (log2-ratios) that have been segmented (i.e., split into chromosomal segments of similar log2-ratios) and called (i.e., a copy number assigned to each segment) on a per-sample basis and adjusts the segmentation so that break-points that are in similar locations across multiple samples are set to be in identical locations. Segmented and called data can be obtained by using the CGHcall package. The resulting dimensionality reduction facilitates downstream analysis in a variety of ways (e.g., reduces severity of multiple hypothesis testing, facilitates clustering and visualization, reduces computer memory requirements).
}
\details{
	\tabular{ll}{
	Package: \tab CGHregions\cr
	Type: \tab Package\cr
	Version: \tab 0.99.0\cr
	Date: \tab 2008-04-29\cr
	License: \tab GPL\cr
	}
}
\author{
	Mark van de Wiel and Sjoerd Vosse
	Maintainer: Mark van de Wiel <mark.vdwiel@vumc.nl>
}
\references{
	Mark A. van de Wiel and Wessel N. van Wieringen (2007). CGHregions: Dimension Reduction for Array CGH Data
	with Minimal Information Loss. \emph{Cancer Informatics, 2}, 55-63.
}
\keyword{ package }

