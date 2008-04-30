\name{CGHregions}
\alias{CGHregions}
\title{ Dimension Reduction for Array CGH Data with Minimal Information Loss }
\description{
  Dimension Reduction for Array CGH Data with Minimal Information Loss
}
\usage{
	CGHregions(input, averror=0.01)
}
\arguments{
  \item{input}{ An object of class \code{\link[CGHbase:class.cghCall]{cghCall}}, a character string or a dataframe. See details for information on the latter two.}
  \item{averror}{ Maximal information loss allowed. }
 }
\details{
    Please read the article and the supplementary information for detailed information on the algorithm.
    
    If the input is not an object of class \code{\link[CGHbase:class.cghCall]{cghCall}} it should be either a dataframe or a tabseparated textfile (textfiles must contain a header). The first three columns should contain the name, chromosome and position in bp for each array target respectively. The chromosome and position column must contain numbers only. Following these is a column with log2 ratios for each of your samples. If the input type is a textfile, missing values should be represented as 'NA' or an empty field.
    
    The algorithm reduces the call matrix to a smaller matrix that contains regions rather than individual clones. 
    The regions consist of consequtive clones the signatures of which are very much alike. The dimension reduction 
    is potentially for testing and clustering puposes. The amount of information lost by this dimension 
    reduction is controlled by \code{averror}. The larger \code{averror}, the less regions will result.
    
}
\value{
    This function returns an object of class \code{\link[CGHbase:class.cghRegions]{cghRegions}}
}
\references{
	Mark A. van de Wiel and Wessel N. van Wieringen (2007). CGHregions: Dimension Reduction for Array CGH Data
	with Minimal Information Loss. \emph{Cancer Informatics, 2}, 55-63.
}
\author{ 
	Mark van de Wiel and Sjoerd Vosse
	Maintainer: Mark van de Wiel <mark.vdwiel@vumc.nl>
}
\examples{
	data(WiltingCalled)
	result <- CGHregions(WiltingCalled)
}
\keyword{ misc }