\name{rocplus}
\alias{rocplus}
\title{ROC, Precision Recall, Convex Hull and other plots}
\description{
 Creates single and multiple ROC plots, precision-recall, Bayes plots and others. The plots
 may be either rough or smooth with confidence limits. Threshold values, optimum points and 
 other interesting points may be marked on the plots. The convex hull for multiple
 curves with blending fractions for interesting points may be obtained. The axes may have
 normal or lognormal scaling.
}
\usage{
	rocplus(control,treatment,title="",
	texture=c("smooth","rough","both"),
	plot.type=c("roc","senbayes","specbayes","threshbayes","precrecall"),
	axes.type=c("prob","norm","lognorm"),	
	points=FALSE,
	bands=FALSE,
	prevalence,
	skew,
	optimum.weight,
	markedPoints,
	convex.hull,
	conf.level=0.95) 
}


\arguments{
	\item{control}{A vector, matrix, or list: see details.}
	\item{treatment}{A vector matrix or list: see details.}
	\item{title}{The title for the graph.}
	\item{texture}{A smooth plot or not.}
	\item{plot.type}{ROC or Bayes with Sensitivity or Bayes with Specificity or Bayes with thresholds or
			precision-recall}
	\item{axes.type}{Probability, Normal, Lognormal.}
	\item{points}{If \code{TRUE}, then threshold values will be marked on the curve.}
	\item{bands}{Show confidence bands on plot.}
	\item{prevalence}{Used by Bayes plots: the default is 0.5 if it is not specified. If it is specified 
			the optimum point will appear on the plots.}
	\item{skew}{is the ratio of population sizes for treatment over control. It sets prevalence to skew/(1+skew),
          and is available only for precision-recall calculations. If missing, prevalence 
		  will be set to length(control)/(length(control)+length(treatment) for precision-recall calculations.}
	\item{optimum.weight}{This will multiply treatment values when calculating the optimum point. 
		If prevalence is not set, setting  optimum.weight will make an optimum point appear on the plots.}
	\item{markedPoints}{(1) a vector or scalar of threshold values to be marked on the curve; or 
		(2) When convex.hull is TRUE, a vector of 1-Specificity values marking points on a convex hull}
	\item{convex.hull}{When TRUE, a convex hull will be drawn for multiple plots.}
	\item{conf.level}{A probability for significance tests and confidence intervals.}
}
\value{
	Either a ROC, a Bayes or a precision-recall plot will be drawn. In addition, a list will be 
		ouput with elements depending on which parameters are set. The elements of the list are:

	(1) AUC -- the AUC with a p-value and confidence limits.

	(2) optimum -- The optimum point with its location and a confidence interval on the vertical axis value.

	(3) TIES -- The percentage of ties in the combined control-treatment set.

	(4) Marked.points -- If there are threshold values, the list will contain a table of these, 
		their locations and confidence intervals on the vertical axis values. 

	(5) Marked.points -- When doing multiple plots and when	convex.hull is true, the list will contain
		a table indicating the way curves are blended to obtain the points.
}

\details{
	\code{control} and \code{treatment} may be vectors of differing lengths. A single 
		curve will be drawn on the plot.

	\code{control} and \code{treatment} may be matrices with the same number of columns, but
		the row lengths of the two may be different. Multiple curves will be drawn on the plot, one for 
		each column. If the control columns, have names, these will be drawn to mark the curves.

	\code{control} and \code{treatment} may be lists with the same number of elements in each:
		for example \code{control = list(a=1:3,b=c(2.7,8,3,9.1,10.5)}. The elements from each
		list will be successively processed as \code{control} and \code{treatment} to produce
		multiple curves on the plot. The curves will be labeled with the names from the 
		\code{control} list. Lists allow multiple curves to be draw when \code{control}
		and \code{treatment} have do not have the same number of values.

	The Bayes plots, show a Bayesian estimate on the vertical axis and a sensitivity, 
		1-Specificity or threshold value on the horizontal. For sensitivity, A, 1-specificity, B, and 
		prevalence, p, the Bayesian estimate is (Ap)/((Ap)+B(1-p). This is the ``precision'' 
		part of a precision-recall plot. The recall is of course 1-specificity. If p is the 
		prevalence of a disease in a population, then the Bayes estimate is an estimate of 
		the probability that an	an individual who exceeds the the corresponding threshold, 
		will have the disease. 
		

	The optimum point is the point on the curve whose likelihood, L, is equal to or greater than
		(1-prevalence)/prevalence which makes the optimum point the point at which the sum of the
		two errors, (1-sensitivity and 1-specificity), are minimized. See the vignette for 
		details. If an optimum.weight, w, is present, this will be taken as a weight on the treatment errors 
		making them worth w times the control errors. In this case the optimum point will be the 
		point at which the sum of (1-sensitivity)*optimum.weight and  (1-specificity) is minimized.
	
	Precision-recall sets prevalence to skew/(1+skew). It is in fact a relabeling of a Bayes-sensitivity
		plot. Users may find it convenient to extract random samples from the treatment and 
		control populations, and then uses skew to adjust the plots so that they are similar
		to those that use all data from both populations. If skew is not specified for a 
		precision-recall plot, it will be set, as given above, using the lengths of the two
		input vectors, which makes the plot the same as one in which precision is calculated in
		the usual way.
	
	Multiple plots with a convex hull and marked points will output a table with six columns. The
		row labels show the names of the curves blended for each marked point. The values alpha
		and 1-alpha are the blending proportions. The resulting Sensitivity and 1-Specificity
		are shown as well as the threshold values for the two curves. To achieve the blended 
		value in practice, one may sample the two curves, using the blending values as frequencies
		and the appropriate threshold values as cut-off criteria.	
	
	}

\note{
There is a vignette with more details and an example. To access it, type 

vignette("rocplus")

}
\examples{

	aA<-c(0,1,7,8,10,14,16,17,20,31)
	aB<-c(6,13,15,18,19,21,22,25,28,34)

	# A bare bones ROC plot
	rocplus(aA,aB)
	# A simple ROC plot with threshold values
	rocplus(aA,aB,points=TRUE)
	# The same with 95% confidence bands
	rocplus(aA,aB,points=TRUE,bands=TRUE)
	# or
	rocplus(aA,aB,po=TRUE,ban=TRUE)
	# The same plot, but rough not smooth
	rocplus(aA,aB,points=TRUE,bands=TRUE,texture="rough")
	# A plot with an optimum point
	rocplus(aA,aB,points=TRUE,bands=TRUE,prevalence=0.3)
	# And one with treatment errors weighted by four
	# Note: because prevalence is not specified, it is take as 0.5.
	rocplus(aA,aB,points=TRUE,bands=TRUE,optimum.weight=4)
	# Here it is again with prevalence specified
	rocplus(aA,aB,points=TRUE,bands=TRUE,prevalence=0.3,optimum.weight=4)


	# A Bayes plot with sensitivity on the x-axis. By default
	# prevalence is taken as 0.5
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="senbayes")
	# The same with a realistic prevalence
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="senbayes",prevalence=0.03)
	# One with 1-specificity on the x-axis
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="specbayes",prevalence=0.03)
	# A Bayes probability plot with thresholds on the x-axis.
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="threshbayes")

	# A Precision-Recall plot: aA and aB are the same length so skew=1
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="precrecall")
	# The same for skewed populations
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="precrecall",skew=4)
	


	# Changing the axes to normal
	rocplus(aA,aB,points=TRUE,bands=TRUE,plot.type="senbayes",axes.type="norm")

	# With two marked ponts
	rocplus(aA,aB,points=TRUE,bands=TRUE,markedPoints=c(10,22))

	# Compare a smooth plot with a rough plot
	rocplus(aA,aB,texture="both")

	# A multiple plot. 
	ym<-matrix(c(102,59,61,90,55,18,46,35,69,49,136,71,25,66,125,74,95,55,100,51,67,147,28,130),12,2)
	xm<-2*(ym[1:10,]-10)
	colnames(xm)<-c("One","Two")
	rocplus(xm,ym,"A multiple plot",points=FALSE) 

	# A multiple plot with a convex hull and marked points
	rocplus(xm,ym,"A multiple plot with a convex hull",convex.hull=TRUE,markedPoints=c(.2,.5,.05))

	# Using lists
	ym<-list(alpha=ym[,1],beta=ym[,2])
	xm<-list(alpha=xm[,1],beta=xm[,2])
	rocplus(xm,ym,"A multiple plot using lists") 

	# Distributions with different variances: Note the upward bend at the right.
	ym<-rnorm(100,2,4)
	xm<-rnorm(100)
	rocplus(xm,ym,"Differing variances",points=TRUE,bands=TRUE)

}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2011). rocplus \emph{rocplus}. The R project for statistical computing \url{http://www.r-project.org/}
}


\keyword{htest}

