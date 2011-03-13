# file rocplus.R
# copyright (C) 2011 by Robert E. Wheeler
#


"rocplus"<-function (control,treatment,title="",
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
	conf.level=0.95) {

# control, treatment contain vectors, matrices or lists
#  lists should be of the form list(a=control1,b=control2, etc) -- this allows control 
#		and treatment to be of different lengths
# title:  the title of the graph
# plot.type: ROC or Bayes with Sensitivity or Bayes with Specificity or Bayes with thresholds or
#			precision-recall
# axes.type: Probability, Normal, Lognormal
# points: will show labeled points on the plot
# bands: will put confidence bands about the line
# prevalence: Used by Bayes plots: the default is 0.5 if it is not specified. If it is specified 
#			  the optimum point will appear on the plots.
# skew: is the ratio of population sizes for treatment over control. It sets prevalence to skew/(1+skew),
#         and is available only for precision-recall calculations. If missing, prevalence 
#		  will be set to #control/(#control+#treatment) for precision-recall calculations.
# optimum.weight: will multiply treatment values when calculating the optimum point. It will cause
#		  the optimum point to appear on the plot even if prevalence is not set.
# markedPoints: (1) a vector or scalar of threshbayes values to be marked on the line,
#				(2) When convex.hull is TRUE, a vector of 1-Specificity values marking points on a convex hull
# convex.hull: when TRUE, a convex hull will be drawn for multiple plots.
# conf.level: used for AUC and the bands

# **************
# peJohnson
# extended pJohnson. Returns probs when X is in range and 0 when the argument is below the minimum or
#	1 when it is above
# *************

peJohnson<-function(X,parm) {
	if (parm$type=="SL" | parm$type=="SB") {
		minm<-parm$xi
		if (parm$type=="SL") {
			rng<-(((X-parm$xi)/parm$lambda)>1e-8)
				# a neg lambda should produce a 1 when the varible is out of range
			sch<-if(parm$lambda>0) !DECREASING else DECREASING
			pVal<-rep(as.numeric(sch),length(X))
			if (is.element(TRUE,rng))
				pVal[rng]<-pJohnson(X[rng],parm,lower.tail=DECREASING)
		} else {
			maxm<-parm$lambda+parm$xi
			rngH<-(X>=maxm)
			rngL<-(X<minm)
			rngM<-(X>minm & X<=maxm)
			pVal<-rep(0,length(X))
			pVal[rngH]<-as.numeric(DECREASING) # when lhs is empty,nothing happens
			pVal[rngL]<-as.numeric(!DECREASING)
			if (is.element(TRUE,rngM))
				pVal[rngM]<-pJohnson(X[rngM],parm,lower.tail=DECREASING)
		}
	
	} else {
		pVal<-pJohnson(X,parm,lower.tail=DECREASING)
	}	
	pVal
}




# **************
# get.adj
# Sets offset placement for labels
# **************
get.adj<-function() {
	switch(PLOT.TYPE,
		 roc = c(-.5,.5),
		 senbayes = c(1.5,1),
		 precrecall = c(1.5,1),
		 specbayes = c(-.5,-.5),
		 threshbayes = c(-.5,.5)
	)
}

# **************
# findLabelPts
# finds nPts equally spaced points on curve for labeling. Returns NULL if cant find enough pts
# **************
findLabelPts<-function(px,py,X,nPts){
	noLabel<-FALSE
	mv<-(is.na(px) | is.na(py))
	if (any(mv)) {  # must have more than nPts that are NA  
		if (sum(!mv)>nPts) {
			X<-X[!mv]
			px<-px[!mv]
			py<-py[!mv]
		} else
			noLabel<-TRUE
	}

	if (!noLabel) {
		n<-length(px)
		dx<-diff(px)
		dy<-diff(py)
		dd<-sqrt(dx^2+dy^2)
		ln<-sum(dd)
		stp<-ln/(nPts-1)

		inds<-findInterval((0:(nPts-2))*stp,cumsum(dd)) # Too close pts if last interval (nPts-2 to nPts-1) used
		
		pd<-c(TRUE,(diff(px[inds])>.02)) # get rid of too close points (mainly for convex hull)
		
		list(Y=(X[inds])[pd],px=(px[inds])[pd],py=(py[inds])[pd])
	} else NULL

}



# *************
# labelPoints
# Labels points on the plot. Returns NULL if the curve can't be labeled
# **************
labelPoints<-function(px,py,X) {
	lp<-findLabelPts(px,py,X,10)
	if (!is.null(lp)) {
		Y<-lp$Y
		px<-lp$px
		py<-lp$py

		if (!MULTIPLE) {  # can't have anything but prob for multiple axes
			if (AXES.TYPE=="norm") {
				if (PLOT.TYPE!="threshbayes")
					px<-qnorm(px)
				py<-qnorm(py)
			} 
			if (AXES.TYPE=="lognorm") {
				if (PLOT.TYPE!="threshbayes")
					px<-qlnorm(px)
				py<-qlnorm(py)
			} 
		}
		L<-signif(Y,SIGDIGITS)
		text(px,py,labels=L,adj=get.adj(),cex=.8)
		points(px,py)
	}
}


# **********************
# labelMultipleCurves
# In multiple plots, puts a label on the curve at a different point for each plot
#   Some curves may not have enough points to label
# **********************

labelMultipleCurves<-function(pC,pT,X){
	ls<-findLabelPts(pC,pT,X,(NPLOTS+4)) #to make sure pts not at edges
	if (!is.null(ls)) {
		px<-ls$px[PLOTNO+2] # NULL may be returned, leading to no point
		py<-ls$py[PLOTNO+2] # It's not worth fixing -- complicates the code for very special cases.

		if (AXES.TYPE=="norm") {
			if (PLOT.TYPE!="threshbayes")
				px<-qnorm(px)
			py<-qnorm(py)
		} 
		if (AXES.TYPE=="lognorm") {
			if (PLOT.TYPE!="threshbayes")
				px<-qlnorm(px)
			py<-qlnorm(py)
		} 

		
		if (!is.null(COLNAMES)) { # null can't happen, but leave it just in case I change my mind
			points(px,py*.995,pch="^",cex=1) # use .995 to move tip of ^ into contact with curve
			text(px,py,labels=sprintf("(%s)",COLNAMES[PLOTNO]),pos=1,cex=1)
		}
	}
}



# *************
# AUC
# area under the curve for treatment offset of delta
# Wilcoxon rank minus n(n+1)/2 is identical to Mann-Whitney
# *************
AUC<-function(delta=0){
	rk <- rank(c(TREATMENT+delta,CONTROL))
    auc<-sum(rk[seq_along(TREATMENT)]) - NTREATMENT*(NTREATMENT+1) / 2
	auc/(NCONTROL*NTREATMENT)
}

# **********
# translate
# translates values to Bayes if necessary
# ********** 
translate<-function(x,y,X) {
	if (PLOT.TYPE!="roc") {
		A<-y*PREVALENCE
		B<-x*(1-PREVALENCE)
		py<-A/(A+B)
		px<-switch(PLOT.TYPE, # horizontal axis values
			roc = x,
			senbayes = y, 
			precrecall = y,   
			specbayes = x,
			threshbayes = X
		)
		mv<-is.nan(py)
		if (any(mv)) {
			py[mv]<-NA
			px[mv]<-NA
		}
	} else {
		py<-y
		px<-x
	}

	list(px=px,py=py)
}

# *************
# translateDist
# translates axis distribution if needed
# *************
translateDist<-function(px,py) {
	if (AXES.TYPE=="norm") {
		if (PLOT.TYPE!="threshbayes")
			px<-qnorm(px)
		py<-qnorm(py)
	}
	if (AXES.TYPE=="lognorm") {
		if (PLOT.TYPE!="threshbayes")
			px<-qlnorm(px)
		py<-qlnorm(py)
	}
	list(px=px,py=py)
}



# *************
# bands.smooth
# Returns confidence bands or limits, based on AUC for  a smooth plot.
# ************
bands.smooth<-function(X) {
	ret.val<-NULL


	ci<-WMW$conf.int  # use wilcox.test 
	ci<-ci-WMW$estimate

	if (DECREASING)
		ci<-rev(ci)

	pmTreatment<-JohnsonFit(TREATMENT+ci[2])

	pControl<-peJohnson(X,PARMCONTROL)
	pTreatment<-peJohnson(X,pmTreatment)



	tv<-translate(pControl,pTreatment,X) # get axes values as probabilites for printing
	prx<-tv$px
	pry<-tv$py

	tv<-translateDist(prx,pry) # get values for plotting
	px<-tv$px
	py<-tv$py

	upper<-py
	upperp<-pry


	pmTreatment<-JohnsonFit(TREATMENT+ci[1])

	pControl<-peJohnson(X,PARMCONTROL)
	pTreatment<-peJohnson(X,pmTreatment)

	tv<-translate(pControl,pTreatment,X) # get axes values as probabilites for printing
	prx<-tv$px
	pry<-tv$py

	tv<-translateDist(prx,pry) # get values for plotting
	px<-tv$px
	py<-tv$py

	lower<-py
	lowerp<-pry

	list(ptsl=lower,ptsu=upper,ptsp=c(lowerp,upperp))

}



# **************
# findPts
# finds points in X at or adjacent to pts -- used only for rough plots. Assumes X is ordered.
# *************
findPts<-function(pts,X,px,py) {
	if (!DECREASING) X<-rev(X)
	indx<-findInterval(pts,X)
	indx[indx==0]<-NA  # findInterval may produce 0's
	if (!DECREASING) indx<-(length(X)+1)-indx
	list(px=px[indx],py=py[indx],indx=indx)
}

# *************
# bands.rough
# Returns  confidence bands or limits, based on AUC for a rough plot.
# ************

bands.rough<-function(Y) {

	ret.val<-NULL

	ci<-WMW$conf.int
	ci<-ci-WMW$estimate

	if (DECREASING)
		ci<-rev(ci)

	st<-sort(c(CONTROL,TREATMENT+ci[2]),index.return=TRUE,decreasing=!DECREASING)
	ix<-st$ix
	X<-st$x

	rTreatment<-c(rep(0,NCONTROL),rep(1/NTREATMENT,NTREATMENT))[ix]
	rControl<-c(rep(1/NCONTROL,NCONTROL),rep(0,NTREATMENT))[ix]
	rTreatment<-cumsum(rTreatment)
	rControl<-cumsum(rControl)

	fpts<-findPts(Y,X,rControl,rTreatment)
	rControl<-fpts$px
	rTreatment<-fpts$py


	tv<-translate(rControl,rTreatment,X) # get axes values as probabilites for printing
	prx<-tv$px
	pry<-tv$py

	tv<-translateDist(prx,pry) # get values for plotting
	px<-tv$px
	py<-tv$py


	upper<-py
	upperp<-pry

	st<-sort(c(CONTROL,TREATMENT+ci[1]),index.return=TRUE,decreasing=!DECREASING)
	ix<-st$ix
	X<-st$x

	rTreatment<-c(rep(0,NCONTROL),rep(1/NTREATMENT,NTREATMENT))[ix]
	rControl<-c(rep(1/NCONTROL,NCONTROL),rep(0,NTREATMENT))[ix]
	rTreatment<-cumsum(rTreatment)
	rControl<-cumsum(rControl)

	fpts<-findPts(Y,X,rControl,rTreatment) 
	rControl<-fpts$px
	rTreatment<-fpts$py

	tv<-translate(rControl,rTreatment,X) # get axes values as probabilites for printing
	prx<-tv$px
	pry<-tv$py

	tv<-translateDist(prx,pry) # get values for plotting
	px<-tv$px
	py<-tv$py



	lower<-py
	lowerp<-pry


	list(ptsl=lower,ptsu=upper,ptsp=c(lowerp,upperp))

}








# ************
# AUCTableOut
# prepares the stat output table for AUC and its confidence limits
# **************
AUCTableOut<-function(conf.level) {

	AUC<-round(WMW$statistic/(NCONTROL*NTREATMENT),3)
	pValue<-round(WMW$p.value/2,4) # wilcox.test puts out two sided p-values by default. Can't change
								   # because two sided confidence intervals are wanted,

	AUCci<-c(0,0)
	ci<-WMW$conf.int-WMW$estimate

		# confidence limits on AUC
	AUCci[1]<-AUC(ci[1])
	AUCci[2]<-AUC(ci[2])
	
	if (DECREASING) { # treatment and control are reversed in AUC() and Wilcox.Test()
		AUC<-1-AUC
		AUCci<-1-rev(AUCci)
	}

	AUCci<-signif(AUCci,SIGDIGITS)

	tab<-matrix(c(AUC,pValue,AUCci),1,4)
	pl<-(1-conf.level)/2
	pu<-1-pl
		
	colnames(tab)<-c("AUC","p-value",sprintf("%.3f%s",pl,"%"),sprintf("%.3f%s",pu,"%"))
	rownames(tab)<-"Mann-Whitney"
	tab
}




# ******************
# multParkPts
# Marks points on convex hull
# *****************
multMarkPts<-function(ch) {
	mkp<-MARKEDPOINTS
	mn<-min(CHULL[ch,1])
	mx<-max(CHULL[ch,1])
	mkp<-pmax(mkp,mn)
	mkp<-pmin(mkp,mx)
	nMark<-length(mkp)

	pt<-findInterval(mkp,CHULL[ch,1]) # the first column is in increasing order
	lng<-pt==length(ch)
	pt[lng]<-pt[lng]-1

	A<-CHULL[ch[pt],1]
	B<-CHULL[ch[pt+1],1]
	alpha<-round((B-mkp)/(B-A),2)
		# names
	A<-CHULL[ch[pt],4]
	B<-CHULL[ch[pt+1],4]
	nms<-sprintf("{%s | %s}",COLNAMES[A],COLNAMES[B])

		#values
	v1<-round(alpha*CHULL[ch[pt],1]+(1-alpha)*CHULL[ch[pt+1],1],2)
	v2<-round(alpha*CHULL[ch[pt],2]+(1-alpha)*CHULL[ch[pt+1],2],2)
		# thresholds
	t1<-signif(CHULL[ch[pt],3],SIGDIGITS)
	t2<-signif(CHULL[ch[pt+1],3],SIGDIGITS)

	points(v1,v2,pch=23,cex=1.3)

	alpha1<-(1-alpha)
	tab<-matrix(c(alpha,alpha1,v1,v2,t1,t2),nMark,6)
	colnames(tab)<-c("alpha","1-alpha","1-Specificity","Sensitivity","Left Threshold", "Right Threshold")
	rownames(tab)<-nms

	list(MarkedPts=tab)
}

# *************
# markPoints
# marks points specified by user on a curve and returns their values in a table for printing
# ************
markPoints<-function(px,py,prx,pry) {
	L<-signif(MARKEDPOINTS,SIGDIGITS)
	text(px,py,labels=L,cex=1.3,adj=get.adj())
	points(px,py,pch=23,cex=1.3)
	
		# get points scaled as probs for printing
	if (TEXTURE=="smooth")
		pts<-bands.smooth(MARKEDPOINTS)
	else
		pts<-bands.rough(MARKEDPOINTS)

	lg<-length(pts$ptsp)


	pl<-(1-conf.level)/2
	pu<-1-pl



	xlab<-switch(PLOT.TYPE,
		roc="1-Specificity",
		specbayes="1-Specifity",
		senbayes="Sensitivity",
		precrecall="Recall",
		threshbayes = "Threshold")

	lab<-"Bayes Pr. Est."

	ylab<-switch(PLOT.TYPE,
		roc="Sensitivity",
		specbayes=lab,
		precrecall="Precision",
		senbayes=lab,
		threshbayes = lab)

	nMark<-length(MARKEDPOINTS)
	if (lg!=2*nMark) {
		tab<-matrix(c(signif(MARKEDPOINTS,SIGDIGITS),round(c(prx,pry),3),rep(NA,nMark)),nMark,4)
		colnames(tab)<-c("value",xlab,ylab,"Can't calculate limits")
	} else {
		tab<-matrix(c(MARKEDPOINTS,round(c(prx,pry),2),round(pts$ptsp,2)),nMark,5)
		colnames(tab)<-c("value",xlab,ylab,sprintf("%.3f%s",pl,"%"),sprintf("%.3f%s",pu,"%"))
	}
	tab
	
}


#**************
# makeMultOpt
# Finds, draws and outputs a table for the optimum point on the convex hull
#**************
makeMultOpt<-function(chull) {
	
	x<-chull[,1]
	y<-chull[,2]
	X<-chull[,3]
	Curve<-chull[,4]
	pc<-0.05
	dx<-diff(x) # uses differences to find slope
	dy<-diff(y)
	n<-length(x)
	
		# Find a single optimum
	p1<-PREVALENCE
	p2<-OPTIMUM.WEIGHT
	dst<-abs(dy*p1*p2-dx*(1-p1))
	dst[is.na(dst)]<-1e6
	
	dst[1:trunc(pc*n)]<-1e6	 # don't include extreme tails in search for min
	dst[trunc((1-pc)*n):n]<-1e6

	op<-min(dst)
	xm<-match(op,dst)+1

	if (xm>=n) xm<-xm-1
	if (xm<=2) xm<-xm<-1


	alpha<-.5

		# names
	A<-Curve[xm]
	B<-Curve[xm+1]
	nms<-sprintf("{%s | %s}",COLNAMES[A],COLNAMES[B])

		#values
	v1<-round(alpha*x[xm]+(1-alpha)*x[xm+1],2)
	v2<-round(alpha*y[xm]+(1-alpha)*y[xm+1],2)
		# thresholds
	t1<-signif(X[xm],SIGDIGITS)
	t2<-signif(X[xm+1],SIGDIGITS)

	points(v1,v2,pch=23,cex=1.3)
	points(v1,v2,pch=19,cex=1.3)

	alpha1<-(1-alpha)
	tab<-matrix(c(alpha,alpha1,v1,v2,t1,t2),1,6)
	colnames(tab)<-c("alpha","1-alpha","1-Specificity","Sensitivity","Left Threshold", "Right Threshold")
	rownames(tab)<-nms

	list(optimum=tab)
}

# *************
# optPoint
#	Finds the optimum point with a value in X
# It does this by finding the point or points with slope equal to ((1-prevalence)/prevalence)
# *************
optPoint<-function(px,py,X) {
	pc<-0.05
	f1<-dJohnson(X,PARMTREATMENT)
	f2<-dJohnson(X,PARMCONTROL)
	
	f1[is.nan(f1)]<-NA  # derivative calc may produce NAN's
	f2[is.nan(f2)]<-NA

	n<-length(X)
	
		# Find a single optimum
	p1<-PREVALENCE
	p2<-OPTIMUM.WEIGHT
	dst<-abs(f1*p1*p2-f2*(1-p1))
	dst[is.na(dst)]<-1e6

	dst[1:trunc(pc*n)]<-1e6	 # don't include extreme tails in search for min
	dst[trunc((1-pc)*n):n]<-1e6

	xm<-match(min(dst),dst)

	optx<-px[xm]
	opty<-py[xm]
	optimumValue<-X[xm]

	list(optimumValue=optimumValue,optx=optx,opty=opty)
}
# **************
# makeOptimum
# outputs optimum table and prints optimum point on plot
# *************
makeOptimum<-function(pControl,pTreatment,X) {
	ls<-optPoint(pControl,pTreatment,X)
	ov<-ls$optimumValue
	optx<-ls$optx
	opty<-ls$opty

	tv<-translate(optx,opty,ov) # get axes values as probabilites for printing
	px<-tv$px
	py<-tv$py

	tv<-translateDist(px,py) # get values for plotting
	prx<-tv$px
	pry<-tv$py


	return.value<-NULL

	points(prx,pry,pch=19,cex=1.3) # a filled circle
	text(prx,pry,labels=signif(ov,SIGDIGITS),adj=get.adj(),cex=1.3)

	# get optimum pts scaled as probs for printing
	if (TEXTURE=="smooth")
		pts<-bands.smooth(ov)
	else
		pts<-bands.rough(ov)
	
	lg<-length(pts$ptsp)

	pl<-(1-conf.level)/2
	pu<-1-pl


	xlab<-switch(PLOT.TYPE,
		roc="1-Specificity",
		specbayes="1-Specifity",
		senbayes="Sensitivity",
		precrecall="Recall",
		threshbayes = "Threshold")

	lab<-"Bayes Pr. Est."

	ylab<-switch(PLOT.TYPE,
		roc="Sensitivity",
		specbayes=lab,
		senbayes=lab,
		precrecall="Precision",
		threshbayes = lab)

	nOpt<-length(ov) # there can be more than one optimum
	if (lg!=2*nOpt) {  # no confidence limits here
		value<-matrix(c(signif(ov,SIGDIGITS),round(c(prx,pry),3),rep(NA,nOpt)),nOpt,4)
		colnames(value)<-c("value",xlab,ylab,"Can't calculate limits")
	} else {
		value<-matrix(c(signif(ov,SIGDIGITS),round(c(prx,pry),3),round(pts$ptsp,2)),nOpt,5)
		colnames(value)<-c("value",xlab,ylab,sprintf("%.3f%s",pl,"%"),sprintf("%.3f%s",pu,"%"))
	}
	rownames(value)<-rep("Optimum Value",nOpt)	

	list(optimum=value)
	
}



# ************
# smoothCurve
# Uses Johnson curves to draw a smooth curve, and return parameters for printing
# ***********
smoothCurve<-function(points,bands,markPtsFlag) {

	X<-sort(c(CONTROL,TREATMENT),decreasing=!DECREASING)
	pControl<-peJohnson(X,PARMCONTROL)
	pTreatment<-peJohnson(X,PARMTREATMENT)

	tv<-translate(pControl,pTreatment,X) 
	px<-tv$px
	py<-tv$py

	tv<-translateDist(px,py)
	prx<-tv$px
	pry<-tv$py

	retVal<-NULL
	if (!MULTIPLE)
		retVal<-list(AUC=AUCTableOut(conf.level))

	if (TEXTURE=="smooth" || TEXTURE=="both") {
		if (MULTIPLE && CONVEX.HULL) {
			mv<-!(is.na(px) | is.na(py)) # axes.type must be "prob" here, soe prx not used.
			smv<-sum(mv)
			retVal<-c(list(CHULL=matrix(c(px[mv],py[mv],X[mv],rep(PLOTNO,smv)),smv,4)),retVal)
		}
		if (MULTIPLE & CONVEX.HULL) { # so that only the convex hull will show -- ie suppress individual plots
			if (PLOTNO==1) start.plot(0,0)
		} else {
			if (PLOTNO==1) {
				start.plot(prx,pry)
			} else {
				lines(prx,pry)
			}
		}
	}

	if (OPTIMUM.POINT && !MULTIPLE) {
		retVal<-c(retVal,makeOptimum(pControl,pTreatment,X))
	}

	if (MULTIPLE & !CONVEX.HULL) {
		labelMultipleCurves(px,py,X)
	}

	if (bands) {
		pts<-bands.smooth(X)
		mv<-(is.na(pts$ptsl) | is.na(pts$ptsu))
		pz<-prx[!mv]
		pty<-pts$ptsl[!mv]
		lines(pz,pty,lty="77")
		pty<-pts$ptsu[!mv]
		lines(pz,pty,lty="77")
	}

	if (points & !CONVEX.HULL) {
		labelPoints(px,py,X) 
	}

	if (markPtsFlag & !MULTIPLE) {
		mkp<-MARKEDPOINTS
		mn<-min(X)
		mx<-max(X)
		mkp<-pmax(mkp,mn)
		mkp<-pmin(mkp,mx)

		pC<-peJohnson(mkp,PARMCONTROL)
		pT<-peJohnson(mkp,PARMTREATMENT)


		tv<-translate(pC,pT,mkp) 
		prx<-tv$px
		pry<-tv$py

		tv<-translateDist(prx,pry)
		px<-tv$px
		py<-tv$py

			
		retVal<-c(retVal,list(Marked.points=markPoints(px,py,prx,pry)))
	}
	retVal

}


# ************
# roughCurve
# Uses Johnson curves to draw a rough curve, and return parameters for printing
# ***********
roughCurve<-function(points,bands,markPtsFlag){

	st<-sort(c(CONTROL,TREATMENT),index.return=TRUE,decreasing=!DECREASING)
	ix<-st$ix
	X<-st$x

	rTreatment<-c(rep(0,NCONTROL),rep(1/NTREATMENT,NTREATMENT))[ix]
	rControl<-c(rep(1/NCONTROL,NCONTROL),rep(0,NTREATMENT))[ix]
	rTreatment<-cumsum(rTreatment)
	rControl<-cumsum(rControl)


	tv<-translate(rControl,rTreatment,X) 
	px<-tv$px
	py<-tv$py

	tv<-translateDist(px,py)
	prx<-tv$px
	pry<-tv$py

	retVal<-NULL
	if (!MULTIPLE && texture=="rough")
		retVal<-list(AUC=AUCTableOut(conf.level))

	if (TEXTURE=="rough") {
		if (MULTIPLE && CONVEX.HULL) {
			mv<-!(is.na(px) | is.na(py)) # axes.type must be "prob" here, soe prx not used.
			smv<-sum(mv)
			retVal<-c(list(CHULL=matrix(c(px[mv],py[mv],X[mv],rep(PLOTNO,smv)),smv,4)),retVal)
		}
		if (MULTIPLE & CONVEX.HULL) {
			if (PLOTNO==1) start.plot(0,0)
		} else {
			if (PLOTNO==1) {
				start.plot(prx,pry)
			} else {
				lines(prx,pry)
			}
		}
	} else 
		if (TEXTURE=="both") {
		lines(prx,pry)	
	}

	if (OPTIMUM.POINT & !MULTIPLE) {
		retVal<-c(retVal,makeOptimum(rControl,rTreatment,X))
	}

	if (MULTIPLE & !CONVEX.HULL)
		labelMultipleCurves(rControl,rTreatment,X)
		
	if (bands) {
		pts<-bands.rough(X)
		mv<-(is.na(pts$ptsl) | is.na(pts$ptsu))
		pz<-prx[!mv]
		pty<-pts$ptsl[!mv]
		lines(pz,pty,lty="77")
		pty<-pts$ptsu[!mv]
		lines(pz,pty,lty="77")
	}


	if (points & !CONVEX.HULL)
		labelPoints(px,py,X) 


	if (markPtsFlag & !MULTIPLE) {
		mkp<-MARKEDPOINTS
		mn<-min(X)
		mx<-max(X)
		mkp<-pmax(mkp,mn)
		mkp<-pmin(mkp,mx)

		p<-findPts(mkp,X,rControl,rTreatment)
		pxt<-p$px
		pyt<-p$py

		tv<-translate(pxt,pyt,mkp) 
		prx<-tv$px
		pry<-tv$py

		tv<-translateDist(prx,pry)
		px<-tv$px
		py<-tv$py



		retVal<-c(retVal,list(Marked.points=markPoints(px,py,prx,pry)))
	}

	retVal
}

# *************
# start.plot
# Starts the plot with px and py
# ***************
start.plot<-function(px,py) {

	lab<-"Bayes Probability Estimate"

	if (PLOT.TYPE=="senbayes") {
		xlab<-"Sensitivity"
		ylab<-lab
	} else 
	if (PLOT.TYPE=="precrecall") {
		xlab<-"Recall"
		ylab<-"Precision"

	} else 
	if (PLOT.TYPE=="specbayes") {
		xlab<-"1-Specificity"
		ylab<-lab
	} else 
	if (PLOT.TYPE == "threshbayes") {
		xlab<-"Threshold"
		ylab<-lab
	} else {
		xlab<-"1-Specificity"
		ylab<-"Sensitivity"
	}



	if (AXES.TYPE=="norm") {
		xlim<-ylim<-c(-3,3)
		if (PLOT.TYPE!="roc") {
			if (PLOT.TYPE == "threshbayes") {
				xlim<-NULL
				plot(px,py,type="l",main=title,xlab=xlab,yaxt="n",ylab=ylab,xlim=xlim,ylim=ylim)
			} else 
				plot(px,py,type="l",main=title,xlab=xlab,xaxt="n",yaxt="n",ylab=ylab,xlim=xlim,ylim=ylim)
		} else {
			plot(px,py,type="l",main=title,xlab=xlab,xaxt="n",yaxt="n",ylab=ylab,xlim=xlim,ylim=ylim)
			lines(xlim,ylim,lty=3) # diagonal line
		}
		title(sub="Normal distribution scaling")
		if (PLOT.TYPE != "threshbayes")
			axis(1,at=seq(-3,3,by=1),label=round(pnorm(c(seq(-3,3,by=1))),2))
		axis(2,at=seq(-3,3,by=1),label=round(pnorm(c(seq(-3,3,by=1))),2))
	} else 
	if (AXES.TYPE=="lognorm") {
		xlim<-ylim<-c(0,6)
		if (PLOT.TYPE!="roc") {
			if (PLOT.TYPE == "threshbayes") {
				xlim<-NULL
				plot(px,py,type="l",main=title,xlab=xlab,yaxt="n",xlim=xlim,ylim=ylim,ylab=ylab)
			}
			else 
				plot(px,py,type="l",main=title,xlab=xlab,xaxt="n",yaxt="n",xlim=xlim,ylim=ylim,ylab=ylab)
		} else {
			plot(px,py,type="l",main=title,xlab=xlab,ylab=ylab,xaxt="n",yaxt="n",xlim=xlim,ylim=ylim)
			lines(xlim,ylim,lty=3) # diagonal line
		}
		title(sub="Log-normal distribution scaling")
		ax<-c(0,.5,1,2,3,4,5,6)
		al<-round(plnorm(ax),2)
		if (PLOT.TYPE != "threshbayes")
			axis(1,at=ax,label=al)
		axis(2,at=ax,label=al)

	} else {
		xlim<-ylim<-c(0,1)
		if (PLOT.TYPE!="roc") {
			if (PLOT.TYPE == "threshbayes") xlim<-NULL
			plot(px,py,type="l",main=title,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
		} else {
			plot(px,py,type="l",main=title,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
			lines(xlim,ylim,lty=3) # diagonal line
		}
	}
}





#************** Program Start *****************************************************
# Capitolized words are global

returnValue<-NULL

PLOT.TYPE<-match.arg(plot.type)
TEXTURE<-match.arg(texture)
AXES.TYPE<-match.arg(axes.type)

SIGDIGITS<-3
MULTIPLE<-(is.matrix(control) || is.list(control))

CONVEX.HULL<-if (!missing(convex.hull)) convex.hull else FALSE

if (CONVEX.HULL && PLOT.TYPE !="roc")
	stop("Only ROC plots may be made when convex.hull is true.")
if (CONVEX.HULL && AXES.TYPE != "prob")
	stop("axes.type must be 'prob' when convex.hull is true.")
if (PLOT.TYPE!="precrecall" & !missing(skew))
	returnValue<-c(returnValue,list(ERROR="skew ignored"))

doSmooth<-doRough<-FALSE

if (!missing(prevalence)) {
	if (prevalence>(1-1e-6) || prevalence < 1e-6)
		stop("prevalence should be in the range (0,1)")
	PREVALENCE<-prevalence
	OPTIMUM.POINT<-TRUE
} else {
	OPTIMUM.POINT<-FALSE
	PREVALENCE<-0.5
	if (PLOT.TYPE!="roc")
		returnValue<-c(returnValue,list(Note="prevalence set to 0.5"))
}

if (!missing(optimum.weight)){
	if (optimum.weight<=0)
		stop("optimum.weight must be positive")
	OPTIMUM.WEIGHT<-optimum.weight
	OPTIMUM.POINT<-TRUE
} else {
	OPTIMUM.WEIGHT<-1
}

if (PLOT.TYPE=="precrecall") { # skew has meaning only for precision recall
	if (missing(skew)) {
		nc<-length(control)
		nt<-length(treatment)
		PREVALENCE<-nc/(nc+nt)
		returnValue<-c(returnValue,list(Note=sprintf("skew set to %.3f",nc/nt)))
	}
	else {
		if (skew<=0)
			stop("skew must be positive")
		PREVALENCE<-skew/(1+skew)	
	}
	OPTIMUM.WEIGHT<-1
}




markPtsFlag<-(!missing(markedPoints)) 
if (markPtsFlag)
	MARKEDPOINTS<-markedPoints

if (TEXTURE=="smooth" || TEXTURE=="both") 
	doSmooth<-TRUE
if (TEXTURE=="rough" || TEXTURE=="both")
	doRough<-TRUE


doingList<-FALSE
if (MULTIPLE) {
	CHULL<-NULL
	AXES.TYPE<-"prob"
	if (is.list(control)) {
		if (!is.list(treatment))
			stop("Both control and treatment must be of the same class")
		NPLOTS<-length(control)
		nT<-length(treatment)
		if (NPLOTS!=nT)
			stop("Both control and treatment must be the same length")
		COLNAMES<-names(control)
		doingList<-TRUE
	} else {
		NPLOTS<-ncol(control)
		nT<-ncol(treatment)
		if (nT != NPLOTS) 
			stop("Dimensions of control and treatment do not match")
		COLNAMES<-colnames(control)
	}
	if (is.null(COLNAMES)) COLNAMES<-1:NPLOTS # name the curves if necessary
	bands<-FALSE
} else {
	if (!is.numeric(control) || !is.numeric(treatment))
		stop("control and treatment must be of the same class")
	NPLOTS<-1
}

# ************** Begin looping
for (PLOTNO in 1:NPLOTS) {
	if (PLOTNO>1) {
		bands<-FALSE
	}

	if (MULTIPLE) {
		if (doingList) {
			CONTROL<-control[[PLOTNO]]
			CONTROL<-CONTROL[!is.na(CONTROL)]
			TREATMENT<-treatment[[PLOTNO]]
			TREATMENT<-TREATMENT[!is.na(TREATMENT)]
		} else {
			CONTROL<-control[!is.na(control[,PLOTNO]),PLOTNO]
			TREATMENT<-treatment[!is.na(treatment[,PLOTNO]),PLOTNO]
		}
	} else {
		CONTROL<-control[!is.na(control)]
		TREATMENT<-treatment[!is.na(treatment)]
	}

	NCONTROL<-length(CONTROL)
	NTREATMENT<-length(TREATMENT)

	PARMCONTROL<-JohnsonFit(CONTROL)
	PARMTREATMENT<-JohnsonFit(TREATMENT)


		# use the upper tail if the threshbayes values are increasing
	DECREASING<-(median(CONTROL)>median(TREATMENT))

	WMW<-suppressWarnings(wilcox.test(TREATMENT,CONTROL,pval=TRUE,conf.int=TRUE,
			conf.level=conf.level))

	if (doSmooth) {
		val<-smoothCurve(points,bands,markPtsFlag)
		if (MULTIPLE && CONVEX.HULL) { # collect data for convex hull calculation
			CHULL<-rbind(CHULL[,1:4],val$CHULL)
			val<-val[-1]
		}
		returnValue<-c(returnValue,val)
	}

	if (doRough) {
		val<-roughCurve(points,bands,markPtsFlag)
		if (MULTIPLE && CONVEX.HULL) {
			CHULL<-rbind(CHULL[,1:4],val$CHULL)
			val<-val[-1]
		}
		returnValue<-c(returnValue,val)
	}

} # ******************** end looping

if (MULTIPLE && CONVEX.HULL) { 
	ch<-chull(CHULL[,1:2])

	ch<-ch[-(1:2)] # convex hull closure points usually

	ix<-sort(CHULL[ch,1],index.return=TRUE)$ix # curves may have different threshold sets
	ch<-ch[ix]

	lines(CHULL[ch,1:2],lty="FAFA")

	if (points) {
		labelPoints(CHULL[ch,1],CHULL[ch,2],CHULL[ch,3]) 
	}


	if (OPTIMUM.POINT)
		returnValue<-c(returnValue,makeMultOpt(CHULL[ch,]))

	if (markPtsFlag)
		returnValue<-c(returnValue,multMarkPts(ch))

	title(sub="Convex Hull Plot")

}

if (!MULTIPLE) {
	nTies<-(NCONTROL+NTREATMENT)-length(unique(c(CONTROL,TREATMENT)))
	returnValue$TIES=sprintf("%.1f%s",100*nTies/(NCONTROL+NTREATMENT),"%")
}

	if (length(returnValue)>0)
		returnValue

}
