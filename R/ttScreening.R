irwsva.build2<- 
function (dat, mod, mod0 = NULL, n.sv, B = 5)  
{ 
    n <- ncol(dat) 
    m <- nrow(dat) 
    if (is.null(mod0)) { 
        mod0 <- mod[, 1] 
    } 
    Id <- diag(n) 
    resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%  
        t(mod)) 
    uu <- eigen(t(resid) %*% resid) 
    vv <- uu$vectors 
    ndf <- n - dim(mod)[2] 
    pprob <- rep(1, m) 
    one <- rep(1, n) 
    Id <- diag(n) 
    df1 <- dim(mod)[2] + n.sv 
    df0 <- dim(mod0)[2] + n.sv 
    rm(resid) 
    cat(paste("Iteration (out of", B, "):")) 
    for (i in 1:B) { 
        mod.b <- cbind(mod, uu$vectors[, 1:n.sv]) 
        mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv]) 
        ptmp <- f.pvalue(dat, mod.b, mod0.b) 
        pprob.b <- (1 - as.function(get("edge.lfdr",envir=environment(sva)))(ptmp)) 
        mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv]) 
        mod0.gam <- cbind(mod0) 
        ptmp <- f.pvalue(dat, mod.gam, mod0.gam) 
        pprob.gam <- (1 - as.function(get("edge.lfdr",envir=environment(sva)))(ptmp)) 
        pprob <- pprob.gam * (1 - pprob.b) 
        dats <- dat * pprob 
        dats <- dats - rowMeans(dats) 
        uu <- eigen(t(dats) %*% dats) 
        cat(paste(i, " ")) 
    } 
    # Patch code 
    if(any(dats!=0)) {sv = fast.svd(dats, tol = 0)$v[, 1:n.sv]
    } else {sv=svd(dats)$v[, 1:n.sv]
		print("error in fast.svd(); svd() applied instead")}

    retval <- list(sv = sv, pprob.gam = pprob.gam, pprob.b = pprob.b,  
        n.sv = n.sv,message) 
    return(retval) 
} 

sva2=function (dat, mod, mod0 = NULL, n.sv = NULL, method = c("irw",  
    "two-step"), vfilter = NULL, B = 5, numSVmethod = "be")  
{ 
    if (is.null(n.sv)) { 
        n.sv = num.sv(dat, mod, method = numSVmethod, vfilter = vfilter) 
    } 
    if (!is.null(vfilter)) { 
        if (vfilter < 100 | vfilter > dim(dat)[1]) { 
            stop(paste("The number of genes used in the analysis must be between 100 and",  
                dim(dat)[1], "\n")) 
        } 
        tmpv = rowVars(dat) 
        ind = which(rank(-tmpv) < vfilter) 
        dat = dat[ind, ] 
    } 
    if (n.sv > 0) { 
        cat(paste("Number of significant surrogate variables is: ",  
            n.sv, "\n")) 
        method <- match.arg(method) 
        if (method == "two-step") { 
            return(twostepsva.build(dat = dat, mod = mod, n.sv = n.sv)) 
        } 
        if (method == "irw") { 
            return(irwsva.build2(dat = dat, mod = mod, mod0 = mod0, n.sv = n.sv, B = B)) 
        } 
    } 
    else { 
        print("No significant surrogate variables") 
        return(list(sv = 0, pprob.gam = 0, pprob.b = 0, n.sv = 0)) 
    } 
} 


ttScreening<- function(y=y,x1=x1,x2=x2,type=c("numeric","factor"),interaction=c(TRUE,FALSE),iterations=100,
					cv.cutoff=85,n.sv=NULL,train.alpha=0.05,test.alpha=0.1, 
					FDR.alpha=0.05,Bon.alpha=0.05,percent=(2/3), method =c("TT","FDR","Bonferroni"),
					linear= c("robust","ls"),vfilter = NULL, B = 5, numSVmethod = "be",rowname=NULL){

	if (is.null(rowname) == TRUE){
	rownames(y) <- seq(1:nrow(y))
	} else {rownames(y) <- rownames(y)}	

	if(is.null(x2) == FALSE){x2<-as.factor(x2)}else{interaction = 'FALSE'}
	if(type == "factor"){x1<-as.factor(x1)}

	if (is.null(x2) == TRUE){
		
		y.omit <- which(is.na(y), arr.ind=TRUE)[,1]		#gives rows of missing cpg data
		x1.omit <- which(is.na(x1), arr.ind=TRUE) 		#gives missing covariate data
		if (length(y.omit)<1  && length(x1.omit) < 1){
		data=y
		x1=x1 }
		if(length(y.omit)<1 && length(x1.omit) > 1){
		data <- y[,-x1.omit]
		x1 <- x1[-x1.omit]}
		if(length(y.omit)>1 && length(x1.omit) > 1){
		data <- y[-y.omit,-x1.omit]
		x1 <- x1[-x1.omit]}
		if(length(y.omit)>1 && length(x1.omit) < 1){
		data <- y[-y.omit,]
		x1 <- x1}
		
	}else{

		y.omit <- which( is.na(y), arr.ind=TRUE)[,1] 		#gives rows of missing cpg data
		x1.omit <- which( is.na(x1), arr.ind=TRUE)		#gives missing covariate data
		x2.omit <- which( is.na(x2), arr.ind=TRUE)		#gives missing snp data
		if (length(y.omit) < 1 && length(x1.omit) < 1 && length(x2.omit) < 1){
		data=y
		x1 = x1
		x2 = x2}
		if(length(y.omit) > 1 && length(x1.omit) < 1 && length(x2.omit) < 1){
		data=y[-y.omit,]
		x1 = x1
		x2 = x2}
		if(length(y.omit)<1 && length(x1.omit) > 1 | length(x2.omit) > 1){
		omit <- as.numeric(na.omit(c(x1.omit,x2.omit)))
		data <- y[,-omit]
		x1 <- x1[-omit]
		x2 <- x2[-omit]}
		if(length(y.omit)>1 && length(x1.omit) > 1 | length(x2.omit) > 1){
		omit <- as.numeric(na.omit(c(x1.omit,x2.omit)))
		data <- y[-y.omit,-omit]
		x1 <- x1[-omit]
		x2 <- x2[-omit]}

		}

	if(is.null(x2) == TRUE && type == "factor"){
	if(length(levels(x1)) > 3){stop("the level of x1 can not be more than 3")
	geterrmessage()}}


	if(is.null(x2) == FALSE && type == "factor"){
	if( length(levels(x2)) > 3 & length(levels(x1)) > 2){stop("x1 and x2 are both categorical variables. x1 has to be binary and x2 cannot have more than 3 levels to avoid singular design matrix")
	geterrmessage()}}
	
	if(as.numeric(max(data)<1 && min(data)>0) == 1){edata <- log2(data/(1-data))}else{
		   edata<-data}

	train.length<-test.length<-TT.output<-FDR.output<-Bon.output<-final.temp<-freq.temp<-cpg.select<-NULL
	selection<-matrix(rep(0,nrow(edata)*iterations),nrow=nrow(edata),ncol=iterations)
	pvalue.matrix<-matrix(rep(NA,nrow(edata)*iterations),nrow=nrow(edata),ncol=iterations)
	
	

	if(is.null(x2) == TRUE){
	
	if(interaction == FALSE && type == "numeric"){
		mod<-model.matrix(~x1, data=data.frame(x1))[,1:2]}

	if(interaction == FALSE && type == "factor"){
		mod<-model.matrix(~as.factor(x1), data=data.frame(x1))[,1:(length(unique(x1)))]}
			
		mod0<-model.matrix(~1,data=data.frame(x1))
		
		svobj = sva2(as.matrix(edata),mod,mod0,n.sv=n.sv)
		n.sv.temp<-svobj$n.sv
		temp=data.frame(svobj$sv)
		location<-NULL
		if(n.sv.temp > 0){
			for (w in 1:n.sv.temp){
				if(length(unique(temp[,w])) == 2){
				location<-c(location,w)}else{location<-location}
			}
		}
		if(is.null(location)){svobj$sv<-temp
		}else {svobj$sv<- as.numeric(temp[,-location])
		print(paste(length(location),"surrogate variables were not included in model fitting due to their lack of information"))
		svobj$n.sv<-ncol(svobj$sv)}
		if(n.sv.temp == 0 | is.null(svobj$n.sv) == TRUE){
		modSv = mod
		mod0Sv = mod0
		svobj$n.sv <- 0}else{
		modSv = cbind(mod,as.matrix(svobj$sv))
		mod0Sv = cbind(mod0,as.matrix(svobj$sv))}


if(method == "TT"){
	for (i in 1:iterations){
		set.seed(i*30)
		x <- 1:ncol(edata)
		length <- round(percent * ncol(edata))  
		train <- sample(x, size=length, replace=FALSE)
		test <- x[-train]

		iter<-0
		repeat{
		iter<-iter+1
        	train <- sample(x, size=length, replace=FALSE)
	  	test <- x[-train]
        	if(is.matrix(try(solve(t(modSv[train,])%*%modSv[train,]),silent=TRUE))==TRUE & is.matrix(try(solve(t(modSv[test,])%*%modSv[test,]),silent=TRUE))==TRUE){break}    
		if(iter == 10){stop("sample size too small or solve(t(mod) %*% mod) is singular")
		geterrmessage()}
		}
			
		train.temp<-which(f.pvalue(as.matrix(edata[,train]),modSv[train,],mod0Sv[train,]) <= train.alpha)
		train.length<-c(train.length,length(train.temp))	#number of sites selected at training level
		if(length(train.temp) < 1){
			test.temp<-NULL}
			test.temp<-train.temp[(f.pvalue(matrix(edata[train.temp,test],nrow=length(train.temp),ncol=length(test)),modSv[test,],mod0Sv[test,]) <= test.alpha)]
		test.length<-c(test.length,length(test.temp))		#number of sites selected at testing level
		if(length(test.temp) < 1){selection[,i]<- 0}
		if(length(test.temp) == 1){selection[test.temp,i]<- 1
						  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(rbind(edata[test.temp,],rep(0,ncol(edata)))),modSv,mod0Sv)[1]}
		if(length(test.temp) > 1){selection[test.temp,i]<- 1
				  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(edata[test.temp,]),modSv,mod0Sv)}
		}}	
	}
		

	if(is.null(x2) == FALSE && interaction == FALSE){
		
	if(interaction == FALSE && type == "numeric"){
		mod <- model.matrix(~x1 + as.factor(x2), data=data.frame(x1,x2))[,1:(length(unique(x2))+1)]
		mod0 <- model.matrix(~x1, data=data.frame(x1))[,1:2]}

	if(interaction == FALSE && type == "factor"){
		mod <- model.matrix(~as.factor(x1) + as.factor(x2), data=data.frame(x1,x2))[,1:(length(unique(x1))+length(unique(x2))-1)]
		mod0 <- model.matrix(~as.factor(x1) , data=data.frame(x1))[,1:(length(unique(x1)))]}

		
		svobj = sva2(as.matrix(edata),mod,mod0,n.sv=n.sv)
		n.sv.temp<-svobj$n.sv
		temp=data.frame(svobj$sv)
		location<-NULL
		if(n.sv.temp > 0){
			for (w in 1:n.sv.temp){
				if(length(unique(temp[,w])) == 2){
				location<-c(location,w)}else{location<-location}
			}
		}
		if(is.null(location)){svobj$sv<-temp
		}else {svobj$sv<- as.numeric(temp[,-location])
		print(paste(length(location),"surrogate variables were not included in model fitting due to their lack of information"))
		svobj$n.sv<-ncol(svobj$sv)}
		if(n.sv.temp == 0 | is.null(svobj$n.sv) == TRUE){
		modSv = mod
		mod0Sv = mod0
		svobj$n.sv <- 0}else{
		modSv = cbind(mod,as.matrix(svobj$sv))
		mod0Sv = cbind(mod0,as.matrix(svobj$sv))}



if(method == "TT"){
	for (i in 1:iterations){
		set.seed(i*30+6*i)
		x <- 1:ncol(edata)
		length <- round(percent * ncol(edata))  
		train <- sample(x, size=length, replace=FALSE)
		test <- x[-train]
		
		iter<-0
		repeat{
		iter<-iter+1
        	train <- sample(x, size=length, replace=FALSE)
	  	test <- x[-train]
        	if(is.matrix(try(solve(t(modSv[train,])%*%modSv[train,]),silent=TRUE))==TRUE & is.matrix(try(solve(t(modSv[test,])%*%modSv[test,]),silent=TRUE))==TRUE){break}    
		if(iter == 10){stop("sample size too small or solve(t(mod) %*% mod) is singular")
		geterrmessage()}
		}


		train.temp<-which(f.pvalue(as.matrix(edata[,train]),modSv[train,],mod0Sv[train,]) <= train.alpha)
		train.length<-c(train.length,length(train.temp))	#number of sites selected at training level
		if(length(train.temp) < 1){
			test.temp<-NULL}
			test.temp<-train.temp[(f.pvalue(matrix(edata[train.temp,test],nrow=length(train.temp),ncol=length(test)),modSv[test,],mod0Sv[test,]) <= test.alpha)]
		test.length<-c(test.length,length(test.temp))		#number of sites selected a testing level
		if(length(test.temp) < 1){selection[,i]<- 0}
		if(length(test.temp) ==1){selection[test.temp,i]<- 1
						  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(rbind(edata[test.temp,],rep(0,ncol(edata)))),modSv,mod0Sv)[1]}
		if(length(test.temp) > 1){selection[test.temp,i]<- 1
				  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(edata[test.temp,]),modSv,mod0Sv)}
		}}
	}



	if(is.null(x2) == FALSE && interaction == TRUE){
	
	if(interaction == TRUE && type == "numeric"){
		mod <- model.matrix(~x1 + as.factor(x2)+ x1:as.factor(x2), data=data.frame(x1,x2))[,1:(2*length(unique(x2)))]
		mod0 <- model.matrix(~x1 + as.factor(x2), data=data.frame(x1,x2))[,1:(length(unique(x2))+1)]}

	if(interaction == TRUE && type == "factor"){
		mod <- model.matrix(~as.factor(x1) + as.factor(x2)+ as.factor(x1):as.factor(x2), data=data.frame(x1,x2))[,1:(length(unique(x1))*length(unique(x2)))]
		mod0 <- model.matrix(~as.factor(x1) + as.factor(x2), data=data.frame(x1,x2))[,1:(length(unique(x1))+length(unique(x2))-1)]}
		
		svobj = sva2(as.matrix(edata),mod,mod0,n.sv=n.sv)
		n.sv.temp<-svobj$n.sv
		temp=data.frame(svobj$sv)
		location<-NULL
		if(n.sv.temp > 0){
			for (w in 1:n.sv.temp){
				if(length(unique(temp[,w])) == 2){
				location<-c(location,w)}else{location<-location}
			}
		}
		if(is.null(location)){svobj$sv<-temp
		}else {svobj$sv<- as.numeric(temp[,-location])
		print(paste(length(location),"surrogate variables were not included in model fitting due to their lack of information"))
		svobj$n.sv<-ncol(svobj$sv)}
		if(n.sv.temp == 0 | is.null(svobj$n.sv) == TRUE){
		modSv = mod
		mod0Sv = mod0
		svobj$n.sv <- 0}else{
		modSv = cbind(mod,as.matrix(svobj$sv))
		mod0Sv = cbind(mod0,as.matrix(svobj$sv))}



if(method == "TT"){
	for (i in 1:iterations){
		set.seed(i*30+6*i)
		x <- 1:ncol(edata)
		length <- round(percent * ncol(edata))  
		train <- sample(x, size=length, replace=FALSE)
		test <- x[-train]

		iter<-0
		repeat{
		iter<-iter+1
        	train <- sample(x, size=length, replace=FALSE)
	  	test <- x[-train]
        	if(is.matrix(try(solve(t(modSv[train,])%*%modSv[train,]),silent=TRUE))==TRUE & is.matrix(try(solve(t(modSv[test,])%*%modSv[test,]),silent=TRUE))==TRUE){break}    
		if(iter == 10){stop("sample size too small or solve(t(mod) %*% mod) is singular")
		geterrmessage()}
		}


		train.temp<-which(f.pvalue(as.matrix(edata[,train]),modSv[train,],mod0Sv[train,]) <= train.alpha)
		train.length<-c(train.length,length(train.temp))	#number of sites selected at the training level
		if(length(train.temp) < 1){
			test.temp<-NULL}
			test.temp<-train.temp[(f.pvalue(matrix(edata[train.temp,test],nrow=length(train.temp),ncol=length(test)),modSv[test,],mod0Sv[test,]) <= test.alpha)]
		test.length<-c(test.length,length(test.temp))		#number of sites selected at the testing level
		if(length(test.temp) < 1){selection[,i]<- 0}
		if(length(test.temp) ==1){selection[test.temp,i]<- 1
						  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(rbind(edata[test.temp,],rep(0,ncol(edata)))),modSv,mod0Sv)[1]}
		if(length(test.temp) > 1){selection[test.temp,i]<- 1
				  pvalue.matrix[test.temp,i]<- f.pvalue(as.matrix(edata[test.temp,]),modSv,mod0Sv)}
		}}
	}


if(method == "TT"){
		cutoff <- ((cv.cutoff/100)*iterations)
		cpg.select<-rowSums(selection)
		freq.temp<-rowSums(selection)[rowSums(selection) >= cutoff]
		final.temp<-which(rowSums(selection) >= cutoff)
}		

	if(method == "FDR" | method == "Bonferroni"){
		if(linear == "robust"){
		one<-tryCatch(lmFit(edata, design=modSv, method="robust"), warning=function(w)finally=print("robust regression did not converge; LS regression applied instead"))
		}else{one<-0}
		if(linear=="robust" & is.character(one) == TRUE){linear = "ls"}
		lmfit<-eBayes(lmFit(edata, design=modSv, method=linear))$p.value

	if(is.null(x2) == TRUE){
		if(type == "numeric" || length(levels(x1))==2){lmfit.min<-lmfit[,2]}else{lmfit.min<-apply(lmfit[,2:3],1,min)}
		lmFitFDR.rob<-which(p.adjust(lmfit.min,method="fdr")<= FDR.alpha)
		lmFitBon.rob<-which(p.adjust(lmfit.min,method="bonferroni")<= Bon.alpha)
		}

	if(interaction == FALSE && is.null(x2) == FALSE){
		if(length(levels(x2)) == 2){lmfit.min<-lmfit[,3]}else{lmfit.min<-apply(lmfit[,3:4],1,min)}
		lmFitFDR.rob<-which(p.adjust(lmfit.min,method="fdr")<= FDR.alpha)
		lmFitBon.rob<-which(p.adjust(lmfit.min,method="bonferroni")<= Bon.alpha)
		}

	if(interaction == TRUE){
		if(length(unique(x2))==2){lmfit.min<-lmfit[,4]}else{lmfit.min<-apply(lmfit[,5:6],1,min)}
		lmFitFDR.rob<-which(p.adjust(lmfit.min,method="fdr")<= FDR.alpha)
		lmFitBon.rob<-which(p.adjust(lmfit.min,method="bonferroni")<= Bon.alpha)
		}
	}else{lmFitFDR.rob<-lmFitBon.rob<-NULL}

if(method == "TT"){
	if(length(final.temp) > 1){
		if(linear == "robust"){
		one<-tryCatch(lmFit(edata[final.temp,], design=modSv, method="robust"), warning=function(w)finally=print("robust regression did not converge; LS regression applied instead"))
		}else{one<-0}
		if(linear=="robust" & is.character(one)==TRUE){linear = "ls"}

		tt<-eBayes(lmFit(edata[final.temp,],design=modSv,method=linear))
		TT.output<-data.frame(rownames(edata)[final.temp],freq.temp,tt$coefficients,tt$p.value)}
	if(length(final.temp) == 1){
		tt<- summary(glm(t(as.matrix(edata[final.temp,]))~ modSv[,-1],family=gaussian))
		TT.output<-data.frame(rownames(edata)[final.temp],freq.temp,t(tt$coefficients[,1]),t(tt$coefficients[,4]))}
	if(length(final.temp) < 1){
		TT.output<- t(rep("NA",(ncol(modSv)*2 + 2)))}
	
	colnames(TT.output)<-c("Row_Ind_Select","Selection Prop","Int_Coeff",rep("Var_Coeff",ncol(mod)-1),
						rep("SV_Coeff",svobj$n.sv),"Int_Pvalue",rep("Var_Pvalue",ncol(mod)-1),
						rep("SV_Pvalue",svobj$n.sv))
}


	if(method == "FDR" | method == "Bonferroni"){
	if(length(lmFitFDR.rob) > 1){
		ff<-eBayes(lmFit(edata[lmFitFDR.rob,],design=modSv,method=linear))
		FDR.output<-data.frame(rownames(edata)[lmFitFDR.rob], ff$coefficients, ff$p.value)}
	if(length(lmFitFDR.rob) == 1){
		ff<- summary(glm(t(as.matrix(edata[lmFitFDR.rob,]))~modSv[,-1],family=gaussian))
		FDR.output<-data.frame(rownames(edata)[lmFitFDR.rob], t(ff$coefficients[,1]), t(ff$coefficients[,4]))}
	if(length(lmFitFDR.rob) < 1){
		FDR.output<- t(rep("NA",(ncol(modSv)*2 + 1)))}

	if(length(lmFitBon.rob) > 1){
		bb<-eBayes(lmFit(edata[lmFitBon.rob,],design=modSv,method=linear))
		Bon.output<-data.frame(rownames(edata)[lmFitBon.rob],bb$coefficients, bb$p.value )}
	if(length(lmFitBon.rob) == 1){
		bb<- summary(glm(t(as.matrix(edata[lmFitBon.rob,]))~modSv[,-1],family=gaussian))
		Bon.output<-data.frame(rownames(edata)[lmFitBon.rob], t(bb$coefficients[,1]) , t(bb$coefficients[,4]) )}
	if(length(lmFitBon.rob) < 1){
		Bon.output <- t(rep("NA",(ncol(modSv)*2 + 1)))}
	
	colnames(FDR.output)<-colnames(Bon.output)<-c("Row_Ind_Select","Int_Coeff",rep("Var_Coeff",ncol(mod)-1),
						rep("SV_Coeff",svobj$n.sv),"Int_Pvalue",rep("Var_Pvalue",ncol(mod)-1),
						rep("SV_Pvalue",svobj$n.sv))
	}


output=list(TT.cpg=rownames(edata)[final.temp],train.cpg=train.length,test.cpg=test.length,
			selection=selection,pvalue.matrix=pvalue.matrix,TT.output=TT.output, FDR.output=FDR.output, Bon.output=Bon.output,
			FDR.cpg=rownames(edata)[lmFitFDR.rob],Bon.cpg=rownames(edata)[lmFitBon.rob])

output
}
