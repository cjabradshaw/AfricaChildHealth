aicW <- function(model.list, finite = TRUE, null.model = NULL, N, order = TRUE){
	if(hasArg(N)) { N <- N }
        else N <- NULL
	if(is.null(null.model)) {
		null.model <- model.list[[length(model.list)]]
	}
	LLlist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			LLlist <- c(LLlist, as.numeric(logLik(model.list[[index]])))
		}
		else {LLlist <- c(LLlist, as.numeric(logLik(model.list[[index]])))}
	}
	aiclist <- vector()
	for(index in 1:length(model.list)){
		aicname <- "AIC"
		if(finite==T) {
			aiclist <- c(aiclist, AIC.finite.lm(model.list[[index]], N))
			aicname <- "AICc"}
		else {aiclist <- c(aiclist, AIC.lm(model.list[[index]]))}
	}
	biclist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			biclist <- c(biclist, BIC.lm(model.list[[index]], N))
		}
		else {biclist <- c(biclist, BIC.lm(model.list[[index]], N))}
	}
	dflist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			dflist <- c(dflist, as.numeric(attributes(logLik(model.list[[index]]))$df))
		}
		else {dflist <- c(dflist, length(coef(model.list[[index]])))}
	}
	delist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			delist <- c(delist, lmer.de(model.list[[index]], null.model, N))
		}
		else {delist <- c(delist, lmer.de(model.list[[index]], null.model, N))}
	}
	modform <- vector()	
	if(is.null(names(model.list))){
		for(index in 1:length(model.list)){
			modform <- c(modform, formula(model.list[[index]]))
		}
	}	
	if(!is.null(names(model.list))){
		models <- names(model.list)
		delta.i <- aiclist - min(aiclist)
		aicw <- exp(delta.i * -0.5)/sum(exp(delta.i * -0.5))
		delta.BIC <- biclist - min(biclist)
    bicw <- exp(delta.BIC * -0.5)/sum(exp(delta.BIC * -0.5))
		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(biclist,3), round(delta.BIC, 3), round(bicw,4), round(delist, 2)), ncol = 9)
		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "wAICc", "BIC", "dBIC", "wBIC", "%DE")
		# rownames(return.matrix) <- sapply(strsplit(models, NULL), function(x) paste(x, collapse = " + "))
		rownames(return.matrix) <- models
	}
	else { 
		delta.i <- aiclist - min(aiclist)
		aicw <- exp(delta.i * -0.5)/sum(exp(delta.i * -0.5))
		delta.BIC <- biclist - min(biclist)
    bicw <- exp(delta.BIC * -0.5)/sum(exp(delta.BIC * -0.5))
		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(biclist,3), round(delta.BIC, 3), round(bicw,4), round(delist, 2)), ncol = 9)
#		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(delta.BIC, 3), round(delist, 2)), ncol = 7)
		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "wAICc", "BIC", "dBIC", "wBIC", "%DE")
#		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "weight", "dBIC", "%DE")
		rownames(return.matrix) <- modform
	}
	if(order) {
	return.matrix <- return.matrix[order(return.matrix[, 4]), ]
	}
	return(return.matrix)
}

AIC.finite.lm <- function(model, N){
	LL = as.numeric(logLik(model))
	p = as.numeric(attributes(logLik(model))$df)
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
	AICc <- AIC(logLik(model)) + (2 * p * (p + 1))/(N - p - 1)
	return(AICc)
}

AIC.lm <- function(object) {
	ret <- AIC(logLik(object))
	return(ret)
}

BIC.lm <- function(model, N) {
	LL = as.numeric(logLik(model))
	p = as.numeric(attributes(logLik(model))$df)
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
	BIC <- (-2 * LL) + (p * log(N))
	return(BIC)
}

lmer.de <- function(model, null.model, N) {
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
    if(class(null.model)[1] == "mer"){
        null.dev <- as.vector(null.model@deviance["ML"])
	    model.lr <- as.vector(null.model@deviance["ML"]) - as.vector(model@deviance["ML"])
    }
    else {
		if(class(model)[1] == "mer"){
			null.dev <- deviance(null.model)
		    model.lr <- (deviance(null.model)) - as.vector(model@deviance["ML"])
		}
		else {
			null.dev <- deviance(null.model)
		    model.lr <- deviance(null.model) - deviance(model)
		}
	}
	de <- 100 * (model.lr/null.dev)
    return(de)
}
