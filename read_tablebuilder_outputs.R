### Read in TableBuilder outputs for the geographical hierarchy

tablebuilder_outputs <- list() # contains the TableBuilder files matching the given geographical hierarchy

## <<< Perth example
filename_MB <- "MB_by_HCFMD.csv"
filename_SA1 <- "SA1_by_HCFMD.csv"
filename_SA2 <- "SA2_by_HCFMD.csv"
# MB level Tablebuilder outputs
xdata <- read.csv(filename_MB,skip=11,header = FALSE)
raw_outputs <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("raw_outputs <- cbind(raw_outputs,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(filename_MB,skip=9,header = FALSE,nrows = 1)
xcolnames <- unlist(c("MB",ydata[-1][1:(length(raw_outputs[1,])-1)]))
colnames(raw_outputs) <- xcolnames
raw_outputs <- as.data.frame(raw_outputs)
raw_outputs <- raw_outputs[-which(is.na(raw_outputs$MB)),]
raw_outputs <- raw_outputs[match(as.numeric(MB_polygons$MB_CODE21),raw_outputs$MB),]
raw_outputs <- as.matrix(raw_outputs)
classes <- raw_outputs[,2:(dim(raw_outputs)[2]-1)]
totals <- raw_outputs[,(dim(raw_outputs)[2])]
tablebuilder_outputs[[1]] <- list('classes'=classes,'totals'=totals)
# SA1 level Tablebuilder outputs
xdata <- read.csv(filename_SA1,skip=11,header = FALSE)
raw_outputs <- as.numeric(xdata$V1)
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("raw_outputs <- cbind(raw_outputs,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(filename_SA1,skip=9,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA1",ydata[-1][1:(length(raw_outputs[1,])-1)]))
colnames(raw_outputs) <- xcolnames
raw_outputs <- as.data.frame(raw_outputs)
raw_outputs <- raw_outputs[-which(is.na(raw_outputs$SA1)),]
raw_outputs <- raw_outputs[match(as.numeric(SA1_polygons$SA1_CODE21),raw_outputs$SA1),]
raw_outputs <- as.matrix(raw_outputs)
classes <- raw_outputs[,2:(dim(raw_outputs)[2]-1)]
totals <- raw_outputs[,(dim(raw_outputs)[2])]
tablebuilder_outputs[[2]] <- list('classes'=classes,'totals'=totals)
# SA2 level Tablebuilder outputs
xdata <- read.csv(filename_SA2,skip=11,header = FALSE)
raw_outputs <- as.character(xdata$V1) #SA2 units come as NAMES!
for (i in 2:(dim(xdata)[2]-1)) {
  eval(parse(text=paste0("raw_outputs <- cbind(raw_outputs,as.numeric(xdata$V",i,"))")))
}
ydata <- read.csv(filename_SA2,skip=9,header = FALSE,nrows = 1)
xcolnames <- unlist(c("SA2",ydata[-1][1:(length(raw_outputs[1,])-1)]))
colnames(raw_outputs) <- xcolnames
raw_outputs <- as.data.frame(raw_outputs)
raw_outputs <- raw_outputs[match(SA2_polygons$SA2_NAME21,raw_outputs$SA2),]
raw_outputs$SA2 <- as.numeric(SA2_polygons$SA2_CODE21)
raw_outputs <- matrix(as.numeric(unlist(raw_outputs)),nrow=nrow(raw_outputs),ncol=ncol(raw_outputs))
classes <- raw_outputs[,2:(dim(raw_outputs)[2]-1)]
totals <- raw_outputs[,(dim(raw_outputs)[2])]
tablebuilder_outputs[[3]] <- list('classes'=classes,'totals'=totals)
## Perth example >>>

K <- dim(tablebuilder_outputs[[1]]$classes)[2] # number of categories in TableBuilder query
