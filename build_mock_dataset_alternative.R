build_mock_dataset <- function() {
  load("./outputs/current_cell_count100.dat") 
  fiducial <- current_cell_count

  N_Pop <- sum(fiducial[[1]]$totals)
  
  microdata <- list()
  microdata$ordering <- sample(1:N_Pop,N_Pop)
  microdata$sensitivity <- 1:N_Pop
  microdata$error_scale <- (1.0+(1:N_Pop)/N_Pop)/sqrt(5)
  microdata$error_draw <- rnorm(N_Pop,0,microdata$error_scale)                                
  
  N_MB <- geographical_hierarchy[[1]]$N
  H <- length(fiducial[[1]]$classes[1,])
  x_indices <- as.numeric(matrix(1:N_MB,nrow=N_MB,ncol=M))
  y_indices <- as.numeric(t(matrix(1:H,nrow=M,ncol=N_MB)))
  x_indices <- cbind(x_indices,as.numeric(fiducial[[1]]$classes))
  y_indices <- cbind(y_indices,as.numeric(fiducial[[1]]$classes))
  xx <- list()
  for (i in 1:length(x_indices[,1])) {xx[[i]] <- rep(x_indices[i,1],x_indices[i,2])}
  xx <- unlist(xx)
  yy <- list()
  for (i in 1:length(y_indices[,1])) {yy[[i]] <- rep(y_indices[i,1],y_indices[i,2])}
  yy <- unlist(yy)
  xxyy <- cbind(xx,yy)
  microdata$ordering_MB <- xxyy[microdata$ordering,]
  
  N_SA1 <- geographical_hierarchy[[2]]$N
  x_indices <- as.numeric(matrix(1:N_SA1,nrow=N_SA1,ncol=M))
  y_indices <- as.numeric(t(matrix(1:H,nrow=M,ncol=N_SA1)))
  x_indices <- cbind(x_indices,as.numeric(t(do.call(rbind,(lapply(fiducial[[2]]$classes,as.numeric))))))
  y_indices <- cbind(y_indices,as.numeric(t(do.call(rbind,(lapply(fiducial[[2]]$classes,as.numeric))))))
  xx <- list()
  for (i in 1:length(x_indices[,1])) {xx[[i]] <- rep(x_indices[i,1],x_indices[i,2])}
  xx <- unlist(xx)
  yy <- list()
  for (i in 1:length(y_indices[,1])) {yy[[i]] <- rep(y_indices[i,1],y_indices[i,2])}
  yy <- unlist(yy)
  xxyy <- cbind(xx,yy)
  microdata$ordering_SA1 <- xxyy[microdata$ordering,]
  
  N_SA2 <- geographical_hierarchy[[3]]$N
  x_indices <- as.numeric(matrix(1:N_SA2,nrow=N_SA2,ncol=M))
  y_indices <- as.numeric(t(matrix(1:H,nrow=M,ncol=N_SA2)))
  x_indices <- cbind(x_indices,as.numeric(t(do.call(rbind,(lapply(fiducial[[3]]$classes,as.numeric))))))
  y_indices <- cbind(y_indices,as.numeric(t(do.call(rbind,(lapply(fiducial[[3]]$classes,as.numeric))))))
  xx <- list()
  for (i in 1:length(x_indices[,1])) {xx[[i]] <- rep(x_indices[i,1],x_indices[i,2])}
  xx <- unlist(xx)
  yy <- list()
  for (i in 1:length(y_indices[,1])) {yy[[i]] <- rep(y_indices[i,1],y_indices[i,2])}
  yy <- unlist(yy)
  xxyy <- cbind(xx,yy)
  microdata$ordering_SA2 <- xxyy[microdata$ordering,]
  
  MB_aggregate_errors <- fiducial[[1]]$classes*0
  for (i in 1:N_MB) {
    for (j in 1:H) {
      in_bin <- which(microdata$ordering_MB[,1]==i & microdata$ordering_MB[,2]==j)
      sens_in_bin <- microdata$sensitivity[in_bin]
      error_in_bin <- microdata$error_draw[in_bin]
      if (length(in_bin) > 5) {
        MB_aggregate_errors[i,j] <- sum(error_in_bin[sort.list(sens_in_bin,decreasing=TRUE)][1:5])
      } else if (length(in_bin) > 0 & length(in_bin) <=5) {
        MB_aggregate_errors[i,j] <- sum(error_in_bin)
      } else {
        MB_aggregate_errors[i,j] <- 0
      }
    }
  }
  
  perturbed <- fiducial
  perturbed[[1]]$classes <- round(fiducial[[1]]$classes+MB_aggregate_errors)
  perturbed[[1]]$classes <- perturbed[[1]]$classes*(perturbed[[1]]$classes>2)
  
  MB_aggregate_errors <- fiducial[[1]]$totals*0
  for (i in 1:N_MB) {
    in_bin <- which(microdata$ordering_MB[,1]==i)
    sens_in_bin <- microdata$sensitivity[in_bin]
    error_in_bin <- microdata$error_draw[in_bin]
    if (length(in_bin) > 5) {
      MB_aggregate_errors[i] <- sum(error_in_bin[sort.list(sens_in_bin)][1:5])
    } else if (length(in_bin) > 0 & length(in_bin) <=5) {
      MB_aggregate_errors[i] <- sum(error_in_bin)
    } else {
      MB_aggregate_errors[i] <- 0
    }
  }
  
  perturbed[[1]]$totals <- round(fiducial[[1]]$totals+MB_aggregate_errors)
  perturbed[[1]]$totals <- perturbed[[1]]$totals*(perturbed[[1]]$totals>2)
  
  SA1_aggregate_errors <- t(do.call(rbind,(lapply(fiducial[[2]]$classes,as.numeric))))*0
  for (i in 1:N_SA1) {
    for (j in 1:H) {
      in_bin <- which(microdata$ordering_SA1[,1]==i & microdata$ordering_SA1[,2]==j)
      sens_in_bin <- microdata$sensitivity[in_bin]
      error_in_bin <- microdata$error_draw[in_bin]
      if (length(in_bin) > 5) {
        SA1_aggregate_errors[i,j] <- sum(error_in_bin[sort.list(sens_in_bin)][1:5])
      } else if (length(in_bin) > 0 & length(in_bin) <=5) {
        SA1_aggregate_errors[i,j] <- sum(error_in_bin)
      } else {
        SA1_aggregate_errors[i,j] <- 0
      }
    }
  }
  
  perturbed[[2]]$classes <- round(fiducial[[2]]$classes+SA1_aggregate_errors)
  perturbed[[2]]$classes <- perturbed[[2]]$classes*(perturbed[[2]]$classes>2)
  
  SA1_aggregate_errors <- as.numeric(fiducial[[2]]$totals)*0
  for (i in 1:N_SA1) {
    in_bin <- which(microdata$ordering_SA1[,1]==i)
    sens_in_bin <- microdata$sensitivity[in_bin]
    error_in_bin <- microdata$error_draw[in_bin]
    if (length(in_bin) > 5) {
      SA1_aggregate_errors[i] <- sum(error_in_bin[sort.list(sens_in_bin)][1:5])
    } else if (length(in_bin) > 0 & length(in_bin) <=5) {
      SA1_aggregate_errors[i] <- sum(error_in_bin)
    } else {
      SA1_aggregate_errors[i] <- 0
    }
  }
  
  perturbed[[2]]$totals <- round(fiducial[[2]]$totals+SA1_aggregate_errors)
  perturbed[[2]]$totals <- perturbed[[2]]$totals*(perturbed[[2]]$totals>2)
  
  
  SA2_aggregate_errors <- t(do.call(rbind,(lapply(fiducial[[3]]$classes,as.numeric))))*0
  for (i in 1:N_SA2) {
    for (j in 1:H) {
      in_bin <- which(microdata$ordering_SA2[,1]==i & microdata$ordering_SA2[,2]==j)
      sens_in_bin <- microdata$sensitivity[in_bin]
      error_in_bin <- microdata$error_draw[in_bin]
      if (length(in_bin) > 5) {
        SA2_aggregate_errors[i,j] <- sum(error_in_bin[sort.list(sens_in_bin)][1:5])
      } else if (length(in_bin) > 0 & length(in_bin) <=5) {
        SA2_aggregate_errors[i,j] <- sum(error_in_bin)
      } else {
        SA2_aggregate_errors[i,j] <- 0
      }
    }
  }
  
  perturbed[[3]]$classes <- round(fiducial[[3]]$classes+SA2_aggregate_errors)
  perturbed[[3]]$classes <- perturbed[[3]]$classes*(perturbed[[3]]$classes>2)
  
  SA2_aggregate_errors <- as.numeric(fiducial[[3]]$totals)*0
  for (i in 1:N_SA2) {
    in_bin <- which(microdata$ordering_SA2[,1]==i)
    sens_in_bin <- microdata$sensitivity[in_bin]
    error_in_bin <- microdata$error_draw[in_bin]
    if (length(in_bin) > 5) {
      SA2_aggregate_errors[i] <- sum(error_in_bin[sort.list(sens_in_bin)][1:5])
    } else if (length(in_bin) > 0 & length(in_bin) <=5) {
      SA2_aggregate_errors[i] <- sum(error_in_bin)
    } else {
      SA2_aggregate_errors[i] <- 0
    }
  }
  
  perturbed[[3]]$totals <- round(fiducial[[3]]$totals+SA2_aggregate_errors)
  perturbed[[3]]$totals <- perturbed[[3]]$totals*(perturbed[[3]]$totals>2)
  
  return(list('perturbed'=perturbed,'fiducial'=fiducial))
}
