#!/usr/bin/env Rscript
require(Matrix)
library(rARPACK)
library("optparse")

option_list = list(
  make_option(c("-f", "--sparse_binary_matrix"), type="character", default=NULL,
              help="sparse matrix format (from pourRNA)", metavar="character"),
  make_option(c("-g", "--binary_matrix"), type="character", default=NULL,
              help="binary matrix format (from barriers or pourRNA)", metavar="character"),
  make_option(c("-c", "--text_matrix"), type="character", default=NULL,
              help="binary matrix format (from barriers or pourRNA)", metavar="character"),
  make_option(c("-s", "--states_file"), type="character", default=NULL,
              help="barriers-like output as input file."),
  make_option(c("-i","--initial_state"), type="integer", default=NULL,
              help="state which is populated at the beginning of the Markov process.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sparse_binary_matrix) && is.null(opt$binary_matrix) && is.null(opt$text_matrix)){
  print_help(opt_parser)
  stop("Specify at least the input file!", call.=FALSE)
}

if(is.null(opt$initial_state)){
  print_help(opt_parser)
  stop("Specify at least the start state!", call.=FALSE)
}

startIndex=opt$initial_state

startIndex=as.double(startIndex)

fileNameStructures=NULL
if(!is.null(opt$states_file)){
  fileNameStructures = opt$states_file
}

read_sparse_binary_rate_matrix <- function(rate_matrix_file) {
  bytesToRead <- file.info(rate_matrix_file)$size;
  newdata = file(rate_matrix_file, "rb")
  n_states = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
  bytesToRead = bytesToRead - 4
  rate_matrix=Matrix(0, nrow=n_states, ncol=n_states, sparse=TRUE)
  while(bytesToRead > 0){
    state_from = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
    bytesToRead = bytesToRead - 4
    number_of_transitions = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
    bytesToRead = bytesToRead - 4
    for(j in 1:number_of_transitions){
      state_to = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
      bytesToRead = bytesToRead - 4
      rate_from_to = readBin(newdata, what="double" , size = 8, n = 1, endian = "little")
      bytesToRead = bytesToRead - 8
      rate_matrix[state_from+1, state_to+1] = rate_from_to
    }
  }
  close(newdata)
  return (rate_matrix)
}

read_binary_rate_matrix <- function(rate_matrix_file) {
  bytesToRead <- file.info(rate_matrix_file)$size;
  newdata = file(rate_matrix_file, "rb")
  n_states = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
  rate_matrix=Matrix(0, nrow=n_states, ncol=n_states, sparse=TRUE)
  for(i in 1:n_states){
    for(j in 1:n_states){
      rate_from_to = readBin(newdata, what="double" , size = 8, n = 1, endian = "little")
      if(rate_from_to != 0.0){
        rate_matrix[i, j] = rate_from_to
      }
    }
  }
  close(newdata)
  return (rate_matrix)
}

rate_matrix = NULL
fileNameMatrix = NULL
if (!is.null(opt$sparse_binary_matrix)){
  fileNameMatrix=opt$sparse_binary_matrix
  rate_matrix = read_sparse_binary_rate_matrix(fileNameMatrix)
}
if (!is.null(opt$binary_matrix) && is.null(fileNameMatrix)){
  fileNameMatrix=opt$binary_matrix
  rate_matrix = read_binary_rate_matrix(fileNameMatrix)
}
if(!is.null(opt$text_matrix) && is.null(fileNameMatrix)){
  fileNameMatrix=opt$text_matrix
  rate_matrix = as.matrix(read.table(fileNameMatrix)) # barriers text matrix
}

print("read matrix")

if(!is.null(fileNameStructures)){
  structures=read.table(fileNameStructures, stringsAsFactors=FALSE, skip=1)
  structures=structures[,2]
}

n=nrow(rate_matrix)
#Q = as.matrix(rate_matrix)
Q = rate_matrix
#Q=t(Q)
#typeof(Q[0,0])
#compute the diagonal
for(i in 1:n){
  Q[i,i]=0
  Q[i,i]=-sum(Q[i,])
}

n_eigenvalues = n
str(le <- eigs(Q, k = n_eigenvalues, which = "LM", opts = list(retvec = TRUE), maxitr=1000, tol=1e-10)) # or dsyMatrix
eValues = le$values
#eValues
eVectorsL = le$vectors
#eVectorsL #apply(eVectorsL, 2, as.numeric)
#eVectorsL

eVectorsR = t(solve(eVectorsL))
#eVectorsR

#Q = as.matrix(Q)
#ev=eigen(Q)
#eValues=ev$val
#eValues
#eVectorsL=ev$vec
#eValues
#eVectorsL

print("eigenvalues computed")
#eVectorsR=t(solve(ev$vec))

#equilibr. density
#which.max(eValues)
#eVectorsL[which.max(eValues),]

print("eigenvectors computed")

## initial distribution
p0 <- rep(0, 1, n)
p0[startIndex]=1
time <- c(0.001 %o% 10^(seq(0,20,0.1)))  #log time steps

V=t(eVectorsR)

Pt=c()
plotTime=c()
startP0=as.vector(p0)
for(i in time){  
  Rt = (eVectorsL %*% diag(exp(eValues*i)) %*% V)
  t1 <- startP0 %*% Rt 
  if(is.nan(sum(t1)) || abs(sum(t1)) < 0.90 || abs(sum(t1)) > 1.10){
    print("break!")
    print(i)
    break
  }
  if(length(Pt) == 0){
    Pt = as.vector(t1)
  } else {
    Pt=cbind(Pt,as.vector(t1))
  }
  plotTime=rbind(plotTime,i)
}

Pt_t = t(Pt)
colnames(Pt_t) <- structures
print(structures)
#Pt
pdf(paste0(fileNameMatrix,".pdf"))
#col_set <- rainbow(10) #palette('default')
col_set = c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray", "green", "orange")
#palette('default')[1:8]
#print(col_set)

all_colors <- rep("black", 1, n)
if(!is.null(fileNameStructures)){
  
  sums <- rep(0, 1, n)
  selected_structures = rep(0, 1, 10)
  best_colors = rep(0, 1, 10)
  for(i in 1:n){
    sums[i] = sum(Pt[i,])
  }
  print(sums)
  tmp_sums = sums
  best_indices <- rep(0, 1, 10)
  
  for(i in 1:10){
    #v_max = max(tmp_sums)
    i_max = which.max(tmp_sums) #match(v_max, sums)
    tmp_sums = tmp_sums[-i_max]
    
    tmp_i_max = i_max
    to = i-1
    if(to>0){
    for(j in 1:to){
      if(best_indices[j] < i + i_max){
        tmp_i_max = tmp_i_max + 1
      }
    }
    }
    i_max = tmp_i_max
    #print(i_max)
    best_indices[i] = i_max
    best_colors[i] = col_set[i]
    selected_structures[i] = structures[best_indices[i]]
    all_colors[i_max] = col_set[i]
  }  
  print(selected_structures)
  
  #legend("topright",structures,lty=1, bg='white', , col = col_set)
}
print(all_colors)
matplot(plotTime,Pt_t,type="l",xlim=c(min(plotTime),max(plotTime)), log="x", lty=1, col = all_colors)
legend("topright",selected_structures,lty=1, bg='white', col = best_colors)

dev.off()





