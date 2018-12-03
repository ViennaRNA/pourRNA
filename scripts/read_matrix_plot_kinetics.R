require(Matrix)

#rate_matrix_file="/home/quin/Desktop/btest/parkin/test_sparse.bin"
args=commandArgs(trailingOnly=TRUE)
fileNameMatrix=args[1]
startIndex=args[2]
startIndex=as.double(startIndex)

#fileNameStructures=args[1]
#fileNameMatrix=args[2]



read_sparse_binary_rate_matrix <- function(rate_matrix_file) {
  bytesToRead <- file.info(rate_matrix_file)$size;
  #print("toread")
  #print(bytesToRead)
  newdata = file(rate_matrix_file, "rb")
  #varnames = readBin(newdata, character(), n=3)
  n_states = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
  bytesToRead = bytesToRead - 4
  #print(n_states)
  rate_matrix=Matrix(0, nrow=n_states, ncol=n_states, sparse=TRUE)
  while(bytesToRead > 0){
    state_from = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
    #print(state_from)
    bytesToRead = bytesToRead - 4
    number_of_transitions = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
    #print(number_of_transitions)
    bytesToRead = bytesToRead - 4
    for(j in 1:number_of_transitions){
      state_to = readBin(newdata, what="integer", size = 4, n = 1, endian = "little")
      bytesToRead = bytesToRead - 4
      rate_from_to = readBin(newdata, what="double" , size = 8, n = 1, endian = "little")
      bytesToRead = bytesToRead - 8
      #print(state_from)
      #print(state_to)
      #print(rate_from_to)
      rate_matrix[state_from+1, state_to+1] = rate_from_to
    }
  }
  close(newdata)
  #output=list("structures" = varnames, "matrix" = a)
  #output$structures
  return (rate_matrix)
}

rate_matrix = read_sparse_binary_rate_matrix(fileNameMatrix)
#rate_matrix

#structures=read.table(fileNameStructures, stringsAsFactors=FALSE, skip=1)
#structures=structures[,2]
#structures

#rate_matrix=as.matrix(read.table(fileNameMatrix)) # barriers text matrix

n=nrow(rate_matrix)
Q = as.matrix(rate_matrix)
#Q=t(Q)
#typeof(Q[0,0])
#compute the diagonal
for(i in 1:n){
  Q[i,i]=0
  Q[i,i]=-sum(Q[i,])
}

ev=eigen(Q)
eValues=ev$val
eValues
eVectorsL=ev$vec
eValues
eVectorsL

eVectorsR=t(solve(ev$vec))
eVectorsR

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

#Pt
pdf(paste0(fileNameMatrix,".pdf"))
matplot(plotTime,t(Pt),type="l",xlim=c(min(plotTime),max(plotTime)), log="x", lty=1)
#legend("topright",structures,col=1:n,lty=1)
dev.off()





