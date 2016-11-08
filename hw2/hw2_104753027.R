######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
library("Biostrings",verbose=F,quietly=T)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<- as.integer(args[i+1]) 
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-as.integer(args[i+1])
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# read fasta file
ff <- readAAStringSet(i_f)


seq_name = names(ff)
sequence = paste(ff)
# print(paste("sequences  :", sequence))
# print("=====================================================")
# print(paste("sequences 1 :", sequence[1]))
# print(paste("sequences 2 :", sequence[2]))

# aln length
aln_length <- nchar(sequence[1])

# read score file
s_m<-read.table(s_f)
s_m<-as.matrix(s_m)

# nr <- as.numeric(nchar(sequence[1])+1)
# nc <- as.numeric(nchar(sequence[2])+1)
nc <- nchar(sequence[1])+1
nr <- nchar(sequence[2])+1

# dynamic programing global alignment
# initial matrix
#-----------------------------------------------------------
print(paste("nc : ",nc))
print(paste("nr : ",nr))
mat <- matrix(0, nrow = nr, ncol = nc)
mat_tb <- matrix("0", nrow = nr, ncol = nc)

#print(mat)
mat[1,1]<- 0
cnt <-0
for(i in 2:nc)
{cnt <- cnt-1
 mat[1,i]<-cnt
 mat_tb[1,i] <- "Del"
}
  
cnt<-0
for(i in 2:nr)
{cnt <- cnt-1
 mat[i,1]<-cnt
 mat_tb[i,1] <- "Ins"
}
# 
for(i in 2:nr){
  for(j in 2:nc){
    diago <- 0
    if(substring(sequence[2], i-1, i-1) == substring(sequence[1], j-1, j-1)) {
      diago <- mat[i-1,j-1]+2
    }else{
      diago <- mat[i-1,j-1]-1
    }
    x <-c(mat[i-1,j]-1, mat[i,j-1]-1, diago)
    mat[i,j] <- max(x)
    if(which(x==max(x))==1 ){
      mat_tb[i,j] <- "Ins"
    }else if(which(x==max(x))==2){
      mat_tb[i,j] <- "Del"
    }else{
      mat_tb[i,j] <- "Sub"
    }
  }
}


# print ("------------------------matrix ---------------------")
# print(mat)
# print(mat_tb)
#-----------------------------------------------------------

# matrix trace-back
#-----------------------------------------------------------
# alg1,alg2 are global alignment result
alg1 <- ""
alg2 <- ""
i <- 61
j <- 58

print(mat[0,47])
print(mat_tb[0,47])
while(!(i==1 && j==1)){
  if (mat_tb[i,j]=="Sub")    #SUBSTITUTION
  { i<- i-1
    j<- j-1
    alg1 <- paste(alg1, substring(sequence[1], i, i), sep="")
    alg2 <- paste(alg2, substring(sequence[2], j, j), sep="")
  }
  else if (mat_tb[i,j]=="Ins") 	#DELETION
  { i<- i-1
    j<- j
    alg1 <- paste(alg1, substring(sequence[1], i, i), sep="")
    alg2 <- paste(alg2, "-", sep="")
  }
  else if (mat_tb[i,j]=="Del") 	#INSERTION
  { i<- i
    j<- j-1
    alg1 <- paste(alg1, "-", sep="")
    alg2 <- paste(alg2, substring(sequence[2], j, j), sep="")
  }

  print(paste(i, j, sep=" "))

}

print("here")
# reverse string
alg1 <- paste(rev(substring(alg1,1:nchar(alg1),1:nchar(alg1))),collapse="")
alg2 <- paste(rev(substring(alg2,1:nchar(alg2),1:nchar(alg2))),collapse="")


print("here2")
print(alg1)
print(alg2)

sequence[1] <- alg1
sequence[2] <- alg2

# 
aln_score<-0
for(i in 1:aln_length)
{
  a<-substring(sequence[1], i, i)
  b<-substring(sequence[2], i, i)
  
  if((a != "-")&&(b != "-"))
  {
    print(paste(a, "-", b, "=", s_m[a,b]))
    aln_score = aln_score + s_m[a,b]
  }
  else{
    aln_score = aln_score + g_o
  }
}

print(aln_score)

# output
#write.fasta(as.list(sequence),seq_name,ff)
write.fasta(sequences=sequence,names=seq_name,file.out="result.fasta")
writeXStringSet(ff, o_f)