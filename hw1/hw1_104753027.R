# read PAM1 from data
pam1<-read.table("C:/Users/user/Desktop/hw1/pam1.txt")


# check PAM1 data
dim(pam1)
str(pam1)
#---------------

base<-pam1/10000

typeof(base)
#basematrix <- do.call(rbind, base)
basematrix <- data.matrix(base)
#basematrix
dim(basematrix)
typeof(basematrix)


# Read in the data
#x <- scan("./pam1.txt", what="", sep="\n")
# Separate elements by one or more whitepace
#y <- strsplit(x, "[[:space:]]+")

#dim(x)
#y

pam250 <- basematrix
for (i in 1:249) {
  pam250 <- pam250%*%basematrix
}

#pam250
pam250_out <- round(pam250, digits = 2)*100
pam250_out

write.table(pam250_out, file = "C:/Users/user/Desktop/hw1/pam250.txt", sep = "\t")

