# mergeBlast.r 
# @authors: Isavannah Reyes & Dylan Sosa
# 2 May 2017

setwd("~/Documents/Semester 4.2/Bioinformatics/Lab/parsedBlastData")
filenames = list.files(path = "~/Documents/Semester 4.2/Bioinformatics/Lab/parsedBlastData", pattern = "*csv")
# set your own wd dude/filenames path

c = read.csv(filenames[1],header=T)
c = c[order(c$query),]


# create matrix of all blast output 
library("rowr")
for (i in 2:length(filenames)){
    p <- read.csv(filenames[i], header = T, stringsAsFactors = F) # parsed blast output
    p <- p[order(p$query),]
    c <- cbind.fill(c,p,fill = NA) # concatenated blasts
}

# remove unnecessary columns 
toMatch <- c("query.","hit","X")
toRemove <- grep(paste(toMatch,collapse="|"),colnames(c),value = TRUE)
#colnames(c)
c_final = c[,!(names(c) %in% toRemove)]

# change column names 
library("rlist")
filenames <- list.prepend(filenames,"Cre")
filenames <- lapply(filenames, substr, start = 1 ,stop = 12)
colnames(c_final) <- filenames

# convert factors to numeric
#c_final <- as.numeric(as.character(c_final))

row.names(c_final) <- NULL

c3<-data.frame(c_final, row.names=NULL)
row.names(c3) <- c()
row.names(c_final)

c4 <- as.matrix(c_final)
rownames(c4) <- rep("", nrow(c4))

x <- data.frame(a=1:5, b= 1:5)
x2 <- as.matrix(x)
rownames(x2) <- rep("", nrow(x2))

# change NaN and NA to -10
c_final[is.na(c_final)] <- -1
c_final[is.nan(c_final)] <- -1

gz1 <- gzfile("parsedBlast.csv.gz", "w")
write.csv(c_final, gz1)
close(gz1)
