# Rscript code.r filename errorth

args <- commandArgs(trailingOnly = TRUE)


print(paste("Plotting diversity for file ",args[1],sep=""))
data = read.csv(args[1],header=TRUE)
errorth = as.numeric(args[2])

#seqs=unique(as.vector(data$seq))

jpeg(filename = paste(substr(args[1],1,nchar(args[1])-4),"_errorth_",args[2],".jpeg",sep=""), width =1000, height = 400,quality=100,units = "px", pointsize = 12, bg="white",res = NA) 

datatoplotFirst = data[which(data$First>errorth & data$First<(1-errorth)),]
plot(as.numeric(rownames(datatoplotFirst)),datatoplotFirst$First,pch=20,xlim=c(0,nrow(data)+300),ylim=c(-0.2,1.1),xlab="Position in the genome",ylab="Proportion of each nucleotide",main="Proportion of nucleotide diversity.",xaxt='n',yaxt='n')


#data0 <- data[which(is.element(data$pos,c(1,500,1000,2000))),]
#tickes0 <- rownames(data0) 
#labs0 <- data0$pos
axis(side = 1) #, at = tickes0, labels= labs0)
axis(side = 2, at = c(0,0.25,0.5,0.75,1), labels= c("0","0.25","0.5","0.75","1"))

#vlines <- c(rownames(data[which(data$pos==1),]),rownames(data[nrow(data),]))
#for(vline in vlines) abline(v=vline,col="blue",lty=2)

datatoplotSecond = data[which(data$Second>errorth),]
points(as.numeric(rownames(datatoplotSecond)),datatoplotSecond$Second,col="red",pch=20)
datatoplotThird = data[which(data$Third>errorth),]
points(as.numeric(rownames(datatoplotThird)),datatoplotThird$Third,col="blue",pch=20)
datatoplotFour = data[which(data$Four>errorth),]
points(as.numeric(rownames(datatoplotFour)),datatoplotFour$Four,col="green",pch=20)

#for(i in 1:length(seqs)) text(tickes[i],1.08,strsplit(seqs[i],"_")[[1]][2]) 

legend(min(as.numeric(rownames(data))),0.7, pch=20,c("First call","Second call","Third call","Four call"), col=c("black","red","blue","green")) #-700

datatoplotCov = data[which(data$cov>0),]

#points(as.numeric(rownames(datatoplotCov)),rep(-0.2,nrow(datatoplotCov)),col="blue",pch=".")
m=(-0.999)/(-max(data$cov))
points(as.numeric(rownames(data)),m*(data$cov),col="blue",pch=".")

axis(side = 4, at = c(0,0.25,0.5,0.75,1), labels= as.character(max(data$cov)*c(0,0.25,0.5,0.75,1)))
mtext("Coverage", side=4, line=3)
#mtext("Coverage", side = 4, las=0)
#, line = 3) 
#, cex = par("cex.lab"))

abline(h=0,lty=2)
abline(h=1,lty=2)



dev.off()
