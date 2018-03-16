library(readr)

# data <- read_csv("Stevens_stats.txt", col_names = FALSE)
# gini = as.numeric(data[1,])
# reads = as.numeric(data[2,])
# cistrans = as.numeric(data[3,])
# ploidy = as.numeric(rep(0,length(data[1,])))

gini = c(0.515,0.475,0.401,0.503,0.473,0.433,0.478,0.427)
reads = c(156839,116149,52979,135300,105870,106416,13250,61580)
# cistrans = c(82.506,89.243,68.135,89.041,84.780,90.851,69.840,86.051)
# ploidy = rep(0,8)

data <- read_csv("ND_stats.csv", col_names = FALSE)
gini = c(gini, as.numeric(data[1,]))
reads = c(reads, as.numeric(data[2,]))
cistrans = c(cistrans, as.numeric(data[3,]))
ploidy = c(ploidy, as.numeric(rep(1,length(data[1,]))))

data <- read_csv("NH_stats.csv", col_names = FALSE)
gini = c(gini, as.numeric(data[1,]))
reads = c(reads, as.numeric(data[2,]))
cistrans = c(cistrans, as.numeric(data[3,]))
ploidy = c(ploidy, as.numeric(rep(0,length(data[1,]))))

data <- read_csv("R_stats.csv", col_names = FALSE)
gini = c(gini, as.numeric(data[1,]))
reads = c(reads, as.numeric(data[2,]))
cistrans = c(cistrans, as.numeric(data[3,]))
ploidy = c(ploidy, as.numeric(rep(1,length(data[1,]))))

aty <- axTicks(2)
labels <- sapply(aty,function(i)
  as.expression(bquote(10^ .(i)))
)
axis(2,at=aty,labels=labels)

par(mar=c(5.1, 4.5, 1.5, 1.5))
plot(reads, gini, pch=20, cex=0.6, axes=FALSE, log="x",
     cex.lab=1.5, xlab="Number of contacts", ylab="GiniQC")
labels <- sapply(c(2,3,4,5,6),function(i)
  as.expression(bquote(10^ .(i)))
)
axis(1, at=c(1e2, 1e3, 1e4, 1e5, 1e6), labels=labels, cex.axis=1.5)
axis(2, cex.axis=1.5)
box()
cor.test(log(reads),gini)
plot(cistrans, gini, pch=20, cex=0.6, cex.lab=1.5, cex.axis=1.5,
     xlab=expression(Percentage~of~contacts~'in'~italic(cis)), ylab="GiniQC", xlim=c(0,100))
cor.test(cistrans,gini)

haploid_gini = vector(mode="numeric")
diploid_gini = vector(mode="numeric")
for(i in 1:length(ploidy)){
  if(ploidy[i]==0)
    haploid_gini = c(haploid_gini, gini[i])
  else
    diploid_gini = c(diploid_gini, gini[i])
}
par(mar=c(2.5, 4.5, 1.5, 1.5))
boxplot(haploid_gini, diploid_gini, names=c("Haploid cells", "Diploid cells"),
        ylab="GiniQC", cex.lab=1.5, cex.axis=1.5)
cor.test(ploidy, gini)

cor.test(ploidy, reads)

glm(formula = gini ~ ploidy + log(reads) + cistrans)

