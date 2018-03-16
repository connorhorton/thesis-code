cell1_contacts = c(100000,10000,110000,120000,130000,140000,150000,20000,30000,40000,50000,60000,70000,80000,90000)
cell1_gini = c(0.468,0.513,0.470,0.487,0.495,0.507,0.511,0.425,0.414,0.413,0.422,0.464,0.454,0.447,0.450)
cell2_contacts = c(100000,10000,110000,20000,30000,40000,50000,60000,70000,80000,90000)
cell2_gini = c(0.463,0.478,0.470,0.432,0.431,0.407,0.401,0.419,0.441,0.443,0.448)
cell3_contacts = c(10000,20000,30000,40000,50000)
cell3_gini = c(0.447,0.401,0.391,0.390,0.397)
cell4_contacts = c(100000,10000,110000,120000,130000,20000,30000,40000,50000,60000,70000,80000,90000)
cell4_gini = c(0.482,0.490,0.484,0.488,0.502,0.490,0.464,0.434,0.457,0.443,0.465,0.467,0.456)
cell5_contacts = c(100000,10000,20000,30000,40000,50000,60000,70000,80000,90000)
cell5_gini = c(0.470,0.436,0.424,0.401,0.407,0.410,0.433,0.437,0.437,0.455)
cell6_contacts = c(100000,10000,20000,30000,40000,50000,60000,70000,80000,90000)
cell6_gini = c(0.431,0.426,0.370,0.353,0.388,0.369,0.382,0.391,0.399,0.411)
cell8_contacts = c(10000,20000,30000,40000,50000,60000)
cell8_gini = c(0.462,0.405,0.389,0.386,0.418,0.426)

# final_gini = c(0.515,0.475,0.401,0.503,0.473,0.433,0.478,0.427)
# final_contacts = c(156839,116149,52979,135300,105870,106416,13250,61580)

all_contacts = c(cell1_contacts, cell2_contacts, cell3_contacts, cell4_contacts,
                 cell5_contacts, cell6_contacts, cell8_contacts)
all_gini = c(cell1_gini, cell2_gini, cell3_gini, cell4_gini, cell5_gini, cell6_gini, cell8_gini)


labels=expression(2%*%10^4, 6%*%10^4, 1%*%10^5, 1.4%*%10^5)
plot(cell1_contacts, cell1_gini, col="red", pch=19, axes=FALSE,
     xlim=c(min(all_contacts), max(all_contacts)), ylim=c(min(all_gini),max(all_gini)),
     xlab="Number of contacts", ylab="GiniQC", cex.lab=1.5, cex.axis=1.5)
axis(1, at=c(2e4, 6e4, 1e5, 1.4e5), labels=labels, cex.axis=1.5)
axis(2, cex.axis=1.5)
box()

points(cell2_contacts, cell2_gini, col="orange", pch=19)
points(cell3_contacts, cell3_gini, col="green", pch=19)
points(cell4_contacts, cell4_gini, col="cyan", pch=19)
points(cell5_contacts, cell5_gini, col="blue", pch=19)
points(cell6_contacts, cell6_gini, col="violet", pch=19)
points(cell8_contacts, cell8_gini, col="brown", pch=19)
legend("bottomright", "(x,y)", legend=c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5", "Cell 6", "Cell 8"),
       fill=c("red", "orange", "green", "cyan", "blue", "violet", "brown"), cex=1.5)
cor.test(all_contacts, all_gini)


# points(final_contacts, final_gini, pch=5)

#add the actual values for comparison? yes. that would be good.