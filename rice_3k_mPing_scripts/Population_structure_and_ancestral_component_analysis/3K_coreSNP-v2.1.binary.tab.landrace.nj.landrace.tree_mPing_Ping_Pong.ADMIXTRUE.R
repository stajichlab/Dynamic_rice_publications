library("plotrix")
library("ape")
#tree = read.tree(file="3K_coreSNP-v2.1.binary.tab.landrace.nj.tree")
#tree$edge.length[tree$edge.length<0]<-0
sorted_strains <- read.table("base_filtered_v0.7.pruneddata_1M_20kb_0.8.5.Q.sorted.merge.txt", sep="\t")
sorted_list <- as.vector(sorted_strains$V1)

x = read.table(file="rice_line_ALL_3000.anno.landrace.list", sep='\t', header=1)
y1 = setNames(x[,7], x[,1])
y1 = y1[match(gsub("'", '', sorted_list), names(y1))]
y2 = setNames(x[,8], x[,1])
y2 = y2[match(gsub("'", '', sorted_list), names(y2))]
y3 = setNames(x[,9], x[,1])
y3 = y3[match(gsub("'", '', sorted_list), names(y3))]

#sample_colors = setNames(x[,2], x[,1])
#sample_colors = sample_colors[match(gsub("'", '', tree$tip.label), names(sample_colors))]
#sample_colors = as.vector(sample_colors)

pdf("3K_coreSNP-v2.1.binary.tab.landrace.nj.landrace.tree_mPing_Ping_Pong.ADMIXTRUE.pdf", width=7, height=9)
layout(matrix(c(1,2,3,4),1,4),c(0.4, 0.2, 0.2, 0.3))
par(mar=c(4,1,2,4))
#edge_colors=rep("black", length(tree$edge[,2]))

#https://ecomorph.wordpress.com/2014/10/09/phylogenetic-trees-in-r-4/
#edge_num includes all the edge of internal edge or termial edge.
#the latter is what we need.
#edge_num = tree$edge[,2]
#for (i in 1:length(edge_num)){
#     if (edge_num[i] <= length(sample_colors)){
#         if (!is.na(sample_colors[edge_num[i]])){
#             edge_colors[i] = sample_colors[edge_num[i]]
#         }
#     }
#}

##plot tree##########################################################################################################
#plot(tree, edge.color=edge_colors, show.tip.label = FALSE, cex = 0.5, underscore=TRUE)
#leg_inf = cbind(as.vector(x[,2]), as.vector(x[,6]))
#leg_inf = unique(leg_inf)
#leg_inf = leg_inf[order(leg_inf[,2]),]
#xrange  = par("xaxp")
#yrange  = par("yaxp")
#par(xpd=TRUE) #set this legend can be plot into margin area
#legend(x=xrange[2]*0.8, y=yrange[2]*0.99, substr(leg_inf[,2], 1, 4), fill=leg_inf[,1], border=FALSE, bty='n')
#legend(x=xrange[2]*0.6, y=yrange[2]*0.99, leg_inf[,2], fill=leg_inf[,1], border=FALSE, bty='n', cex=1.4)
##plot tree#########################################################################################################

##plot ADMIXTRUE##########################################################################################################
#tbl <- read.table("core_v0.7.pruneddata3.7.Q.reordered.txt")
#ord = tbl[order(tbl$V1,tbl$V2,tbl$V3,tbl$V4,tbl$V5,tbl$V6,tbl$V7,tbl$V8),]
#barplot(t(as.matrix(tbl)), 
#              space = c(0.2),
#              col=rainbow(7),
#              xlab="Individual #", 
#              ylab="Ancestry",
#              border=NA,
#              horiz=TRUE)
#tbl <- read.table("core_v0.7.pruneddata3.8.Q.reordered.txt", sep="\t")
#color <- read.table('core_v0.7.pruneddata3.8.Q.color.txt', sep="\t")
#barplot(t(as.matrix(tbl)), space = c(0.2), col=as.vector(color$V3), xlab="", ylab="", border=NA, horiz=TRUE, axes=FALSE)


par(mar=c(4,5,2,4))
i = 5
prefix = 'base_filtered_v0.7.pruneddata_1M_20kb_0.8.'
filename <- paste(prefix, i, ".Q.reordered.txt", sep="")
tbl <- read.table(filename, sep="\t")
colname  <- paste(prefix, i, ".Q.color.txt", sep="")
color <- read.table(colname, sep="\t")
barplot(t(as.matrix(tbl)), space = c(0.2), col=as.vector(color$V3), xlab="", ylab="", border=NA, horiz=TRUE, axes=FALSE)
#text(0.25, -200, labels = paste("K = ", i, sep=""), cex = 1.8, col="black", srt = 1, adj = 0, xpd = TRUE)

par(las=2)
mtext("TRJ", side = 2, font=1, at = 200, line= 3, cex=1, adj = 0, col="cornflowerblue")
mtext("TEJ", side = 2, font=1, at = 800, line= 3, cex=1, adj = 0, col="blue")
mtext("ARO", side = 2, font=1, at = 1090, line= 3, cex=1, adj = 0, col="darkorchid")
mtext("AUS", side = 2, font=1, at = 1250, line= 3, cex=1, adj = 0, col="chocolate")
mtext("IND", side = 2, font=1, at = 2500, line= 3, cex=1, adj = 0, col="green")
par(las=1)
mtext("A", side = 3, line = -1, at = -0.1, cex=1)
##plot ADMIXTRUE##########################################################################################################


#par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
#barplot(y,horiz=TRUE,width=1,space=0, xlim=c(-10, 200),  ylim=c(0.5,length(tree$tip.label)),names="", axes=FALSE)
#axis(1, at= c(0, 100, 200), line=1)

#pong
par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
barplot(y3,horiz=TRUE,width=1,space=0, xlim=c(0, 20),  ylim=c(0.5,length(sorted_list)),names="", axes=FALSE)
axis(1, at= c(0, 10, 20), line=0)
mtext("Pong", side=1,font=3, at=10,line=2.5, cex=1, col="black")
mtext("B", side=3, line = -1, at = -1.5, cex=1)

#ping
par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
barplot(y2,horiz=TRUE,width=1,space=0, xlim=c(0, 10),  ylim=c(0.5,length(sorted_list)),names="", axes=FALSE)
axis(1, at= c(0, 5, 10), line=0)
mtext("Ping", side=1,font=3, at=5,line=2.5, cex=1, col="black")
mtext("C", side=3, line = -1, at = -1.5, cex=1)

#mping
par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
cut <- 100
#new scale 300 to 600 cresponding to 120 and 150. So the factor is 300/30=10
#new value should be (503-140)/7.5 + 120
for(i in 1:length(y1)){
   if (y1[i] > 110){
       y1[i] <- (y1[i] - cut)/10 + 100
   }
} 

xx <- barplot(y1,horiz=TRUE,width=1,space=0, xlim=c(0, 150),  ylim=c(0.5,length(sorted_list)),names="", axes=FALSE)
axis(1, at=c(0, 30, 60, 90, 120, 150), labels=c(0, 30, 60, 90, 300, 600), line=0)
mtext("mPing", side=1,font=3, at=70,line=2.5, cex=1, col="black")
#break y
#break y
b <- 100
axis.break(1, b, style="slash")
rect(b, 0.2, b+2, max(xx)+0.8, border=FALSE, col='white')
mtext("D", side=3, line = -1, at = -1.5, cex=1)

#ping
#par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
#barplot(y2,horiz=TRUE,width=1,space=0, xlim=c(0, 10),  ylim=c(0.5,length(sorted_list)),names="", axes=FALSE)
#axis(1, at= c(0, 5, 10), line=0)
#mtext("Ping", side=1,font=3, at=5,line=2.5, cex=1, col="black")
#pong
#par(mar=c(4,0.5, 2, 1)) #set left and right to be tight with other plot
#barplot(y3,horiz=TRUE,width=1,space=0, xlim=c(0, 20),  ylim=c(0.5,length(sorted_list)),names="", axes=FALSE)
#axis(1, at= c(0, 10, 20), line=0)
#mtext("Pong", side=1,font=3, at=10,line=2.5, cex=1, col="black")
dev.off()

