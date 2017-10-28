#subsetting data without weight
mydata<-read.table("C:/Users/compaq/Dropbox/Lubies-FaoNepal/Data/VCdata/R_data_VC.csv", header=TRUE,   sep=",", row.names="MYN")
myvars <- c("DC_Name2","PE_D_IN_F_HA", "PE_MF_IN_F_HA")


newdata <- mydata[myvars]
newdata2<-edit(newdata) #conver third colomn to numeric#

newdata3 <- subset(newdata2, PE_MF_IN_F_HA >= 1) #observations with frequencies are selected# how to remove duplicates rows???#

#removing frequnecies from the table
write.table(newdata3,"C:/Users/compaq/Dropbox/Lubies-FaoNepal/Data/VCdata/newdata3.csv", sep=",")

newdata4<-read.table("C:/Users/compaq/Dropbox/Lubies-FaoNepal/Data/VCdata/newdata3.csv", header=TRUE,   sep=",", row.names="row.names")
#select only two colmns wthout frequency
myvars2<-c("DC_Name2","PE_D_IN_F_HA")
newdata5 <- newdata4[myvars2]

el=as.matrix(newdata5)




#netrok map##http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
#Import data to Igrapg and statnet, edgelist, duplicates in rows possible.
library(igraph)

el=as.matrix(newdata5)

g=graph.edgelist(el,directed=FALSE) # turns the edgelist into a 'graph object'

plot.igraph(g)

plot.igraph(g,layout=layout.circle)#change layout
plot.igraph(g,layout=layout.drl)
plot.igraph(g,layout=layout.fruchterman.reingold)

plot.igraph(g,layout=layout.fruchterman.reingold)
plot.igraph(g,layout=layout.kamada.kawai)
plot.igraph(g,layout=layout.lgl)
plot.igraph(g,layout=layout.reingold.tilford)

plot.igraph(g,layout=layout.sphere, vertex.color="red", vertex.size=5)


#http://www.r-bloggers.com/experiments-with-igraph/

ecount(g)#a graph with 47 vertices (nodes) and 472 edges (connections)
vcount(g)
diameter(g)
farthest.nodes(g)# average minimum distance between a pair of vertices in the 3dist


####test#
lc <- largest.cliques(g)
V(g)$label <- V(g)$name
g.lc <- induced.subgraph(g, lc[[1]])
png(filename = "lc.png")
plot(g.lc, layout=layout.fruchterman.reingold, vertex.color="gray60", vertex.size = 0, edge.arrow.size = 0.5, edge.color = "gray80")

dev.off()
####
sgc <- spinglass.community(g)
V(g)$membership <- sgc$membership
# found 4 communities 0, 1, 2, 3
V(g) [ membership == 0 ]$color <- "cyan"
V(g) [ membership == 1 ]$color <- "green"
V(g) [ membership == 2 ]$color <- "blue"
V(g) [ membership == 3 ]$color <- "red"
V(g)$size <- 4
V(g) [ name == "Poultry and eggs fomr farms and hatcheries" ]$size <- 20
png(filename = "tls.png", height = 800, width = 800)
plot(g, layout=layout.fruchterman.reingold, vertex.color=V(g)$color, vertex.size = V(g)$size, vertex.label = NA, edge.arrow.size = 0.5)
dev.off()

####


plot.igraph(g,vertex.label=V(g)$name,vertex.size=2,,vertex.label.color="black", vertex.label.font=1,vertex.color="red",edge.color="black")
