
D_meta<- read.csv("D:/IIM A - ABA/Network Analysis/Nodes.csv")
D <- read.csv("D:/IIM A - ABA/Network Analysis/Connections.csv")



#Manage dataset
B<-as.data.frame(table(D)) # Create an edge weight column named "Freq"
B1<-subset(B,Freq>0) # Delete all the edges having weight equal to 0


#Create an igraph object from the data frames
library(igraph)
Pcount <-graph_from_data_frame(B1, directed = FALSE, vertices = D_meta)
E(Pcount)$weight<-E(Pcount)$Freq # Assigning edge attribute to each edge
class(Pcount)



g <- graph.data.frame(B1, directed=TRUE)
# Remove unnecessary margins
par(mar = c(0, 0, 0, 0))
plot(g, layout = layout.fruchterman.reingold, vertex.size = 8,vertex.label.cex =.9,vertex.label.dist=.9, edge.curved=0.1,
     edge.arrow.size = 0.1, vertex.label = V(Pcount)$Name,vertex.label.family ='Calibri',edge.color = 'lightgreen')


#igraph summary
Pcount
gsize(Pcount)
gorder(Pcount)

#4. Attributes
V(Pcount)$Class

#Nodelist
V(Pcount)

#Edgelist
E(Pcount)

#Adjacency matrix
K<-Pcount[c(1:89),c(1:89)]
K

#1. Degree centrality
Pcount_deg<-degree(Pcount,mode=c("All"))
V(Pcount)$degree<-Pcount_deg
mean(V(Pcount)$degree)
which.max(Pcount_deg)


#2. Eigenvector centrality
Pcount_eig <- evcent(Pcount)$vector
V(Pcount)$Eigen<-Pcount_eig
V(Pcount)$Eigen
which.max(Pcount_eig)

#3. Betweenness centrality
Pcount_bw<-betweenness(Pcount, directed = FALSE)
V(Pcount)$betweenness<-Pcount_bw
V(Pcount)$betweenness
which.max(Pcount_bw)
mean(V(Pcount)$betweenness)

#4. closeness centrality
Pcount_cls<-closeness(Pcount)
V(Pcount)$closeness<-Pcount_cls
V(Pcount)$closeness
which.max(Pcount_cls)

plot(V(Pcount)$closeness,type = 'l',ylim = c(0,20))


#calculating correlations between centralities
deg<- c(V(Pcount)$degree)
cls<- c(V(Pcount)$closeness)
eig<- c(V(Pcount)$Eigen)
btw<- c(V(Pcount)$betweenness)

cor_mat_data<-data.frame(deg,cls,eig,btw)
cor_mat<- cor(cor_mat_data)




library(igraph)

#1. Plotting a network with the degree centrality

set.seed(1001)
library(RColorBrewer) # This is the color library
pal<-brewer.pal(length(unique(V(Pcount)$Class)), "Set3") # Vertex color assigned per each class number
plot(Pcount,edge.color = 'lightgreen',vertex.label.cex =.8, 
     vertex.color=pal[as.numeric(as.factor(vertex_attr(Pcount, "Class")))],
     vertex.size = Pcount_deg, edge.width=E(Pcount)$weight,
     vertex.label.family ='Calibri',
     layout = layout.fruchterman.reingold)





meanarrayd<-array(0,dim = 100)
meanarrayc<-array(0,dim = 100)
meanarrayb<-array(0,dim = 100)
meanarraye<-array(0,dim = 100)
for (i in 1:100){
        g_np <- erdos.renyi.game(89, 265, type = "gnm")
        meanarrayd[i]<-mean(centr_degree(g_np)$res)
        meanarrayc[i]<-mean(centr_clo(g_np)$res)
        meanarrayb[i]<-mean(centr_betw(g_np)$res)
        meanarraye[i]<-mean(centr_eigen(g_np)$vector)
}

f<-degree_distribution(g_np)
m<-table(f)
b<-f/sum(f)
plot(b,type='l',xlab='Degree',ylab='Probability of degree',col='blue',main = 'Degree Distribution of Random network')


mean(meanarrayd)
mean(meanarrayb)

g_random <- erdos.renyi.game(89, 265, type = "gnm")
g_random_deg<-degree(g_random,mode=c("All"))
plot(g, vertex.label= NA, edge.arrow.size=0.2,vertex.size = 1, xlab = "Random Network: G(N,L) model")










meanarraydb<-array(0,dim = 100)
meanarraycb<-array(0,dim = 100)
meanarraybb<-array(0,dim = 100)
meanarrayeb<-array(0,dim = 100)
for (i in 1:100){
        g_sfnm<-barabasi.game(89, power = 1, m = 1, out.dist = NULL, out.seq = NULL, 
                              out.pref = FALSE, zero.appeal = 1, 
                              directed = FALSE, algorithm ="psumtree", start.graph = NULL)
        meanarraydb[i]<-mean(centr_degree(g_sfnm)$res)
        meanarraycb[i]<-mean(centr_clo(g_np)$res)
        meanarraybb[i]<-mean(centr_betw(g_np)$res)
        meanarrayeb[i]<-mean(centr_eigen(g_np)$vector)
}

f<-degree.distribution(g_sfnm)
fb<-table(f)
bb<-fb/sum(fb)
plot(bb,xlab = 'Degree',ylab = 'Probability of degree',main='Degree Distribution of Barabasi Albert Model',col='blue',type='l')

g_sfnm<-barabasi.game(89, power = 1.5, m = 1, out.dist = NULL, out.seq = NULL, 
                      out.pref = FALSE, zero.appeal = 1, 
                      directed = FALSE, algorithm ="psumtree", start.graph = NULL)
plot(g_sfnm, vertex.label= NA, edge.arrow.size=0.2,vertex.size = 1, xlab = "BA model")
