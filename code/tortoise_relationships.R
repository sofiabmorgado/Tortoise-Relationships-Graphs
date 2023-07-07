#This project consists of a graph analysis. 
#The dataset is available in Kaggle website and pertains relationships in a tortoise community

library(igraph)

#Import the data
pth <- "C:/Users/sofia/OneDrive/Ambiente de Trabalho/projetoAGG/reptilia-tortoise-network-pv.zip" #Choose the correct path for your data
download.file(pth, destfile = "reptilia-tortoise-network-pv.zip")
# see file names
unzip("reptilia-tortoise-network-pv.zip", list = TRUE)
# unzip
unz <- unzip("reptilia-tortoise-network-pv.zip", "reptilia-tortoise-network-pv.edges")

#Take a look into the data
readLines(unz, n=10)

# skip first line
dat <- read.table(unz, skip=0, sep=","); dat
#os dados foram colhidos durante 7 anos, em 2005, 2008, 2009, 2010, 2011, 2012, 2013
#data2 <- read.graph("reptilia-tortoise-network-pv.edges")
#data <- read_graph("C:/Users/sofia/OneDrive/Ambiente de Trabalho/projetoAGG/reptilia-tortoise-network-pv.edges",
  #format = "edgelist")

#Remove years column
edges <- list()

for (i in 1:104){
  s <- strsplit(dat[i,], " ", fixed=T)
  num <- as.numeric(unlist(s))
  l <- c(num[1], num[2])
  edges <- append(edges, l)
}

edges_final <- as.numeric(unlist(edges))
edges_final

set.seed(1)
g <- make_graph(edges = edges_final, n=30, directed = FALSE)
plot(g, vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)


###SOME STATISTICS

#VERTICES AND EDGES
V(g)
E(g)
#Show edge list
as_edgelist(g)
#Show for each vertex the list of other vertices connected to it
as_adj_list(g)
#Show for each vertex the list edges connected to it
as_adj_edge_list(g)



#AJACENCY MATRIX
adj_m <- as_adjacency_matrix(g)
adj_m
df <- as.data.frame(as.matrix(adj_m))

write.table(df, file = "adj_m.png", sep = " ", row.names = FALSE, col.names = FALSE)



#DEGREE
degree(g)
degree_distribution(g)

max(degree(g))
mean(degree(g))
min(degree(g))
hist(degree(g), xlab = "Nodes Degree", main = "Nodes degree histogram", col="lightblue")

plot(degree_distribution(g), type="l", xlab = "Nodes Degree", ylab = "Nodes Degree Distribution")
plot(degree_distribution(g,cumulative=TRUE), type="l", xlab = "Nodes Degree", ylab = "Comulative Nodes Degree Distribution")
fit_power_law(degree(g))
?fit_power_law

#DENSITY
#For a simple gragh
g_simple <- simplify(g)
E(g_simple)
simple_dens <- edge_density(g_simple)
simple_dens

#Mean for each year
mean_dens_year <- edge_density(g) / 7
mean_dens_year



#PATHS AND WALKS
short_path <- shortest_paths(g, from=10, to=32, output = "both")
short_path


E(g)$color <- "Gray"
E(g)$color[unlist(short_path$epath)] <- "Red"
set.seed(1)
plot(g, vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)


#DIAMETER AND ECCENTRICITY
#Diametro do grafo - max valor das distancias minimas entre todos os vertices
distances(g)
diameter(g)

#Excentricidade - ver o maximo na coluna correspondente de distances(g0)
eccentricity(g) #distancia minima para qualquer ponto tendo em conta qualquer ponto
eccentricity(g) 
?eccentricity



#COMPONENTS
is_connected(g)
components(g, mode = c("weak", "strong"))
fraction_main <- components(g)$csize[1] / sum(components(g)$csize) 
fraction_main
  
decomp <- decompose(g, mode = c("weak", "strong"), max.comps = NA, min.vertices = 0)
sub_g_1 <- subcomponent(g, v=1)
sub_g_1
sub_g_2 <- subcomponent(g, v=19)
sub_g_2
sub_g_3 <- subcomponent(g, v=11)
sub_g_3


V(g)[sub_g_1]$color <- "dodgerblue4"
V(g)[sub_g_2]$color <- "deepskyblue"
V(g)[sub_g_3]$color <- "lightcyan"
E(g)$color <- "gray"
set.seed(1)
plot(g, vertex.size=15,edge.arrow.size=0.01,vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)



#CONECTIVITY
vertex_connectivity(g, 1, 25)
vertex_connectivity(g, 1, 15)

#cut set
min_cut(g, 1, 15)
min_cut(g, 1, 25)




#Laplacian matrix



#Centrality
##Degree centrality
centr_degree(g)

c_scale <- colorRamp(c('lightcyan','deepskyblue','dodgerblue4'))
V(g)$color <- apply(c_scale(degree(g)/max(degree(g))), 1,
                     function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
set.seed(1)
plot(g)



##Centrality eigenvector
centr_eigen(g)
V(g)$color <- apply(c_scale(centr_eigen(g)$vector), 1,
                    function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
set.seed(1)
plot(g)


##Centrality closeness
centr_clo(g)

##Centrality betweeness
centr_betw(g)



#Groups and communities

##Cliques 
cliques(g)
max_cliques(g)
count_max_cliques(g)
largest_cliques(g)
clique_num(g)

#plot largest cliques
set.seed(1)
plot(g,mark.groups =largest_cliques(g),
     vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)

##Cores
?coreness
coreness(g)
set.seed(1)
plot(g,  mark.groups = list(which(coreness(g
)>=7)),
vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
vertex.label.color = "black",
vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)




##Clustering coefficient
transitivity(g, type = "average")
?transitivity
triangulos <- count_triangles(g)
triangulos
?count_triangles
#C = triangulos /caminhos com compimento 2

##Clustering coefficient of a vertex
k_4 = degree(g)[4]
Clust_coef_vert_4 <- 1 * k *(k-1)/2 #, where k is the degree of the node
Clust_coef_vert_4

##Redundancy of two nodes
redundancy <- function(g){
  return(transitivity(g, type="local") * (degree(g) - 1))
}
redundancy(g)



##Similarity

##Jaccard
similarity(g, method = "jaccard")
similarity(g, method = "jaccard")


#Cosine
#Cosine similarity
Cosine_similarity <- function(g){
  l<-length(V(g))
  M<-matrix(rep(0,l*l),ncol=l)
  for(i in 1:l){
    for(j in 1:l){
      M[i,j]<-sum(g[i,]*g[,j])/sqrt(degree(g,i)*degree(g,j))
    }}
  return(M)
}
Cosine_similarity(g)


#Pearson similarity
Pearson_similarity <- function(g){
  l<- length(V(g))
  M<- matrix(rep(0,l*l),ncol=l)
  for(i in 1:l){
    for(j in 1:l){
      M[i,j]<-cor(g[i,],g[j,])
    }}
  return(M)
}
Pearson_similarity(g)

#Hamming distance
Hamming_distance<-function(g){
  l<-length(V(g))
  M<-matrix(rep(0,l*l),ncol=l)
  for(i in 1:l){ for(j in 1:l){
    M[i,j]<- sum((g[i,]*g[j,])^2)
  }}
  return(M)
}
Hamming_distance(g)


#Clusters
#Cluster  optimal
?cluster_optimal
cluster1<-cluster_optimal(g)
cluster1
set.seed(1)
plot(cluster1, g, vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
#Modulariedade - arestas in / arestas 
modularity(g,membership(cluster1))


#CLuster eigen
cluster2<- cluster_leading_eigen(g)
set.seed(1)
plot(cluster2, g,  vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
modularity(g,membership(cluster2))


#Cluster louvain
clust_louv <- cluster_louvain(g)
communities(clust_louv)

# Plotting the degree centrality network with the community detection
set.seed(1)
plot(clust_louv, g, vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
modularity(g,membership(clust_louv))

#Cluster fast greedy: fast-greedy community finding works only on graphs without multiple edges, Invalid value

#CLuster infomap
cluster5<- cluster_infomap(g)
set.seed(1)
plot(cluster5, g,vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
modularity(g,membership(cluster5))

#Cluster: nor possible (ex g1a)


#Betwenness 
set.seed(1)
cluster7<-cluster_edge_betweenness(g)
plot(cluster7, g,vertex.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
modularity(g,membership(cluster7))

#Hierarchical cluster
d<- dist(1- similarity(g,method ="jaccard"))
cc <- hclust(d)
plot(cc, xlab = "", ylab = "")
cluster8 <- cutree(cc, k = 3)
cluster8a<- make_clusters(g,cluster8)
set.seed(1)
plot(cluster8a, g,.size=15, vertex.color = "light blue",edge.arrow.size=0.01,
     vertex.label.color = "black",
     vertex.label.cex = 0.6, edge.width=2, edge.curved=0.2)
modularity(g,membership(cluster8a))
