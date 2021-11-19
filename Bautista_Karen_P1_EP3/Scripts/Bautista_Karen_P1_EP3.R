
#Cargar librerias
library(ggtree)
library(treeio)
library(ape)
library(msa)

#Subir la secuencia
insulina<-readAAStringSet("Secuencias/Insulinas.fasta")

#1. Alineamiento C/ algoritmos
i_clustalW<-msa(insulina,method = "ClustalW") #alineamiento con clustal w
i_clustalW
i_clustalM<-msa(insulina,method = "Muscle") #alineamiento con muscle
i_clustalM

#2. Inferir arboles

### Neighbor Joining
## Esto es para convertir las secuencias ya alineadas
library(seqinr)#para convertir se debe cargar libreria seqinr
msaConvert(i_clustalW, type = "seqinr::alignment") -> i_clustalW
msaConvert(i_clustalM, type = "seqinr::alignment") -> i_clustalM

##Esto es para ponerlo como distancias
d_i_muscle<-dist.alignment(i_clustalM, "identity")
d_i_clustalW<-dist.alignment(i_clustalW, "identity")

##Esto es para sacar las matices
m_d_i_clustalW<-as.matrix(d_i_clustalW)
m_d_i_muscle<-as.matrix(d_i_muscle)

##Arboles filogeneticos

###Arbol Clustal W
TreeClustalW<- nj (m_d_i_clustalW)
plot(TreeClustalW, main="Arbol Filogenetico de Secuencias de Insulina (ClustalW)")

ggtree(TreeClustalW) + geom_tiplab(as_ylab=TRUE, color='royalblue')

###Arbol Muscle
TreeMuscle<- nj (m_d_i_muscle)
plot(TreeMuscle, main="Arbol Filogenetico de Secuencias de Insulina (Muscle)")

ggtree(TreeMuscle) + geom_tiplab(as_ylab=TRUE, color='darkslategrey')


