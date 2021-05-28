
## ALMACENAR DATOS Y EXTRAER METADATOS ##
#############
## La primer muestra:
# Los vectores vacios para almacenar metadatos
muestra <- c()
procedencia <- c()
sexo <- c()
numero <- c()

# La función, definida por el usuario
muestra_1 <- function(){
  m1 <<- read.table(file.choose(), header = T, sep = "") 
  names(m1) <<- c("gene_name", names(m1)[2]) 
  
  muestra[1] <<- names(m1)[2] 
  
  .start <- readline(prompt = "Procederemos a obtener la informacion.") 
  from. <- readline(prompt = "La muestra procede de ileum (I), recto (R) o colon (C) ?") 
  if(from. == "I" | from. == "i"){ 
    procedencia[1] <<- "ileum" 
  }else if(from. == "R" | from. == "r"){ 
    procedencia[1] <<- "recto" 
  }else if(from. == "C" | from. == "c"){ 
    procedencia[1] <<- "colonT" 
  }else{
    print("ERROR") 
  }
  
  sex. <- readline(prompt = "La muestra procede de un especimen macho (M) o hembra (H) ?") 
  if(sex. == "M" | sex. == "m"){
    sexo[1] <<- "macho"
  }else if (sex. == "H" | sex. == "h"){ 
    sexo[1] <<- "hembra"
  }else{
    print("ERROR")
  }
  
  num. <- readline(prompt = "Cual es el numero de muestra ?") 
  .num <- as.numeric(num.) 
  numero[1] <<- .num 
}

# Ejecución
muestra_1()


## Llenado de la base de datos y metadatos:
# La función, definida por el usuario:
muestras <- function(i){
  m <<- read.table(file.choose(), header = T, sep = "")
  m1 <<- cbind(m1, m[ ,2])
  names(m1)[i+1] <<- names(m)[2]
  
  muestra[i] <<- names(m)[2]
  
  .start <- readline(prompt = "Procederemos a obtener la informacion.") 
  from. <- readline(prompt = "La muestra procede de ileum (I), recto (R) o colon (C) ?") 
  if(from. == "I" | from. == "i"){ 
    procedencia[i] <<- "ileum" 
  }else if(from. == "R" | from. == "r"){
    procedencia[i] <<- "recto" 
  }else if(from. == "C" | from. == "c"){ 
    procedencia[i] <<- "colonT" 
  }else{
    print("ERROR") 
  }
  
  sex. <- readline(prompt = "La muestra procede de un especimen macho (M) o hembra (H) ?")
  if(sex. == "M" | sex. == "m"){
    sexo[i] <<- "macho"
  }else if (sex. == "H" | sex. == "h"){
    sexo[i] <<- "hembra"
  }else{
    print("ERROR")
  }
  
  num. <- readline(prompt = "Cual es el numero de la muestra ?") 
  .num <- as.numeric(num.) 
  numero[i] <<- .num 
}

# Ejecución:
n <- 12 #numero de muestras totales
for(i in 2:n){
  muestras(i)
}


## Definicion de los objetos de tipo data.frame:
# La base de datos:
base_datos <- m1
View(base_datos)
head(base_datos)

# Los metadatos:
metadatos <- data.frame(muestra, procedencia, sexo, numero)
View(metadatos)
head(metadatos)
#############

#############

## EL ANALISIS DE EXPRESION DIFERENCIAL ##
#############
BiocManager::install("DESeq2") 
library(DESeq2)


# Transformar la base de datos a una matriz:
matrix_data <- as.matrix(base_datos[ , -1]) 
rownames(matrix_data) <- base_datos$gene_name 
View(matrix_data)

# Renombrar los renglones del data.frame de metadatos:
metadata <- metadatos[ , -1]
rownames(metadata) <- metadatos$muestra
View(metadata)

# Generar el dataset de tipo DESeq:
axd <- DESeqDataSetFromMatrix(countData = matrix_data, 
                              colData = metadata,
                              design = ~procedencia,
                              tidy = F)

# Realizar el analisis:
axdDE <- DESeq(axd, test="LRT", reduced=~1)
#############

#############

## LOS RESULTADOS ##
#############
#Ver las tablas de resultados de disntintas comparativas:
#colonT_vs_ileum
CTvsI <- results(axdDE, contrast=c("procedencia", "colonT", "ileum"))
summary(CTvsI)

#recto_vs_colonT
RvsCT <- results(axdDE, contrast=c("procedencia", "recto", "colonT"))
summary(RvsCT)

#ileum_vs_recto
IvsR <- results(axdDE, contrast=c("procedencia", "ileum", "recto"))
summary(IvsR)

# Ver los datos ordenados y modificar los archivos:
#colonT_vs_ileum
CTvsI <- CTvsI[order(CTvsI$padj),]
head(CTvsI, 10)[, c("log2FoldChange", "padj")]

#recto_vs_colonT
RvsCT <- RvsCT[order(RvsCT$padj),]
head(RvsCT, 10)[, c("log2FoldChange", "padj")]

#ileum_vs_recto
IvsR <- IvsR[order(IvsR$padj),]
head(IvsR, 10)[, c("log2FoldChange", "padj")]
#############

#############

## VISUALIZACION DE LOS RESULTADOS ##
#############
# Ver los resultados via Volcano plot:
par(mfrow=c(1,1)) #arejustar el display de los graficos a 1
# Volcano plot basico
with(CTvsI, plot(log2FoldChange, -log10(pvalue),
                 pch=20, main="Volcano plot: colon transverso vs ileon", xlim=c(-3,3)))

with(RvsCT, plot(log2FoldChange, -log10(pvalue),
                 pch=20, main="Volcano plot: recto vs colon transverso", xlim=c(-3,3)))

with(IvsR, plot(log2FoldChange, -log10(pvalue),
                pch=20, main="Volcano plot: ileon vs recto", xlim=c(-3,3)))


# Agregamos formato de colores de al valor de p y FoldChange 
with(subset(CTvsI , padj<.01 ),
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="blue"))

with(subset(RvsCT , padj<.01 ),
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="blue"))

with(subset(IvsR , padj<.01 ),
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="blue"))


with(subset(CTvsI , padj<.01 & abs(log2FoldChange)>2), # Si quieremos que subtraiga los genes no solo por su valor de p ajustado sino también por el valor de la magnitud de cambio mayor a 2
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="red"))

with(subset(RvsCT , padj<.01 & abs(log2FoldChange)>2), 
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="red"))

with(subset(IvsR , padj<.01 & abs(log2FoldChange)>2), 
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="red"))


# Ver los conteos por gen, entre tratamientos:
#funcion auxiliar dentro de un ciclo, definida por el usuario
plotting <- function(i){
  plotCounts(axdDE, gene = rownames(res_data)[i],
             intgroup = "procedencia")
}

#funcion definida por el usuario
PC <- function(DESeq_res){
  res_data <<- DESeq_res
  par(mfrow=c(2,5)) #ajustar el display
  for(i in 1:10){
    plotting(i)
  }
}

# Graficar de manera automatizada por cada objeto de los resultados
PC(CTvsI) #colonT_vs_ileum

PC(RvsCT) #recto_vs_colonT

PC(IvsR) #ileum_vs_recto
#############

#############

## REDES DE MICROBIOTA ##
#############
# Instalar y cargar los paquetes necesarios:
install.packages("igraph")
install.packages("Hmisc")
install.packages("Matrix")

library(igraph)
library(Hmisc)
library(Matrix)

## Cargar los archvos a la sesion desde el ordeandor:
datos_otus <- read.csv(file.choose(), header=T, row.names = 1) #archivo "outdata.csv"
metadatos_taxa <-read.csv(file.choose(),header=T, row.names = 1) #archivo "otu_taxonomy.csv"

## Filtrar los datos:
datos_otus <- datos_otus[ ,colSums(datos_otus) >= 10]

## Generar la matriz de correlaciones:
cor_data <- rcorr(as.matrix(datos_otus), type="spearman") #objeto de tipo list

# Almacenar y filtrar la informacion en relacion al p-value:
#Se almacenan los nombres de las entradas(columnas) que cumplen con lo que evalua el codigo.
p_valued <- forceSymmetric(cor_data$P) #la funcion define como NA a todas las "auto-correlaciones"

# Extraer los metadatos de las muestras filtradas por p-value:
#Almacenar toda la informacion de los OTUs almacenados en "p_valued".
ft_metadata <- metadatos_taxa[rownames(p_valued), ,
                              drop = FALSE]

## Ajustar los datos:
#Almacenar los OTUs (columnas) que cumplen la condicion asignada.
vcP <- p_valued<0.05 #indicar el valor corte para el subscripting por p-value

#Almacenar los OTUs que pasaron los filtros, y tomar la lista de valores
#del estadistico R para esas columnas.
r_vals <- cor_data$r

#Almacenar los OTUs que cumplan las condiciones del subscripting en el objeto "r_vals".
cRP <- r_vals*vcP

#Almacenar los OTUs que cumplan con la ultima condicion, una correlacion mas alta que el 75%.
cRP <- abs(cRP)>0.75 #indicar el valor corte para el subscripting por correlacion
finale <- cRP*r_vals #hacer el subscripting y almacenar los datos filtrados

## Generar la matriz de adyacencias:
# Guardar el objeto como matriz en "matx"
matx <- as.matrix(finale)

# Renombrar columnas conforme a la informacion taxonomica almacenada en "ft_metadata"
colnames(matx) <- as.vector(ft_metadata$Family)
rownames(matx) <- as.vector(ft_metadata$Family)

## Crear el objeto igraph:
NTWR <- graph_from_adjacency_matrix(matx, mode="undirected",
                            weighted = TRUE, diag = FALSE) #diagonal=FALSE para evitar errores, la diagonal es NA

# Eliminar los nodos que no esten conectados
non_connected <- V(NTWR)[degree(NTWR) == 0] #encontrar y almacenar el index los nodos
NTWR <- delete.vertices(NTWR, non_connected) #eliminarlos 


#############

############# fin