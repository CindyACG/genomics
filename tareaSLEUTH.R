########################################   CINDY A. CARRILLO GARCIA
########## GENOMICA FUNCIONAL ##########   LIC. EN MICROBIOLOGIA
########################################   SEMESTRE: 2021-I

######## TAREA SLEUTH Y KALLISTO ########

### DESCARGAR EL PAQUETE "SLEUTH" PARA R CON BIOCONDUCTOR ...
## Descargar el manager de Bioconductor para R: "BiocManager"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install() #instalarlo

## Instalar la libreria "devtools" de R, a partir de ella puede instalrse "seluth"
BiocManager::install("devtools") #instalacion via bioconductor on BiocManager para R

## Instalar la libreria "seluth" de R: "pachterlab/sleuth"
BiocManager::install("pachterlab/sleuth")
library("sleuth") #cargar a libreria al entorno

### ACCEDER A LOS DATOS A UTILIZAR Y ALMACENARLOS EN EL OBJETO "t2g" ...
## Generar una funcion que permita acceder a los  datos: tx2gene
# La funcion consiste de 4 pasos:
# En el primero, se accede al "mart" del cual vamos a obtener la informacion, nos
#permite conectarnos a BioMart de Ensembl (su base de datos datsets);indicamos que
#se utilize el de Ensembl con el argumento "biomart". Con el argumento "dataset"
#indicamos cual datset del BioMart de Ensembl queremos utilizar, en este caso es
#el de Homo Sapiens, en Ensembl podemos acceder a el con el nombre "hsapiens_gene_ensembl";
#contiene la informacion genetica del humano (por ello "gene").
#Estas instrucciones se almacenan en el objeto "mart". La subfuncion "useMart" de "biomaRt"
#permite acceder a un BioMart e indicar el "mart" al que queremos conectarnos, su
#funcion va mas alla e incluso nos permite acceder a un dataset dntro de este "mart"
#(base de datos).

# En el segundo, se obtiene informacion del dataset al que nos hemos conectado con
#las instrucciones del objeto "mart". La subfuncion "getBM" de "biomaRt" permite buscar
#y tomar la informacion de los atributos de los objetos/observaciones del dataset
#conforme a su estructura del BioMart al que nos conectamos. Permite buscar y almacenar
#estos atributos del dataset en que se ha indicado el argumento "mart", en este caso,
#el objeto "mart" del paso 1. Los atributos que quremos buscar y extraer del dataset
#se encuentran indicados cone el argumento "attributes"; al ser mas de uno, concatenamos
#los nombres de estos (por ello el "c()"). Los atributos que queremos obtener son
#aquellos que nos permitan saber a que gen vamos a analizar posteriormente, por lo
#tanto, los atributos son el ID del transcrito del gen, el ID y el nombre comun del
#gen a modo de poder conocer de cuales se trata para poder realizar investigacion
#posterior o corroborar informacion experimental.
#El output que obtendremos sera un dataset de esta informacion del dataset del
#mart, esta informacion la almacenamos en el objeto "t2g".

#En el tercero, con la funcion "dplyr" podemos modificar los elementos del dataset
#obtenido en el paso 2. Buscamos cambiar el nombre de los atributos (columnas) a
#uno mas reconocible tanto para nosotros como para funciones que seran utilizadas
#delante. Con la subfuncion "rename" de "dplyr" permite renombrar elementos de un
#dataframe o dataset, lo utilizamos para renombrar estas columnas. Se indica dentro
#de la funcion primero el nombre del objeto a modificar, en este caso, el dataset
#"t2g"; en segundo lugar se van nombrando los elementos a renombrar y tras el "="
#se indica el nuevo nombre del elemento.
#Se sobreescribe esta informacion en el miso objeto de modo que el dataset modificado
#conserve l mismo nombre "t2g".

#En el ultimo, se indica que el objeto "t2g" se imprima como output a modo de poder
#visualizarlo. Podemos hacer esto con la funcion "return()" de R. El objeto que
#se obtenga deberia ser el que tiene las columnas con nombres modificados.

tx2gene <- function(){
  #paso 1:
  mart <- biomaRt::useMart(biomart = "ensembl",
                           dataset = "hsapiens_gene_ensembl")
  
  #paso 2:
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                       "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  
  #paso 3:
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id",
                       ext_gene = "external_gene_name")
  
  #paso 4:
  return(t2g)
}

## Almacenar toda la informacion y sus modificaciones en un objeto: t2g
# Indicamos que el resultado de nuestra funcion se almacene en el objeto "t2g"
#a modo de manter un orden y evitar generar objetos de mas en el entorno.
t2g <- tx2gene()

### ACCEDER A LOS ARCHIVOS DE ABUNDANCIAS ...
## Indicar la ruta en el ordenador donde tenemos los archivos de kallisto: base_dir
# Almacenamos tal cual la direccion en el objeto "base_dir" que sera de tipo caracter.
base_dir <- "C:\Users\acer\Downloads\genomica\Archivo"

## Indicar los archivos especificos a utilizar: samples
# En el objeto "samples" almacenamos los nombres de nuestros archivos, sera igualmente
#un objeto de tipo caracter. Utilizamos la funcion "paste0()" a modo de ahorrarnos
#tiempo para escribir los nombres; lo utilizamos dado que los archivos se encuentran
#en carpetas que tienen nombres similares, la primer parte de estos nombres es identca
#por lo que indicamos este al inicio de "paste0()" ("sample"); utilizamos la concatenacion
#del resto del nombre, la arte variable de estos. La funcion "paste0()" a diferencia de
#"paste()" no deja ningun caracter o esacio entre los strings de caracter a unir, en este
#caso particular es util dado que el texto "1", "2" etc se encuentra en seguida del texto
#"sample" en el nombre de nuestras carpetas. Se entiende a "paste0()" como un analogo
#de "paste()" con la indicacion "" en el argumento "sep" dado que no tiene caracter de
#sepparacion.
samples <- paste0("sample", c("1", "2", "3",
                              "10", "11", "12"))

## Obtener y almacenar los directorios de nuestros archivos en un objeto: kal_dirs
# Con la funcion "sapply()" podemos obtener un vector de las direcciones de
#los archivos de kallisto (objetos de "sample") dado el vector de caracter "samples"
#del mismo tamaño y con el orden de "samples". La funcion tiene como input el vector
#y da como output un vector ordenado conforme al input y de su mismi largo; puede
#utilizarse tambien con listas y dataframes como input e inclusive puede generar matrices
#como output.
#El primer argumento de esta funcion es el objeto base, el inpt; en este caso es nuestro
#vector de nombres de los archivos de kallisto (objeto "samples"); el segundo argumento
#permite acceder a los archivos, se indica en "file.path()" la ruta de los mismos
#de modo que puedan ubicarse (objeto "base_dir"). La infoormacion adquirida se guarda
#en el objeto kal_dirs para poder utilizarse posteriormente.
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

### GENERAR UN DATAFRAME PARA EL ANALISIS DE EXPRESION DIFERENCIAL ...
## Hacer un dataframe a partir de la informacion de los archivos de kallisto: s2c
# Guardamos en el objeto "s2c" el dataframe que contenga la informacion sobre
#nuestras muestras analizadas con kallisto. Para generar este objeto utilizamos
#la funcion "data.frame()", para acceder a la informacion utilizamos el argumento
#"path"; el valor de este argumento es el vector de rutas "kal_dirs" generado con
#anterioridad. El argumento "sample" es de tipo tag=value, lo cual permite que
#indiquemos a "sample" como el nombre de la columna que contenga los valores mappeados
#con "path" para los archivos de nuestras muestras cuyos nombres se encuentran en 
#"samples". El argumento "muestras" es igual de tipo tag=value para el dataframe, en este
#caso indicamos los valores de la columna "muestras" manualmente. Esta parte nos permite
#indicar en el dataframe a que tratamiento (control=ctrl o tratamiento1=T1) pertenece la
#muestra a modo de utilizar esto para la comparativa de analisis de expresion.
#Se indica el argumento "StringsAsFactors" de modo que nuestro "string" de caracter
#"muestras" sea considerado como diferencia de grupos (como factores, los de control
#juntos y separados de T1.
s2c <- data.frame(path = kal_dirs, sample = samples,
                  muestras = c("ctrl","ctrl","ctrl", "T1", "T1", "T1"),
                  stringsAsFactors=FALSE)

### INDICAR LAS CONDICIONES PARA EL ANALISIS DE EXPRESION DIFERENCIAL ...
## Generar un archivo de tipo "sleuth", es decir, un grupo de archivos kallisto: so
# Podemos hacer dicho archivo con la funcion "sleuth_prep()", permite juntar un
#grupo de archivos de kallisto via metadatos; ademas de ello permite generar un modelo
#lineal para evaluar la informacion de expresion de cada una. Con esta funcion podemos
#juntar los archivos kallisto y analizarlos como un solo objeto.
#El primer argumento es el conjunto de archivos de kallisto a juntar, en este caso
#"s2c". Este objeto es "sample_to_covariates", este objeto requiere de la informacion
#que ya hemos almacenado en "sc2"; requiere de el "path" del objeto kallisto (su ruta) y
#de el "sample" que es la variable a la que son mappeadas las covariables "muestra".
#El segundo argumento es el modelo experimental a seguir, en este caso de evalua
#el efecto del factor "muestras" que separa a nuestros archivos en grupos; buscamos ver
#el efecto del grupo como diferencia entre muestras "~muestras": ¿lo observado depende
#del grupo?
#El tercer argumento es el dataframe que contiene el ID del gen target que se esta
#evaluando. En este dataframe dicho ID debe encontrarse en una columna llamada "target_id";
#esta columna se encuentra en el objeto "t2g" generado anteriormente y fue renombra a 
#"target_id" tambien anteriormente. Utilizamos como valor de este argumento a dicho objeto.
#El ultimo argumento, "extra_bootstrap_summary", permite indicar en el mismo dataframe
#la informacion referente al bootstrap generado para el calculo de las abundancias, en
#kallisto escritas como "estimated counts"; genera el resumen estadistico del remuestreo
#dado que se ha indicado como TRUE.
#La informacion, el dataframe, es generado y almacenado en el objeto "so".
so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,
                  extra_bootstrap_summary = TRUE)

## Ajustar el modelo lineal con un modelo de error de medicion
# Para ello utilizamos la funcion "seluth_fit()", se genera un modelo que permite
#estimar la variacion del bootstrap y ajustar a este modelo para evaluar el error
#de medicion asociado a este remuestreo. Ajustamos el modelo a esta variacion y es
#por ello que lo sobreescribimos. La funcion permite tomar como input un objeto de
#tipo sleuth por lo que especfica de lla libreria cargada con aterioridad.
so <- sleuth_fit(so)

## Hacer el test de Wald (chi^2) especificando el coeficiente beta
# El analisis se genera para cada transcrito evaluado; la funcion "seluth_wt()"
#permite determinar el valor de beta, la cndicion con base a la cual se realizara
#el test; en este caso de evalua el tratamiento 1. Permite determinar con que "gruping"
#se anlizaran los transcritos, para ello se indican los grupos que no son del control,
#en este caso solo se cuenta con uno, el T1; de contar con mas de uno -fuera del cntrol-
#deberiamos reaizar el analisis para cada grupo.
#El primer argumento es nuestro dataframe con la informacion necesaria para el analisis,
#esta almacenado en el objeto "so" por lo que este es el valor del argumento. El segundo
#argumento es el grouping para determinar el coeficiente beta, en este caso el tratamiento
#1 de nuestro factor de muestras "muestrasT1".
so <- sleuth_wt(so, which_beta="muestrasT1")                   

### CONOCER LOS RESULTADOS DEL ANALISIS DE EXPRESION DIFERENCIAL ...
## Visualizar los resultados en Rstudio con Shiny: FUN "sleuth_live()"
# Podemos ver los resultados del analisis realizado con las condiciones especificadas
#en el paso anterior, "sleuth_live()" nos permite utilizar la app "Shiny" para ver
#en tiempo real y de manera interactiva nuestros resultados. Su argumento es el objeto
#sleuth evaluado con "sleuth_fit()" y "sleuth_wt()", pr lo que fue una ventaja sobre
#-escribir el objeto; indicamos a "so".
sleuth_live(so)

### EVALUAR LOS RESULTADOS Y CONOCER LOS GENES SIGNIFICATIVOS (UP/DOWN) ...
## Acceder a la tabla de resultados del analisis del Sleuth
# Indicamos primero la ruta en que quremos ubicarnos, para ello utilizamos la funcion
#"setwd()"; se indica la ruta del work directory en que queremos estar en el ordenador.
setwd("C:\Users\acer\Downloads\genomica\Archivo")

# Ya que nos encontramos en dicho work directory, buscamos nuestra tabla de resultados,
#la cargamos a RStudio y a almacenamos en el objeto "resultados". Para ello, utilizamos
#la funcion "read.table"; el archivo es de texto y separado por comas. Los argumentos de
#la funcion incluyen (1°) el nombre del archivo/tabla a leer y cargar; (2°) el tipo de
#separador/separacion del archivo de texto, en este caso una coma "," dado que es un
#archivo .csv; (3°) si se han de manter los nombres de as columnas o solo los valores
#en cada una, en este caso es TRUE por lo que se mantendran los nombres.
resultados <- read.table("test_table.csv", sep=",", header=TRUE)

## Separar/fragmentar el dataframe conforme a resultados de valores significativos
significativos <- which(resultados$qval<0.1)
significativos <- resultados[significativos,]

## Separar/fragmentar nuevamente el dataframe conforme a los genes supraexpresados
upregulated <- which(significativos$b>0)
upregulated <- significativos[upregulated,]

## Separar/fragmentar nuevamente el dataframe conforme a los genes subexpresados
downregulated <- which(significativos$b<0)
downregulated <- significativos[downregulated,]


## Exportar en una tabla los genes sub- y supraexpresados

write.table(upregulated,file="C:\Users\acer\Downloads\genomica\Archivo\Upregulated_ctrlvsT1.txt",
            sep="\t")
write.table(downregulated,file="C:\Users\acer\Downloads\genomica\Archivo\Downregulated_ctrlvsT1.txt",
            sep="\t")

## Acceder a el nombre comun (nombre externo) de los genes supra- y subexpresados
upregulated$ext_gene
downregulated$ext_gene
