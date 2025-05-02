#Curso De la Metagenómica a los MAGs

####Instalación de paquetes ####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
install.packages("vegan") 
install.packages("RColorBrewer")
install.packages("ggplot2")

####Importación de archivos####

#Entrar a la ubicación donde está guardado el archivo .biom
setwd("C:ruta/al/archivo/biom")
#Cargar el archivo biom a la sesión de R
Curso_Metagenomica <- import_biom("Metatrans_completo.biom")

#Importamos el archivo de meta datos a la carpeta de trabajo
#el nombre de las muestras debe ser idéntico en ambos archivos para evitar errores
metadata <- read.csv("metadata.csv", row.names = 1) 

#Convertir el objeto a sample data
metadata <- sample_data(metadata)

#Agregar los metadatos al objeto phyloseq
Curso_Metagenomica <- merge_phyloseq(Curso_Metagenomica, metadata)

####Exploración de los datos ####

#Qué clase de objeto
class(Curso_Metagenomica)
#Cuántas anotaciones se realizaron
View(Curso_Metagenomica@tax_table@.Data)

#Eliminar identificadores y cambiar los nombres de las columnas
Curso_Metagenomica@tax_table@.Data <- substring(Curso_Metagenomica@tax_table@.Data, 4)
colnames(Curso_Metagenomica@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

View(Curso_Metagenomica@tax_table@.Data)

#Enlistar todos los phylum
unique(Curso_Metagenomica@tax_table@.Data[,"Phylum"])

#Cuántas anotaciones pertenecen al phylum Actinomycetota
sum(Curso_Metagenomica@tax_table@.Data[,"Phylum"] == "Actinomycetota")

#Clases del phylum Actinomycetota
unique(Curso_Metagenomica@tax_table@.Data[merged_metagenomes@tax_table@.Data[,"Phylum"] == "Actinomycetota", "Class"])

#Conteos por cada anotación
View(Curso_Metagenomica@otu_table@.Data)

####Filtro para más de 4 lecturas####

Curso_Metagenomica_min5 = prune_taxa(taxa_sums(Curso_Metagenomica) > 4, Curso_Metagenomica)

#visualizar las anotaciones y conteos nuevamente
View(Curso_Metagenomica_min5@tax_table@.Data)
View(Curso_Metagenomica_min5@otu_table@.Data)

Curso_Metagenomica
Curso_Metagenomica_min5

####Curvas de rarefaccion#### 

library(vegan)

otu.rare = otu_table(Curso_Metagenomica_min5)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)

otu.rarecurve = rarecurve(otu.rare, step = 10000, label = T)

####índices de diversidad alfa####

library(ggplot2)
#Básico
plot_richness(Curso_Metagenomica_min5, color = "Specie", shape = "Tissue")

#Personalizado
plot_richness(Curso_Metagenomica_min5, x = "Sample", color = "Specie", measures = c("Chao1", "observed", "shannon", "simpson")) + 
  geom_boxplot() + theme_bw() + ggtitle("Add Your Title Here") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Gráficos de abundancia
#Abundancia absoluta

#Gráfico de barras
psTopNOTUs = names(sort(taxa_sums(Curso_Metagenomica_min5), TRUE)[1:100])
pstop.prune = prune_taxa(psTopNOTUs, Curso_Metagenomica_min5)

plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Phylum")

plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#Abundancia relativa

#Transformados para abundancia relativa
#inspeccionar si hay lecturas no identificadas a los diferentes niveles
summary(Curso_Metagenomica_min5@tax_table@.Data=="")

#Eliminar las lecturas no identificadas a nivel de género

Curso_Metagenomica_min5_AR = subset_taxa(Curso_Metagenomica_min5, Genus !="")
summary(Curso_Metagenomica_min5_AR@tax_table@.Data=="")


head(Curso_Metagenomica_min5_AR@otu_table@.Data)
percentages_Gp = transform_sample_counts(Curso_Metagenomica_min5_AR, function(x) x*100/sum(x))
head(percentages_Gp@otu_table@.Data)

glom = tax_glom(percentages_Gp, taxrank = "Genus")
View(glom@tax_table@.Data)

percentages_Gp = psmelt(glom)
str(percentages_Gp)
percentages_Gp = percentages_Gp[,-4]
percentages_Gp =  as.data.frame(percentages_Gp)

#Reino
A_Relativa = ggplot(data=percentages_Gp, aes(x=Sample, y=Abundance, fill=Kingdom))+
  geom_bar(aes(), stat="identity", position="stack")+
  theme(aspect.ratio = 1, element_blank(), axis.title = element_text(size = 20), legend.text = element_text(size = 15))
A_Relativa

#Unir los menos abundantes
library("RColorBrewer")
percentages_Gp2 = percentages_Gp
#Phylum

#convertir a caracter la columna
percentages_Gp2$Phylum <- as.character(percentages_Gp2$Phylum) # Return the Phylum column to be of type character
#agrupar los phylum que representan menos del 5% de la abundancia relativa
percentages_Gp2$Phylum[percentages_Gp2$Abundance < 5] <- "Phyla < 5% abund."
unique(percentages_Gp2$Phylum)

#convertir la columna a factor
percentages_Gp2$Phylum <- as.factor(percentages_Gp2$Phylum)
#configurar los colores
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_Gp2$Phylum)))
#Graficar la abundancia realtiva 
relative_plot <- ggplot(data=percentages_Gp2, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")
relative_plot

#Género
percentages_Gp2$Genus <- as.character(percentages_Gp2$Genus) # Return the Phylum column to be of type character
percentages_Gp2$Genus[percentages_Gp2$Abundance < 5] <- "Genus < 5% abund."
unique(percentages_Gp2$Genus)

percentages_Gp2$Genus <- as.factor(percentages_Gp2$Genus)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_Gp2$Genus)))
relative_plot <- ggplot(data=percentages_Gp2, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")
relative_plot

####Diversidad Beta####

#Análisis de coordendas principales

# Normalizar los datos con base en los conteos
Curso_Metagenomica_norm <- transform_sample_counts(Curso_Metagenomica_min5, function(x) x / sum(x))

# Calcular PCoA con distancia de Bray-Curtis
ordu <- ordinate(Curso_Metagenomica_norm, method = "PCoA", distance = "bray")

# Visualizar (ajusta "TuVariable" al nombre de la variable de agrupación en tus metadatos)
plot_ordination(Curso_Metagenomica_rel, ordu, type = "Sample", 
                color = "Specie", shape = "Tissue") +
  geom_point(size = 5) +
  geom_text(aes(label = Sample), size = 5, vjust = -1) +  # aumenta el tamaño aquí
  theme_bw()

#PERMANOVA

#Extraer la tabla de conteos normalizados
Conteos_norm = Curso_Metagenomica_norm@otu_table@.Data
#Trasponer la tabla de conteos (Sino la función no corre)
Conteos_t= t(Conteos_norm)
#Convertir la tabla en un data frame
Conteos_t = as.data.frame(Conteos_t)

#Ejecutar la PERMANOVA
adonis2(Conteos_t ~ Specie, data = Meta_data, permutations = 9999, method = "bray")
adonis2(Conteos_t ~ Tissue, data = Meta_data, permutations = 9999, method = "bray")
adonis2(Conteos_t ~ Specie*Tissue, data = Meta_data, permutations = 9999, method = "bray")

