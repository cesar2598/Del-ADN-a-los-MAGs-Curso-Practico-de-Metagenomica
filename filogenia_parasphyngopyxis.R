#### Parasphingopyxis ####

# 📚 Cargar paquetes necesarios
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ape)

# 📂 Establecer directorio de trabajo
setwd("C:/Users/57300/Desktop/micromanglar/de_novo_wf_results/Pseudomonadota/Parasphingopyxis")

# 📥 Cargar árbol en formato Newick
tree <- read.tree("parasphingopyxis.tree")

# 🐛 Cargar archivo de nombres taxonómicos GTDB
gtdb <- read.csv("bac120_taxonomy.tsv", sep = "\t") %>%
  select(user_genome, classification) %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(across(everything(), ~str_remove(., "^[a-z]__"))) %>%
  mutate(Species = ifelse(is.na(species), user_genome, paste(genus, species))) %>%
  mutate(Species_parse = str_replace(Species, " ", "~"),
         Species_parse = paste0("italic(", Species_parse, ")"))

# 🧬 Unir información del árbol con taxonomía
tree_data <- tibble(label = tree$tip.label) %>%
  left_join(gtdb, by = c("label" = "user_genome"))

# 📊 Graficar árbol
p <- ggtree(tree, layout = "rectangular") %<+% tree_data +
  geom_tiplab(aes(label = Species_parse),
              parse = TRUE,
              size = 3.5,
              align = TRUE,
              linetype = "dotted",
              linesize = 0.3) +
  scale_x_continuous(name = "Distancia evolutiva basada en el ANI", expand = expansion(mult = c(0.05, 0.5))) +
  scale_y_continuous(name = "Especies bacterianas") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_blank()
  )

# Mostrar árbol
p

# Mostrar árbol
ggsave("arbol_parasphingopyxis.svg", plot = p, width = 12, height = 10, units = "in", dpi = 300)

