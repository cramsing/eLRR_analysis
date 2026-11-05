# --- eLRR tree v1--- #### 

## ---  tree  --- #### 

library(phangorn)
library(ggtree)
library(dplyr)
library(treeio)

tree <- read.iqtree("alignments.treefile")

tree_data <- fortify(tree)
# Check  tip labels first
tree_data$label

# Make sure metadata matches
#meta2 <- read.delim("all_grp_hitsv2.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library(readxl)
meta <- read_xlsx("../Antagonism_review.xlsx", sheet = "tree_metadata")
meta$label <- meta$tree_name

meta <- meta %>% select(label, everything()) #reoders columns so label is first
head(meta)

# Join tree and metadata 
tree_joined <- tree_data %>%
  left_join(meta, by = "label")



tree <- ggtree(tree_joined) +
  geom_tiplab(data = subset(tree_joined, isTip),
              aes(label = name), size = 5, hjust=-.25) +
  geom_tippoint(data = subset(tree_joined, isTip),
                aes(fill = kingdom), shape = 21, size = 4, alpha = 0.8, stroke = 0.8) +
  geom_nodelab(data = subset(tree_joined, !isTip),
               aes(label = UFboot), size =4, nudge_x = 0.05)

#shape 21 = circle with stroke outside and color inside - https://sape.inf.usi.ch/quick-reference/ggplot2/shape
#geom_treescale(fontsize = 10, linesize = 0.5, family = "Arial") + theme_tree()

# - domains ----


# Load pacakges
library(gggenomes)
library(ggplot2)
library(drawProteins)
library(Biostrings)

seqs <- readAAStringSet("../eLRRs.fasta")
# strip anything after a space in the header, keeping just the ID
seq_ids <- sub(" .*", "", names(seqs))
w <- Biostrings::width(seqs)
chaindf <- data.frame(
  seq_id = seq_ids, 
  type        = "CHAIN",
  start       = 1,
  end         = w
)

tq <- read_gff3("gffs/tqxxx.gff")
qst <- read_gff3("gffs/qxxst.gff")
lrrs <- read_gff3("gffs/lrrpred.gff")
interpro <- read_gff3("gffs/interpro.gff")
sp <- interpro %>% filter(name ==  "SIGNAL_PEPTIDE")

head(tq)
head(qst)
head(lrrs)
head(sp)
#select only essential columns
tq <- tq %>% select(seq_id, start, end, strand, type, source, score, sequence)
qst <- qst %>% select(seq_id, start, end, strand, type, source, score, sequence)
lrrs <- lrrs %>% select(seq_id, start, end, strand, type, source, score, sequence)
sp <- sp %>% select(seq_id, start, end, strand, type, source, score)
sp$sequence <- NA
#combined <- rbind(tq, qst, lrrs, sp) #got error there were unequal columns even though this was not true
combined <- bind_rows(tq, qst, lrrs, sp, chaindf)


head(combined)

#Rename columns
combined <- combined %>% rename( 'begin' = 'start', 
                         'description' = 'type',
                         'entryName' = "seq_id")
combined <- combined %>%
  mutate(type = case_when(description == "CHAIN" ~ 'CHAIN',
                          TRUE ~ "DOMAIN" ))

head(combined)


unique(combined$description)
#[1] "TQ..."            "Q[A-Za-z]{2}[TS]" "LRR"              "signal_peptide"   "CHAIN" 
combined <- combined %>% mutate(description = recode(description,
                                                     "TQ..." = "TQxxx",
                                                     "Q[A-Za-z]{2}[TS]" = "QxxS/T",
                                                     "signal_peptide" = "Signal Peptide"))


combined <- combined %>% filter(description != "TQxxx")
                                
# Get tip order from tree
tip_order <- tree_joined %>%
  filter(isTip) %>%
  arrange(y) %>%  
  pull(name) %>%
  unique()                                      

#NEED TO MATCH NAMES!! not there yet!!

#tip_order
#chaindf$seq_id
combined <- combined %>% mutate(entryName = recode(entryName,
                                                     "NbLRRP1_" = "NbLRRP1"))
# Reorder your protein data to match
combined <- combined %>%
  mutate(order = match(entryName, tip_order)) %>%
  arrange(order)

combined$order <- as.numeric(combined$order) #makes sure this column is set to numeric 
combined$end <- as.numeric(combined$end) #makes sure this column is set to numeric 
combined$begin <- as.numeric(combined$begin) #makes sure this column is set to numeric 

# --- Plot --- ####        

library(drawProteins)
# Your domain plot
p <- draw_canvas(combined) %>%
  draw_chains(combined, 
              outline      = "black",  # outline colour of each chain
              fill         = "white") %>%
  draw_domains(combined, label_domains = FALSE) 
p + theme_minimal(base_size = 14)


p2 <- draw_canvas(combined)
p2 <- p2 + ggplot2::geom_segment(
  data = combined,
  mapping = ggplot2::aes(x = begin, xend = end, y = order, yend = order),
  colour = "darkgray",
  linewidth = 1.5,       # thickness of the line itself
  lineend = "round"      # "round", "butt", or "square" end caps
)

# still add chain labels manually since you skipped draw_chains()'s labeling
p2 <- p2 + ggplot2::annotate("text", x = -10,
                           y = combined[combined$type == "CHAIN", ]$order,
                           label = combined[combined$type == "CHAIN", ]$entryName,
                           hjust = 1, size = 4)
p2 <- draw_domains(p2, combined,
                  label_domains = FALSE)
p2 + theme_minimal(base_size = 14)
p2 <- p2 + theme_minimal(base_size = 14)  +        # white background, bigger base text
     theme(
       panel.grid.minor = element_blank(),
       panel.grid.major = element_blank(),
       axis.ticks       = element_blank(),
      axis.text.y      = element_blank(),
       panel.border     = element_blank()
     )


p2

library(patchwork)

combo <-  tree + p2
# 
# tree_plot_clean <- tree + theme(legend.position = "none")
# # 
# domain_plot_clean  <- domain_plot + theme(legend.position = "none")
# 
# combo_clean <-  tree_plot_clean + domain_plot_clean

ggsave(combo, file =  "elrrs.pdf", width=15, height=5)


