# Author: Wyatt J. McDonnell, PhD
# (c) 2022 10x Genomics, Inc.
# Script provided for reproduction of Jaffe et al. __ 2022
# --- Libraries
library(data.table)
library(tidyverse)
library(msa)
library(ggseqlogo)
library(ggrepel)
library(scales)

# --- Data
curves = fread('curves_1.csv')
curvestwo = fread('curves_2.csv')
curvesthree = fread('curves_3.csv')
curvesfour = fread('curves_4.csv')
seqs = fread('seqlogofig.csv')

# --- Figure 1 panel (a)
ggplot(curves %>%
         filter(group %in% c('vh_name_cdrh3_length'), donor != 'any'),
       aes(x=threshold,y=percent,color=dref,shape=donor)) +
  geom_point(shape='circle') +
  geom_path() +
  scale_y_continuous(limits=c(1,100)) +
  xlab('CDRH3 amino acid percent identity') +
  ylab('Probability (%) of shared light chain V gene') +
  labs(color='Donor pair',shape = 'Distance from\nreference') +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c('#E84B50','#00A1DF'), name = '',
                     labels = c('Memory','Naive')) +
  theme(legend.position = c(0.25,0.75),
        legend.direction = 'vertical',
        legend.text = element_text(size=30))

# --- Figure 2 panel (a)
ggplot(curves %>%
         filter(group %in% 'cdrh3_length_diff_vh', donor != 'any'),
       aes(x=threshold,y=percent,color=dref,shape=donor)) +
  geom_point(shape='circle') +
  geom_path() +
  scale_y_continuous(limits=c(1,100)) +
  xlab('CDRH3 amino acid percent identity') +
  ylab('Probability (%) of shared light chain V gene') +
  labs(color='Donor pair',shape = 'Distance from\nreference') +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c('#E84B50','#00A1DF'), name = '',
                     labels = c('Memory','Naive')) +
  theme(legend.position = c(0.25,0.75),
        legend.direction = 'vertical',
        legend.text = element_text(size=30))

# --- Figure 3 panel (a) for weighted vs. simple matrices
ggplot(curvesthree,
       aes(coherence, pairs,
           color = method, label=identity)) +
  geom_point(alpha=0.5) +
  geom_line() +
  scale_y_log10(labels=comma, limits=c(1000,1000000)) +
  geom_label_repel(size = 6,
                   fontface = 'bold') +
  scale_x_continuous(limits = c(50,100)) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(legend.title = element_text(face='bold')) +
  theme(legend.background =
          element_rect(size=0.5,
                       linetype="solid",
                       colour ="black")) +
  xlab('Light chain coherence (%)') +
  ylab('Cell pairs') +
  labs(color='Method') +
  scale_color_manual(values = c('#E84B50','#00A1DF'))

# --- Figure 3 panel (b) for weighted vs. simple matrices
ggplot(curvesfour,
       aes(coherence, cells,
           color = method, label=identity)) +
  geom_point(alpha=0.5) +
  geom_path() +
  scale_y_log10(labels=comma, limits=c(1000,50000)) +
  geom_label_repel(max.overlaps = Inf, segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20,
                   size = 6,
                   fontface = 'bold') +
  scale_x_continuous(limits = c(50,100)) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(legend.title = element_text(face='bold')) +
  theme(legend.background =
          element_rect(size=0.5,
                       linetype="solid",
                       colour ="black")) +
  xlab('Light chain coherence (%)') +
  ylab('Cells in group') +
  labs(color='Method') +
  scale_color_manual(values = c('#E84B50','#00A1DF'))

# --- Figure 3 panel (c)
# color schemes
rasmol = make_col_scheme(
  chars=c('D','E','C','M','K','R','S','T','F','Y','N','Q','G','L','V','I','A',
          'W','H','P'),
  groups=c('Asp/Glu','Asp/Glu','Cys/Met','Cys/Met','Lys/Arg','Lys/Arg','Ser/Thr',
           'Ser/Thr','Phe/Tyr','Phe/Tyr','Asn/Gln','Asn/Gln','Gly','Leu/Val/Ile',
           'Leu/Val/Ile','Leu/Val/Ile','Ala','Trp','His','Pro'),
  cols=c('#E60A0A','#E60A0A','#E6E600','#E6E600','#145AFF','#145AFF','#FA9600',
         '#FA9600','#3232AA','#3232AA','#00DCDC','#00DCDC','#EBEBEB','#0F820F',
         '#0F820F','#0F820F','#C8C8C8','#B45AB4','#8282D2','#DC9682')
)
property = make_col_scheme(
  chars=c('D','E','C','M','K','R','S','T','F','Y','N','Q','G','L','V','I','A',
          'W','H','P'),
  groups=c('Acidic','Acidic','Sulfurous','Sulfurous','Basic','Basic','Hydroxylic',
           'Hydroxylic','Aromatic','Aromatic','Amidic','Amidic','Aliphatic',
           'Aliphatic','Aliphatic','Aliphatic','Aliphatic','Aromatic','Basic',
           'Aliphatic'),
  cols=c('#ec671a','#ec671a','#eeaa03','#eeaa03','#739dd3',
         '#739dd3','#ec95bc','#ec95bc','#aac74b','#aac74b',
         '#263481','#263481','#cd1719','#cd1719','#cd1719',
         '#cd1719','#cd1719','#aac74b','#739dd3','#cd1719')
)
# plot Figure 3 panel (c)
ggplot() +
  geom_logo(
    msa(seqs %>% pull(cdrh3),
        type='protein',
        method = 'Muscle')@unmasked[] %>%
      as.character(),
    col_scheme = property,
    method = 'bits') +
  theme_logo() +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ylab('')

ggplot() +
  geom_logo(
    msa(seqs %>% pull(cdrl3),
        type='protein',
        method = 'Muscle')@unmasked[] %>%
      as.character(),
    col_scheme = property,
    method = 'bits') +
  theme_logo() +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  ylab('')

# --- Extended Figure 5
ggplot(curvestwo, aes(ndn,`percent of data`, color=category)) +
  geom_point() +
  geom_line() +
  xlab('# of inserted bases in CDRH3') +
  ylab('Percent of sequences within sequence type') +
  labs(color='Antibody\nsequence\ntype') +
  theme_classic() +
  scale_y_log10(breaks=c(100,75,50,25,10,1,0.1))
