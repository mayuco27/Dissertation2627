###Code used to prepare all figures and analyses related to:
#Gómez-Gras et al.The Octocoral Trait Database: a global database of trait information for octocoral species 

# Libraries
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(ggplot2)
library(ggtreeExtra)
library(ggtree)
library(dplyr)
library(RRphylo)
library(aplot)
library(rstudioapi)

#We set the directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### Figure 2
##### Species-complete tree with phylogenetic distribution of trait data coverage (i.e., number of trait observations per species)

# We upload the family-resolved Maximum likelihood tree of Octocorallia inferred from 1059 bp alignment of mitochondrial gene mtMutS from McFadden et al. 2022 (iqtree.tre file)
coral_tree <- read.tree("iqtree.tre")

# We root it to add Hexacorallia as an outgroup
tree1 <- ape::root(coral_tree, node = 311)

#To include all species as polytomies, we need first to prune the tree to have only one tip per family.

pr.species<-c("OCT065_Echinoptilu","ANT63","OCT067_Pennatula","ANT78",
              "OCT023_Protoptilu","ANT32_Virgulari","OCT064_Chunellidae","ANT74_Ver",
              "OCT063_Antho","OCT066_Funicu","ANT58_Balt","OCT069_Umbel",
              "ANT59_Kop","ANT178_Dichotell","ANT22_An","ANT288_Ps","ANM078_Hete",
              "ANT85_An","ANM074_An","OCT010_Ant","OCT014_Pa","OCT123_Sibog","OCT153_Hem",
              "ANT37_Pa","OCT039_Sph","OCT144_Tita","OCT145_Scl","ANT56_Diodog","ANT277_Paras",
              "N110_Dasyg","N045_Chr","N088_C","S091_Chrys","T010_Chry","T062_Chr","S005_Chrys","N147_Pseu",
              "OCT116_M","OCT120_Ir","OG012_R","OCT130_R","ANT47_Prim","OCT133_Ca","OCT118_Ca","OCT119_N",
              "OCT053_Or","OCT131_K","OCT052_Jas","OCT029_Kera","OCT079_Kerato","OCT086_Ke",
              "OCT090_Kerato","OCT055_Pe","ANT28_Hel","ANT50_Par","OCT041_Clav","ANT48_Incr",
              "ANT70_Nid","OCT061_C","ANT72_Me","ANT62_Al","ANT19_Al","OCT142_Ele","OCT088_A",
              "OCT150_G","ANT20_Alc","ANT21_Alc","ANT60_Ant","OCT078_stolo","ANT41_Ca","ANT109_Acant",
              "ANT67_Sc","ANT173_S","OCT035_Lo","ANT38_Sar","OCT036_Lo","OCT044_stolon","ANT80_Aq",
              "ANT43_Tenero","OCT057_Lito","OCT146_P","OCT060_Scl","OCT058_C","OCT152_Duv",
              "OCT059_Pse","ANT55_Eu","ANT275_Su","OCT034_Kl","OCT076_Hete","OCT077_S",
              "ANT115_Pr","ANT98_A","ANT171_C","OCT037_R","OCT135_M","OCT134_Swi","ANT188_E","ANT31_Pl",
              "OCT050_Pinn","ANT84_Mu","OCT001_Dac","ANT116_R","ANT54_Bebr","ANT81_T","ANT76_P","OCT018_P","OCT070_Di",
              "OCT019_Ac","OCT021_A","ANT52_Ac","OCT071_Echi","ANT42_C","OCT033_An","OCT121_Call","OCT049","OCT101_Ca",
              "OCT047_An","OCT048_Gor","ANT265_Por","ANT182_Oculi","ANT286_Pa","ANT163_A","ANT229_El","ANT230_Nano",
              "ANT231_Pal","ANT131_A","ANT159_Sto","ANT205_Ali")

ii<-sapply(pr.species,grep,tree1$tip.label)
ii

pr.species<-tree1$tip.label[ii]
pr.species

# We visualize the resulting baseline trees before and after pruning
plotTree(tree1,fsize=0.7,lwd=1,
         ftype="i")
add.arrow(tree1,tip=ii,arrl=0.5)

tree1.noPR<-drop.tip(tree1,pr.species)
plotTree(tree1.noPR,fsize=0.7,lwd=1,
         ftype="i")

# We are only interested in the topology of the tree up to the family level, so we set all species
# within each family as equally related by selecting branch.length = "none". And visualize it.
ggpruned <- ggtree(tree1.noPR, branch.length = "none") + geom_tiplab(align=FALSE,linesize=0) +
  hexpand(.01)+
  geom_text2(aes(subset=!isTip, label=node))

ggpruned

# We save it, since it will be useful to examine nodes later on.
ggsave("pruned.pdf", width = 10,  height = 20,  units = "in", limitsize = FALSE)

# Now we can introduce all accepted species within each family tip. To do so:

#1) We import a df containing the list of all octcoral species for which trait information was collected in the database. 
sp_list1<-read.csv("sp_id.csv")

#2) We remove the few species not accepted yet in WORMS. 
sp_list1<-sp_list1[sp_list1$master_species!="Leptogorgia dictynna",]
sp_list1<-sp_list1[sp_list1$master_species!="Thouarella amundseni",]
sp_list1<-sp_list1[sp_list1$master_species!="Thouarella dolichoespinosa",]
sp_list1<-sp_list1[sp_list1$master_species!="Thouarella pseudoislai",]

#3) Any further species from families with just one species such as Coelogorgiidae,Taiaroidae,Rosgorgiidae,Haimeidae,Discophytidae, Isidoidae,Pseudonephtheidae,
#Corymbophytidae, Sinularidae, Malacanthidae or Skamnariidae needs to be included at this step since one tip is already present in the tree for these families. Family Aulopsammiidae (Lithotelestidae) will be included 
#within Helioporidae, as proposed by McFadden et al. 2022.
#**The family Pleurogorgiadae is not present in the Maximum likelihood tree of Octocorallia inferred from concatenated alignment of 739 UCE and exon
#loci (60% taxon occupancy, edge-trimmed), but the Maximum likelihood tree of Octocorallia inferred from 1059 bp alignment of mtMutS with tree topology
#constrained to match the topology of the target-enriched tree, indicates that this family is located aside Ifalukelidae. Therefore, we keep two tips for ifalukelidae rather than one, 
#leaving one as ifalukelidae (ANT46_Trichogorgia_capensis) and one as Pleurogorgiidae (Plumigorgia sp.), so both families will be represented. 

merger1<-sp_list1 %>%
  mutate(reference = case_when(
    family_molecules %in% "Acanthoaxiidae" ~ "OCT040_Altumia_delicata",
    family_molecules %in% "Subergorgiidae" ~ "ANT65_Annella_mollis",
    family_molecules %in% "Gorgoniidae" ~ "Pacifigorgia_irene",
    family_molecules %in% "Paramuriceidae" ~ "OCT012_Muriceides_cf._hirta",
    family_molecules %in% "Eunicellidae" ~ "ANT117_Eunicella_tricoronata",
    family_molecules %in% "Plexauridae" ~ "OCT073_Muricea_atlantica",
    family_molecules %in% "Pterogorgiidae" ~ "OCT007_Stenogorgia_casta",
    family_molecules %in% "Plexaurellidae" ~ "ANT100_Plexaurella_nutans",
    family_molecules %in% "Primnoidae" ~ "OCT097_Plumarella_sp.",
    family_molecules %in% "Coralliidae" ~ "ANM009_Bathyalcyon_sp.",
    family_molecules %in% "Keratoisididae" ~ "OCT022_Acanella_arbuscula",
    family_molecules %in% "Briareidae" ~ "ANT51_Briareum_hamrum",
    family_molecules %in% "Melithaeidae" ~ "ANT68_Solenocaulon_sp.",
    family_molecules %in% "Spongiodermidae" ~ "ANT276_Homophyton_verrucosum",
    family_molecules %in% "Erythropodiidae" ~ "ANT27_Erythropodium_caribaeorum",
    family_molecules %in% "Keroeididae" ~ "ANT112_Keroeides_sp.",
    family_molecules %in% "Ellisellidae" ~ "ANT83_Ellisella_sp.",
    family_molecules %in% "Isididae" ~ "ANT114_Isis_hippuris",
    family_molecules %in% "Mopseidae" ~ "OCT054_Mopseidae_sp.",
    family_molecules %in% "Nephtheidae" ~ "ANT44_Dendronephthya_sinaiensis",
    family_molecules %in% "Alcyoniidae" ~ "ALC43_Alcyonium_coralloides",
    family_molecules %in% "Helioporidae" ~ "ANT79_Nanipora_cf._kamurai",
    family_molecules %in% "Sarcophytidae" ~ "D449_Sclerophytum_humilis_n._comb.",
    family_molecules %in% "Capnellidae" ~ "OCT151_Pseudodrifa_racemosa",
    family_molecules %in% "Chrysogorgiidae" ~ "B010A_Chrysogorgia_sp.",
    family_molecules %in% "Tubiporidae" ~ "ANT36_Tubipora_sp.",
    family_molecules %in% "Siphonogorgiidae" ~ "ANT53_Siphonogorgia_godeffroyi",
    family_molecules %in% "Euplexauridae" ~ "OCT072_Euplexaura_rhipidalis",
    family_molecules %in% "Ifalukellidae" ~ "ANT46_Trichogorgia_capensis",
    family_molecules %in% "Astrogorgiidae" ~ "ANT66_Astrogorgia_sp.",
    family_molecules %in% "Cerveridae" ~ "ANT108_Azoriella_bayeri",
    family_molecules %in% "Carijoidae" ~ "OCT042_Cryptophyton_goddardi",
    family_molecules %in% "Clavulariidae" ~ "ANT24_Clavularia_inflata",
    family_molecules %in% "Cladiellidae" ~ "ANT69_Cladiella_pachyclados",
    family_molecules %in% "Acrophytidae" ~ "ANT105_Acrophytum_claviger",
    family_molecules %in% "Cornulariidae" ~ "ANT26_Cornularia_pabloi",
    family_molecules %in% "Lemnaliidae" ~ "ANT102_Paralemnalia_thyrsoides",
    family_molecules %in% "Xeniidae" ~ "ANT104_Ovabunda_macrospiculata",
    family_molecules %in% "Acrossotidae" ~ "ANT73_Acrossota_amboinensis",
    family_molecules %in% "Nidaliidae" ~ "OCT062_Nidalia_cf._occidentalis",
    family_molecules %in% "Scleracidae" ~ "OCT149_Scleracis_guadalupensis",
    family_molecules %in% "Aquaumbridae" ~ "OCT148_Elbeenus_lauramartinae",
    family_molecules %in% "Leptophytidae" ~ "ANT39_Leptophyton_benayahui",
    family_molecules %in% "Nephthyigorgiidae" ~ "ANT97_Nephthyigorgia_sp.",
    family_molecules %in% "Paralcyoniidae" ~ "ANT107_Ceeceenus_quadrus",
    family_molecules %in% "Pleurogorgiidae" ~ "OCT051_Plumigorgia_sp.",
    family_molecules %in% "Chelidonisididae" ~ "ANT82_Chelidonisis_aurantiaca_mexicana",
    family_molecules %in% "Anthogorgiidae" ~ "ANT77_Muricella_sp.",
    family_molecules %in% "Dendrobrachiidae" ~ "ANT75_Dendrobrachia_sarmentosa",
    family_molecules %in% "Ideogorgiidae" ~ "ANT61_Ideogorgia_capensis",
    family_molecules %in% "Parasphaerascleridae" ~ "ANT30_Parasphaerasclera_valdiviae",
    family_molecules %in% "Parisididae" ~ "ANT113_Parisis_sp.",
    family_molecules %in% "Victorgorgiidae" ~ "OCT046_Trachythela_rudis",
    family_molecules %in% "Arulidae" ~ "ANT99_Arula_petunia",
    family_molecules %in% "Incrustatidae" ~ "OCT043_Inconstantia_procera",
    family_molecules %in% "Sarcodictyonidae" ~ "ANT174_Telestula_n._sp.",
    family_molecules %in% "Aulopsammiidae" ~ "ANT79_Nanipora_cf._kamurai",
    family_molecules %in% c("Gyrophyllidae","Umbellulidae", "Anthoptilidae", "Balticinidae", "Pennatulidae", "Veretillidae", "Virgulariidae", "Renillidae", "Kophobelemnidae", "Funiculinidae", "Chunellidae", "Echinoptilidae", "Protoptilidae", "Pseudumbellulidae", "Scleroptilidae", "Stachyptilidae") ~ "ANT111_Renilla_muelleri" 
  ))

# There are also some species in the tree whose family assignation is uncertain. Since in Figure 2 we will only be showing species with family-assignation resolved, we remove these species. We do so by ommiting NAs
merger1<-na.omit(merger1)

# We are only interested in reference, and master_species (bind) columns so we subset the df
merger1<-subset(merger1, select = c(master_species, reference))
colnames(merger1) <- c("bind","reference")

# We also need to add the column poly = logical specifying if 'bind' and 'reference' should form a polytomous clade.
merger1$poly = FALSE

#Finally, we need to remove from the df those species that are already in the backbone tree to avoid duplicate values in the merged tree
toremove<-c(
  "Ellisella elongata", "Bathyalcyon robustum", "Telestula batoni",
  "Parisis australis", "Chrysogorgia carolinensis", "Plumarella delicatissima",
  "Tethrisis suzannae", "Tubipora syringa",
  "Solenocaulon tubulosum", "Nephthyigorgia annectens", "Keroeides pallida",
  "Astrogorgia splendens", "Muricella vasseuri", "Plexaurella nutans",
  "Scleracis guadalupensis", "Eunicella tricoronata",
  "Muriceides hirta", "Euplexaura rhipidalis",
  "Isis hippuris", "Stenogorgia casta", "Muricea atlantica",
  "Paralemnalia thyrsoides", "Ovabunda macrospiculata",
  "Cladiella pachyclados", "Annella mollis", "Acrophytum claviger",
  "Leptophyton benayahui", "Trachythela rudis", "Dendronephthya sinaiensis",
  "Pseudodrifa racemosa", "Ceeceenus quadrus",
  "Elbeenus lauramartinae", "Altumia delicata", "Cryptophyton goddardi",
  "Sclerophytum humile",
  "Inconstantia procera", "Nidalia occidentalis", "Siphonogorgia godeffroyi",
  "Azoriella bayeri", "Arula petunia", "Clavularia inflata",
  "Homophyton verrucosum", "Briareum hamrum", "Ideogorgia capensis",
  "Erythropodium caribaeorum", "Parasphaerasclera valdiviae",
  "Chelidonisis aurantiaca", "Trichogorgia capensis", "Acanella arbuscula",
  "Nanipora kamurai", "Dendrobrachia sarmentosa",
  "Renilla muelleri", "Cornularia pabloi", "Alcyonium coralloides",
  "Pleurogorgia plana", "Acrossota amboinensis", "Pacifigorgia irene"
)

merger2 <- merger1[!(merger1$bind %in% toremove),]

# Now we can proceed to add the new tips 
merged<-tree.merger(backbone=tree1.noPR,data=merger2,plot=FALSE)

# We plot it and save it as a complementary file that will be useful later
ggpruned2 <- ggtree(merged, branch.length = "none") + geom_tiplab(align=FALSE,linesize=0) +
  hexpand(.01)+
  geom_text2(aes(subset=!isTip, label=node))
ggpruned2
ggsave("pruned_MOD.pdf", width = 10,  height = 200,  units = "in", limitsize = FALSE)

# We now need to replace the labels of the few tips that are still showing genus information (from McFadden et al. 2022) with
#those species names that are still missing within each family. We also need to replace names from McFadden et al. 2022 with the species names accepted in WORMS.

library(plyr)

oldlabels<-merged$tip.label
my_vec_new1 <- revalue(oldlabels,c(ANT83_Ellisella_sp.="Ellisella elongata",ANM009_Bathyalcyon_sp.="Bathyalcyon robustum",ANT174_Telestula_n._sp.="Telestula batoni",
                                   ANT113_Parisis_sp.="Parisis australis",B010A_Chrysogorgia_sp.="Chrysogorgia carolinensis",OCT097_Plumarella_sp.="Plumarella delicatissima",OCT054_Mopseidae_sp.="Tethrisis suzannae",
                                   ANT11_Cerianthus_membranaceus="Hexacorallia",ANT36_Tubipora_sp.="Tubipora syringa",ANT68_Solenocaulon_sp.="Solenocaulon tubulosum",ANT97_Nephthyigorgia_sp.="Nephthyigorgia annectens",
                                   ANT112_Keroeides_sp.="Keroeides pallida",ANT66_Astrogorgia_sp.="Astrogorgia splendens",ANT77_Muricella_sp.="Muricella vasseuri",
                                   ANT100_Plexaurella_nutans="Plexaurella nutans",OCT149_Scleracis_guadalupensis="Scleracis guadalupensis",ANT117_Eunicella_tricoronata="Eunicella tricoronata",ANT45_Taiaroa_tauhou="Taiaroa tauhou",
                                   OCT012_Muriceides_cf._hirta="Muriceides hirta",OCT072_Euplexaura_rhipidalis="Euplexaura rhipidalis",ANT106_Discophyton_rudyi="Discophyton rudyi",ANT114_Isis_hippuris="Isis hippuris",OCT007_Stenogorgia_casta="Stenogorgia casta",
                                   OCT073_Muricea_atlantica="Muricea atlantica",ANT102_Paralemnalia_thyrsoides="Paralemnalia thyrsoides",ANT104_Ovabunda_macrospiculata="Ovabunda macrospiculata",ANT25_Coelogorgia_palmosa="Coelogorgia palmosa",ANT69_Cladiella_pachyclados="Cladiella pachyclados",
                                   ANT65_Annella_mollis="Annella mollis",ANT105_Acrophytum_claviger="Acrophytum claviger",ANT39_Leptophyton_benayahui="Leptophyton benayahui",ANT40_Corymbophyton_bruuni="Corymbophyton bruuni",OCT140_Pseudonephthea_liltvedi_n._comb="Pseudonephthea liltvedi",
                                   OCT147_Skamnarium_complanatum="Skamnarium complanatum",OCT046_Trachythela_rudis="Trachythela rudis",ANT44_Dendronephthya_sinaiensis="Dendronephthya sinaiensis",OCT151_Pseudodrifa_racemosa="Pseudodrifa racemosa",ANT107_Ceeceenus_quadrus="Ceeceenus quadrus",
                                   ANT64_Sinularia_brassica="Sinularia brassica",OCT148_Elbeenus_lauramartinae="Elbeenus lauramartinae",OCT040_Altumia_delicata="Altumia delicata",OCT042_Cryptophyton_goddardi="Cryptophyton goddardi",D449_Sclerophytum_humilis_n._comb.="Sclerophytum humile",
                                   ANT101_Malacacanthus_capensis="Malacacanthus capensis",OCT143_Rosgorgia_inexspectata="Rosgorgia inexspectata",OCT043_Inconstantia_procera="Inconstantia procera",OCT062_Nidalia_cf._occidentalis="Nidalia occidentalis",ANT53_Siphonogorgia_godeffroyi="Siphonogorgia godeffroyi",
                                   ANT108_Azoriella_bayeri="Azoriella bayeri",ANT99_Arula_petunia="Arula petunia",ANT24_Clavularia_inflata="Clavularia inflata",ANT276_Homophyton_verrucosum="Homophyton verrucosum",ANT51_Briareum_hamrum="Briareum hamrum",
                                   ANT61_Ideogorgia_capensis="Ideogorgia capensis",ANT27_Erythropodium_caribaeorum="Erythropodium caribaeorum",ANT30_Parasphaerasclera_valdiviae="Parasphaerasclera valdiviae",ANT82_Chelidonisis_aurantiaca_mexicana="Chelidonisis aurantiaca",ANT46_Trichogorgia_capensis="Trichogorgia capensis",
                                   OCT022_Acanella_arbuscula="Acanella arbuscula",ANT170_Isidoides_armata="Isidoides armata",ANT79_Nanipora_cf._kamurai="Nanipora kamurai",ANT75_Dendrobrachia_sarmentosa="Dendrobrachia sarmentosa",ANT111_Renilla_muelleri="Renilla muelleri",ANT26_Cornularia_pabloi="Cornularia pabloi",
                                   ALC43_Alcyonium_coralloides="Alcyonium coralloides",OCT051_Plumigorgia_sp.="Pleurogorgia plana",ANT73_Acrossota_amboinensis="Acrossota amboinensis", Pacifigorgia_irene="Pacifigorgia irene"))



# We check in one of them that the new tags are there, and that the old tags are no longer there
"Ellisella elongata" %in% my_vec_new1
'ANT83_Ellisella_sp.' %in% my_vec_new1

# Now we substitute the merged%tip.labels old by the new vector
merged_fin<-merged
merged_fin$tip.label <- my_vec_new1

#Now we can import the database dataset
library("dplyr")

OTD<-read.csv("OctocoralTraits_v2_2.csv")

# Before we proceed with further analyses, we remove families with family assignation uncertain, as well as the species from recent studies that have not yet been accepted in WORMS

OTD1<-OTD[OTD$Family_molecules!="Malacalcyonacea_incertae_sedis",]
OTD2<-OTD1[OTD1$Family_molecules!="Scleralcyonacea_incertae_sedis",]
OTD3<-OTD2[OTD2$Family_molecules!="incertae_sedis",]
OTD4<-OTD3[OTD3$Family_molecules!="uncertain",]
OTD5<-OTD4[OTD4$Family_molecules!="Haimeidae",]
OTD6<-OTD5[OTD5$specie_name!="Leptogorgia dictynna",]
OTD7<-OTD6[OTD6$specie_name!="Thouarella amundseni",]
OTD8<-OTD7[OTD7$specie_name!="Thouarella dolichoespinosa",]
OTD9<-OTD8[OTD8$specie_name!="Thouarella pseudoislai",]

# We also integrate Aulopsammiidae within Helioporidae (as proposed in McFadden et al. 2022)

OTD9$Family_molecules[OTD9$Family_molecules == "Aulopsammiidae"] <- "Helioporidae"

# To prepare Figure 2, we calculate the Number of trait observations per species

count_OTD <- aggregate(data = OTD9,                # Applying aggregate
                       observation_id ~ specie_name+Family_molecules,
                       function(observation_id) length(unique(observation_id)))

count_OTD$observation_id=as.numeric(count_OTD$observation_id)

# We group sea pen families. We do so in a new column that we will call superfamily

# Step 1: Define the levels to be grouped in the case of the superfamily pennatuloidea
levels_to_group <- c("Umbellulidae", "Anthoptilidae", "Balticinidae", "Pennatulidae", "Veretillidae", "Virgulariidae", "Renillidae", "Kophobelemnidae", "Funiculinidae", "Chunellidae", "Echinoptilidae", "Protoptilidae", "Pseudumbellulidae", "Scleroptilidae", "Stachyptilidae", "Gyrophyllidae")

# Step 2: Add a new column with grouped levels
count_OTD1 <- count_OTD %>%
  mutate(superfamily = ifelse(Family_molecules %in% levels_to_group, "Pennatuloidea", Family_molecules))

# In Figure 2, we will only color families harboring >75 species. We create an object that shows the N of species per family.

count_fam <- aggregate(data = count_OTD1,                
                       observation_id ~ superfamily,
                       function(observation_id) length(observation_id))
count_fam <- count_fam[order(count_fam$observation_id, decreasing = TRUE), ]


# We create a vector of the maintained levels only.

maintained_levels_short <- c("Nephtheidae",
                       "Paramuriceidae",
                       "Primnoidae",
                       "Sarcophytidae",
                       "Pennatuloidea",
                       "Gorgoniidae",
                       "Xeniidae",
                       "Melithaeidae",
                       "Alcyoniidae",
                       "Ellisellidae",
                       "Coralliidae",
                       "Chrysogorgiidae",
                       "Plexauridae",
                       "Mopseidae",
                       "Cladiellidae")

# Create a new column based on maintained and grouped levels
count_OTD_grouped_short <- count_OTD1 %>%
  mutate(grouped_family = ifelse(superfamily %in% maintained_levels_short, superfamily, "0_other"))

count_OTD_grouped_short$observation_id=as.numeric(count_OTD_grouped_short$observation_id)

# Some species have too many observations, so they do not fit in the plot. We find those with over 200 observations. Then, we collapse them to an arbitrary number (210), that will be the maximum showed in the plot. 
# These species are:
# Keratoisis grayi (position 1131 - N of trait observations = 4709)
# Muricea californica (position 2362- N of trait observations=437)
# Corallium rubrum (position 558- N of trait observations=283)

#Keratoisis grayi
count_OTD_grouped_short[1131, "observation_id"] <- 210

#Muricea californica
count_OTD_grouped_short[2362, "observation_id"] <- 210

##Corallium rubrum 
count_OTD_grouped_short[558, "observation_id"] <- 210

#Plot Figure 2 (for 15 families represented, >75 species)
#Node Malaca vs Sclera = 3461

tree31_short <- groupClade(merged_fin, c(3461))
p12_short <- ggtree(tree31_short, aes(color=group),branch.length = "none",layout="circular", size=0.1) +
  scale_color_manual(values=c("#979198","black"))+
  #geom_tiplab(align=FALSE,linesize=0) +
  hexpand(.01)

p12_short+ geom_fruit(data=count_OTD_grouped_short,
                axis.params=list(axis="x",text.size=2,limits=c(0,150),title="N of trait observations",title.height = 0.05, title.size=2),
                grid.params=list(vline=FALSE,size=0.15,color="lightgrey"),
                geom=geom_bar,
                mapping=aes(y=specie_name, x=observation_id, fill=grouped_family),
                pwidth=1.2,
                stat="identity",
                orientation="y") +
  scale_fill_manual(values=c("lightgrey","#EDF8B1", "#BD0026", "#41B6C4", "#c36721", "#F3d643",
                             "#253494", "#FDEBB6", "#E31A1C", "#7FCDBB", "#081D58",
                             "#FEB24C", "#225EA8", "#800026", "#C7E9B4" ,"#1D91C0"),
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
  theme(legend.justification = "right",
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.02, "cm"))+
  guides(color = "none")

ggsave("Figure_2.pdf", width = 8,  height =8,  units = "in", limitsize = FALSE)

# Figure 3
library(plyr)

# In this figure, each tip of the tree will represent a family in the plot. We create a vector of family tip names from the species-specific one that we have in the pruned tree.
oldlabelspruned<-tree1.noPR$tip.label
my_vec_new1prun <- revalue(oldlabelspruned,c(ANT83_Ellisella_sp.="Ellisellidae",ANM009_Bathyalcyon_sp.="Coralliidae",ANT174_Telestula_n._sp.="Sarcodictyonidae",
                                             ANT113_Parisis_sp.="Parisididae",B010A_Chrysogorgia_sp.="Chrysogorgiidae",OCT097_Plumarella_sp.="Primnoidae",OCT054_Mopseidae_sp.="Mopseidae",
                                             ANT11_Cerianthus_membranaceus="Hexacorallia",ANT36_Tubipora_sp.="Tubiporidae",ANT68_Solenocaulon_sp.="Melithaeidae",ANT97_Nephthyigorgia_sp.="Nephthyigorgiidae",
                                             ANT112_Keroeides_sp.="Keroeididae",ANT66_Astrogorgia_sp.="Astrogorgiidae",ANT77_Muricella_sp.="Anthogorgiidae",
                                             ANT100_Plexaurella_nutans="Plexaurellidae",OCT149_Scleracis_guadalupensis="Scleracidae",ANT117_Eunicella_tricoronata="Eunicellidae",ANT45_Taiaroa_tauhou="Taiaroidae",
                                             OCT012_Muriceides_cf._hirta="Paramuriceidae",OCT072_Euplexaura_rhipidalis="Euplexauridae",ANT106_Discophyton_rudyi="Discophytidae",ANT114_Isis_hippuris="Isididae",OCT007_Stenogorgia_casta="Pterogorgiidae",
                                             OCT073_Muricea_atlantica="Plexauridae",ANT102_Paralemnalia_thyrsoides="Lemnaliidae",ANT104_Ovabunda_macrospiculata="Xeniidae",ANT25_Coelogorgia_palmosa="Coelogorgiidae",ANT69_Cladiella_pachyclados="Cladiellidae",
                                             ANT65_Annella_mollis="Subergorgiidae",ANT105_Acrophytum_claviger="Acrophytidae",ANT39_Leptophyton_benayahui="Leptophytidae",ANT40_Corymbophyton_bruuni="Corymbophytidae",OCT140_Pseudonephthea_liltvedi_n._comb="Pseudonephtheidae",
                                             OCT147_Skamnarium_complanatum="Skamnariidae",OCT046_Trachythela_rudis="Victorgorgiidae",ANT44_Dendronephthya_sinaiensis="Nephtheidae",OCT151_Pseudodrifa_racemosa="Capnellidae",ANT107_Ceeceenus_quadrus="Paralcyoniidae",
                                             ANT64_Sinularia_brassica="Sinulariidae",OCT148_Elbeenus_lauramartinae="Aquaumbridae",OCT040_Altumia_delicata="Acanthoaxiidae",OCT042_Cryptophyton_goddardi="Carijoidae",D449_Sclerophytum_humilis_n._comb.="Sarcophytidae",
                                             ANT101_Malacacanthus_capensis="Malacacanthidae",OCT143_Rosgorgia_inexspectata="Rosgorgiidae",OCT043_Inconstantia_procera="Incrustatidae",OCT062_Nidalia_cf._occidentalis="Nidaliidae",ANT53_Siphonogorgia_godeffroyi="Siphonogorgiidae",
                                             ANT108_Azoriella_bayeri="Cerveridae",ANT99_Arula_petunia="Arulidae",ANT24_Clavularia_inflata="Clavulariidae",ANT276_Homophyton_verrucosum="Spongiodermidae",ANT51_Briareum_hamrum="Briareidae",
                                             ANT61_Ideogorgia_capensis="Ideogorgiidae",ANT27_Erythropodium_caribaeorum="Erythropodiidae",ANT30_Parasphaerasclera_valdiviae="Parasphaerascleridae",ANT82_Chelidonisis_aurantiaca_mexicana="Chelidonisididae",ANT46_Trichogorgia_capensis="Ifalukellidae",
                                             OCT022_Acanella_arbuscula="Keratoisididae",ANT170_Isidoides_armata="Isidoidae",ANT79_Nanipora_cf._kamurai="Helioporidae",ANT75_Dendrobrachia_sarmentosa="Dendrobrachiidae",ANT111_Renilla_muelleri="Pennatuloidea",ANT26_Cornularia_pabloi="Cornulariidae",
                                             ALC43_Alcyonium_coralloides="Alcyoniidae",OCT051_Plumigorgia_sp.="Pleurogorgiidae",ANT73_Acrossota_amboinensis="Acrossotidae", Pacifigorgia_irene="Gorgoniidae"))



# We check in one of them that the new tags are there and the old are not
"Ellisella elongata" %in% my_vec_new1prun
'Ellisellidae' %in% my_vec_new1prun

# We substitute the old labes by the new vector
finprun<-tree1.noPR
finprun$tip.label <- my_vec_new1prun

# We add a new column with grouped levels to the dataset
OTD10 <- OTD9 %>%
  mutate(superfamily = ifelse(Family_molecules %in% levels_to_group, "Pennatuloidea", Family_molecules))

# We calculate the number of species per family

N_sp_fam<-aggregate(data = OTD10,                # Applying aggregate
                    specie_name ~ superfamily,
                    function(specie_name) length(unique(specie_name)))

# We calculate for each trait, the number of species per family that have data. Since we are not interested in contextual traits, we removed them first

OTD11<-OTD10[OTD10$trait_category!="Contextual",]

Traitcounts<-aggregate(data = OTD11,                # Applying aggregate
                       specie_name ~ trait_name+superfamily+trait_category,
                       function(specie_name) length(unique(specie_name)))

# We group families of the count file, by creating a new column:
# We also merge the dfs so we can calculate the % of species with data for a given family

Traitcounts1 <- merge(Traitcounts, N_sp_fam, by = "superfamily")
Traitcounts1$percentage <- Traitcounts1$specie_name.x/Traitcounts1$specie_name.y*100

# We also need to order Trait count superfamily accordingly to display it in the tree
# and to include a dummy value called Hexacorallia, so it will match the tree
Traitcounts2<-Traitcounts1
Traitcounts2[nrow(Traitcounts2)+1,]<-c("Hexacorallia","Coloniality","Morphological",0,1,0)
Traitcounts2$percentage=as.numeric(Traitcounts2$percentage)

# We also need to add a summary column stating the number of species with info for a given trait, among the total number of species in a given family
heatx <-Traitcounts2
heatx$trait_category=as.factor(heatx$trait_category)
heatx$superfamily=as.factor(heatx$superfamily)
heatx$trait_name=as.factor(heatx$trait_name)
summary<- heatx
summary$specie_name.x=as.numeric(summary$specie_name.x)
summary$specie_name.y=as.numeric(summary$specie_name.y)
data_sum1 <- aggregate(specie_name.x ~ trait_name+trait_category, summary, sum)       
data_sum1$specie_name.y=sum(N_sp_fam$specie_name)

# Some traits are only possible in sea pens. Thus, we adjust the total maximum number of species to the N of sea pens (234), in those traits
data_sum1 <- data_sum1 %>%
  mutate(specie_name.y = ifelse(trait_name == "Symmetry of rachis", 234, specie_name.y))

data_sum1 <- data_sum1 %>%
  mutate(specie_name.y = ifelse(trait_name == "Peduncle length", 234, specie_name.y))

data_sum1 <- data_sum1 %>%
  mutate(specie_name.y = ifelse(trait_name == "Rachis length", 234, specie_name.y))

data_sum1 <- data_sum1 %>%
  mutate(specie_name.y = ifelse(trait_name == "Presence of polyp leaves", 234, specie_name.y))

data_sum1 <- data_sum1 %>%
  mutate(specie_name.y = ifelse(trait_name == "Polyps per polyp leaf", 234, specie_name.y))

data_sum1$percentage<-data_sum1$specie_name.x/data_sum1$specie_name.y*100
data_sum1$superfamily="All"

heatx1<-rbind(heatx,data_sum1)

# We create our tree for Figure 3
treeq <- groupClade(finprun, c(116,71))

# Since we also want to display the number of species per family in the labels of the figure, 
# we create a new column that display both name of the family and number of species
heatx1$family_combined <-paste(heatx1$superfamily, " (",heatx1$specie_name.y, ") ", sep = "")

# And replace the old tip.labels of the tree but the new ones containing both family name and species number
treeq$tip.label <- heatx1$family_combined[match(treeq$tip.label, heatx1$superfamily)]

#We also remove the superfamily column in the df, since we wont need it anymore 

heatx1 <- subset(heatx1, select = -superfamily )

#Now we can create the tree

p1heat <- ggtree(treeq, aes(color=group),branch.length = "none") +
scale_color_manual(values=c("grey","black","black"))+
  theme(legend.position = "none")+
  #geom_tiplab(align=FALSE,linesize=0) +
  scale_y_reverse()+
  coord_flip()+
  hexpand(.01)

# We also want to introduce the corresponding NAs for those sea-pen realted traits in non-sea pens families

specific_traits <- c("Symmetry of rachis", "Peduncle length","Rachis length","Polyp leaves presence","Polyps per polyp leaf")  # Add the specific traits you want
unique_combinations <- expand.grid(family_combined = unique(heatx1$family_combined), trait_name = specific_traits)
unique_combinations$trait_category = "Morphological"
unique_combinations$specie_name.x=NA
unique_combinations$specie_name.y=NA
unique_combinations$percentage=NA
unique_combinations<-unique_combinations%>%
  filter(family_combined != "Pennatuloidea (234) ")
unique_combinations<-unique_combinations%>%
  filter(family_combined != "Hexacorallia (1) ")
unique_combinations<-unique_combinations%>%
  filter(family_combined != "All (3497) ")
unique_combinations<-unique_combinations%>%
  filter(family_combined != "All (234) ")

heatx2= rbind(heatx1,unique_combinations)

# We plot the combined figure: heatmap + tree

p2b<-ggplot(heatx2, aes(family_combined, reorder(trait_name, desc(trait_name)))) +
  geom_tile(aes(fill = percentage), colour = "white") +
  labs(x = "Family", y = "Trait name")+
  #scale_y_discrete(limits = unique(rev(heatx1$trait_name)))+
  facet_grid(vars(trait_category), space="free", scales="free", switch="y") +
  #scale_x_discrete(limits = unique(heatx$family))+
  scale_fill_gradient(low = "white", high = "darkorchid4", na.value="#f3e592")+
  theme_classic()


p3b<-p2b+ theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))

p3b %>% 
  insert_bottom(p1heat, height=0.075) 


# And save it
ggsave("Figure_3.pdf", width = 30,  height = 35,  units = "cm", limitsize = FALSE)

#### Figures 4 and 5

#We load some extra libraries that we will need.
library(tidyverse)
library(tidyr)
library(dplyr)
require(sf)
library(rnaturalearth)
library("meowR")
library(scatterpie)

# We create the world object
world <- map_data("world")

# We remove trait observations without location, such as the case of global estimates
locations<- OTD
locations1<-locations %>% drop_na(latitude)

# We convert all depths to m units
locations2<-locations1[locations1$trait_name=="Water depth",]
locations2$standard_unit=as.factor(locations2$standard_unit)
locations2$value=as.numeric(locations2$value)

locations2 <- locations2 %>%
  mutate(standardized_depth = case_when(
    standard_unit == "m" ~ value,
    standard_unit == "cm" ~ value / 100,
    standard_unit == "ft" ~ value * 0.3048,  # Assuming 1 ft = 0.3048 meters
    standard_unit == "fm" ~ value * 1.8288, 
  ))

# We associate each standardized depth trait measurement with a corresponding "light or depth zone"
# light zones: photic (0-30 m / mesophotic 31 to 150 m / aphotic >150 m))
locations2 <- locations2 %>%
  mutate(light_zone = case_when(
    standardized_depth <= 30 ~ "photic",
    standardized_depth > 30 & standardized_depth <= 150 ~ "mesophotic",
    standardized_depth > 150 ~ "aphotic",
  ))

# We add marine realms / provinces
meow <- sf::read_sf("./MEOW/meow_ecos.shp")

# World map
world_map <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))

# We count observations per Realm
data("regions")
regionalData<-getRegions(locations1$latitude,locations1$longitude)
newdata <-cbind(locations1,regionalData)

#Base map 
kk2 <- ggplot() +
  geom_sf(data = world_map, size = 0.2, fill = "grey95", col = "grey85") +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8)) +
  coord_sf(expand = FALSE)

# Figure_ MAP (panel b)
newdata_viridis<-newdata
newdata_viridis$REALM<- factor(newdata_viridis$REALM, levels = c(levels(newdata_viridis$REALM), "Undefined realm"))
newdata_viridis$REALM[is.na(newdata_viridis$REALM)] <- "Undefined realm"

color_vector2 <- c("#8707A6FE", "#482173FF", "#6F9CDB", "#38598CFF", "#2D708EFF",
                   "#25858EFF", "#2BB07FFF","#00BFC4", "#51C56AFF", "#85D54AFF",
                   "#C2DF23FF", "#FDE725FF")

color_vector <- c("#8707A6FE", "#482173FF", "#6F9CDB", "#38598CFF", "#2D708EFF",
                  "#25858EFF", "#2BB07FFF","#00BFC4", "#51C56AFF", "#85D54AFF",
                  "#C2DF23FF", "#FDE725FF", "darkgrey")
kk2 +  
  geom_sf(data = meow, aes(fill = REALM), size = 0.2, col = 0, alpha = 0.1) +
  geom_point(
    data = newdata_viridis,
    aes(longitude, latitude, color = REALM, shape=REALM),
    alpha = 0.7, size = 0.5) +
  scale_fill_manual(values=color_vector2)+
  scale_color_manual(values=color_vector)+
  scale_shape_manual(values=c(19,19,19,19,19,19,19,19,19,19,19,19,4))+
  coord_sf(expand = FALSE) +
  theme(legend.position = "none")

# We save it
ggsave("Figure4_map.pdf", width = 15,  height =10,  units = "cm", limitsize = FALSE)

##Figure 4c (Barplot panel)

# We calculate the number of trait observations per Realm, adding also the level factor undefined realm.
newdata_bar <-newdata
newdata_bar$REALM<- factor(newdata_bar$REALM, levels = c(levels(newdata_bar$REALM), "Undefined realm"))
newdata_bar$REALM[is.na(newdata_bar$REALM)] <- "Undefined realm"
total_number_obs_Realm <- aggregate(data = newdata_bar,                # Applying aggregate
                                             observation_id ~ REALM,
                                             function(observation_id) length(unique(observation_id)))

total_number_obs_Realm <- total_number_obs_Realm %>%
  dplyr::rename("Total_N_of_trait_observations" = observation_id)

# We will only display the first 1500 trait observations per Realm. 
# Temperate Northern Atlantic (> 1500 observations)
total_number_obs_Realm_limited <- total_number_obs_Realm
total_number_obs_Realm_limited[6, "Total_N_of_trait_observations"] <- 1500

# We plot it

barplot_realms_v2<-ggplot(total_number_obs_Realm_limited, aes(x=factor(REALM, levels = rev(levels(REALM))), y=Total_N_of_trait_observations, fill=REALM)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=color_vector)+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  scale_y_continuous(limits = c(0,1500), expand = c(0, 0)) +
  theme(legend.position = "right")+
  coord_flip()

# We save it
ggsave("Figure_4_barplot.pdf", width = 10,  height =10,  units = "cm", limitsize = FALSE)

# Figure 4 - Panel a

# N of trait observations that are Global estimates 

locations_global <- locations[locations$location_name=="Global estimate",]

glob<-length(unique(locations_global$observation_id))
#87368

#if we consider regional estimates within global estimates
locations_globalv2<- locations[is.na(locations$latitude), ]
glob2<-length(unique(locations_globalv2$observation_id))
#87426


# N of trait observations that have geographical coordinates.

geogr<-length(unique(locations1$observation_id))
# 10325

#We calculate the percentages

per_global <- glob2/(glob2+geogr)*100
per_geogr <- geogr/(glob2+geogr)*100

#per_global_89.44
#per_geogr:10.56

# We create the plot
percentage <- data.frame(data_type = c("Global estimate", "Geographically represented"), 
                       Proportion = c(per_global,per_geogr))

percentage <- percentage %>% 
  mutate(Representation = "data types")

ggplot(percentage, aes(x = Proportion, y = Representation, fill = data_type)) +
  geom_col() +
  scale_fill_manual(values=c("Global estimate"="black","Geographically represented"="grey")) +
  theme_minimal(base_size = 16) +
  ylab("Percentage") +
  theme_void()+
  theme(legend.position = "none")+
  xlab(NULL)

# We save it
ggsave("Figure_4_panel_a.pdf", width = 8,  height =1.5,  units = "cm", limitsize = FALSE)

# Figure 5a. Map

# We upload a file with the coordinates of the marine provinces. 

provinces_coord<-read.csv("provinces_coord.csv")

# Now we need to count how many observations are in each light zone per province
regionalData1<-getRegions(locations2$latitude,locations2$longitude)
newdata1 <-cbind(locations2,regionalData1)
count_light_province <- aggregate(data = newdata1,                # Applying aggregate
                                  observation_id ~ PROVINCE+light_zone,
                                  function(observation_id) length(unique(observation_id)))

# We combine the df with the N of trait observations with the coordinates one
newdata2b<-merge(provinces_coord, count_light_province, by="PROVINCE")

# We calculate the number of total observations per province, for the size of the pies
total_number_light_obs_province <- aggregate(data = newdata1,                # Applying aggregate
                                             observation_id ~ PROVINCE,
                                             function(observation_id) length(unique(observation_id)))
# We combine both df again

newdata2c<-merge(newdata2b, total_number_light_obs_province, by="PROVINCE")


# We we need the df in horizontal format and calculate the percentages
df_transformed2<-spread(newdata2c, light_zone,observation_id.x )
df_transformed2 <- df_transformed2 %>% replace_na(list("aphotic" = 0, "mesophotic"=0,"photic" = 0))
df_transformed2$aphotic_percentage <- df_transformed2$aphotic/df_transformed2$observation_id.y*100
df_transformed2$mesophotic_percentage <- df_transformed2$mesophotic/df_transformed2$observation_id.y*100
df_transformed2$photic_percentage <- df_transformed2$photic/df_transformed2$observation_id.y*100

# We associate colors to the light zones
colors2 <- c("aphotic_percentage" = "#03063E", 
            "mesophotic_percentage" = "#006C82", 
            "photic_percentage" = "#8AD2E9")

# We also need to prepare our radius legend
radius <- log(c(5,50,500,5000))
labeller <- function(x) round(exp(x))
breaks <- scales::breaks_extended(n = 4)(radius)

# We plot it

ggplot() +
  geom_sf(data = world_map, size = 0.2, fill = "gray90", col = "gray90") +
  theme(panel.background = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8)) +
  coord_sf(expand = FALSE) +
  geom_scatterpie(aes(x = longitude, y = latitude, r = log(observation_id.y), fill = factor(cols)),
                  data = df_transformed2, color=NA,
                  cols = c("aphotic_percentage", "mesophotic_percentage", "photic_percentage"))+
  scale_fill_manual(values = colors2) +
   geom_scatterpie_legend(breaks = radius,
                          x = -150.54, y = -50, n = 2,
                         labeller = labeller)+
  theme(legend.position = "none")

# We save it
ggsave("Figure_5_Map.pdf", width = 15,  height =10,  units = "cm", limitsize = FALSE)

# Figure 5b
# We calculate the mean percentage (and SD) of trait observations across all light zones
light_bar <-newdata2c
light_bar<-tidyr::complete(light_bar, light_zone, PROVINCE)
light_bar$percentage <- light_bar$observation_id.x/light_bar$observation_id.y*100
light_bar[is.na(light_bar)] <- 0

light_mean_percentages <- aggregate(percentage ~ light_zone, data=light_bar, FUN = mean)
colnames(light_mean_percentages)[2] <- "mean"
light_sd_percentages<-aggregate(percentage ~ light_zone, data=light_bar, FUN=sd) 
colnames(light_sd_percentages)[2] <- "sd"
light_mean_sd_percentages <- merge.data.frame(light_mean_percentages, light_sd_percentages, all.x=TRUE)

# We plot it
barplot_lights<-ggplot(light_mean_sd_percentages, aes(light_zone, y=mean, fill=light_zone)) +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), position = 'dodge', width = 0.2) +
  geom_bar(stat='identity', position="dodge") +
  scale_fill_manual(values=c("#03063E", "#006C82","#8AD2E9"))+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8))+
  #scale_y_continuous(limits = c(0,1500), expand = c(0, 0)) +
  theme(legend.position = "right")+
  coord_flip()

barplot_lights
#We save it
ggsave("Figure 5_barplot.pdf", width = 10,  height =10,  units = "cm", limitsize = FALSE)

###### Supplementary figures
# Figure S1 - N of traits per species

# We calculate the number of traits per species
unique_traits <- aggregate(trait_name ~ specie_name+Family_molecules, data = OTD9, FUN = function(x) length(unique(x)))
unique_traits1 <- unique_traits %>%
  mutate(superfamily = ifelse(Family_molecules %in% levels_to_group, "Pennatuloidea", Family_molecules))

unique_traits_ordered <-unique_traits1[order(unique_traits1$trait_name), ]

# We color only large families (>75 species), as we did with Figure 2. 
unique_traits1_short <- unique_traits1 %>%
  mutate(grouped_family = ifelse(superfamily %in% maintained_levels_short, superfamily, "0_other"))

# We plot it
p_supp_unique <- ggtree(tree31_short, aes(color=group),branch.length = "none",layout="circular", size=0.1) +
  scale_color_manual(values=c("#979198","black"))+
  #geom_tiplab(align=FALSE,linesize=0) +
  hexpand(.01)

p_supp_unique+ geom_fruit(data=unique_traits1_short,
                      axis.params=list(axis="x",text.size=2,limits=c(0,90),title="N of traits",title.height = 0.05, title.size=2),
                      grid.params=list(vline=FALSE,size=0.15,color="lightgrey"),
                      geom=geom_bar,
                      mapping=aes(y=specie_name, x=trait_name, fill=grouped_family),
                      pwidth=1.2,
                      stat="identity",
                      orientation="y") +
  scale_fill_manual(values=c("lightgrey","#EDF8B1", "#BD0026", "#41B6C4", "#c36721", "#F3d643",
                             "#253494", "#FDEBB6", "#E31A1C", "#7FCDBB", "#081D58",
                             "#FEB24C", "#225EA8", "#800026", "#C7E9B4" ,"#1D91C0"),
                    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)) +
  theme(legend.justification = "right",
        legend.title=element_text(size=7),
        legend.text=element_text(size=6),
        legend.spacing.y = unit(0.02, "cm"))+
  guides(color = "none")

# We save it
ggsave("Figure_S1.pdf", width = 8,  height =8,  units = "in", limitsize = FALSE)

#### Figure S2
# We select only uncertain families
OTD1b<-OTD[OTD$Family_molecules=="Malacalcyonacea_incertae_sedis",]
OTD2b<-OTD[OTD$Family_molecules=="Scleralcyonacea_incertae_sedis",]
OTD3b<-OTD[OTD$Family_molecules=="incertae_sedis",]
OTD4b<-OTD[OTD$Family_molecules=="uncertain",]
OTD5b<-OTD[OTD$Family_molecules=="Haimeidae",]
OTDB<-rbind(OTD1b,OTD2b,OTD3b,OTD4b,OTD5b)
OTDB$genus <- sapply(strsplit(OTDB$specie_name, " "), function(x) x[1])

# We calculate the number of species per genus

N_sp_genus<-aggregate(data = OTDB,                
                    specie_name ~ genus,
                    function(specie_name) length(unique(specie_name)))

# We calculate, for each trait, the number of species per genus that have data. Since we are not interested in contextual traits, we removed them first

OTD11B<-OTDB[OTDB$trait_category!="Contextual",]

Traitcounts_genus<-aggregate(data = OTD11B,                # Applying aggregate
                       specie_name ~ trait_name+genus+trait_category,
                       function(specie_name) length(unique(specie_name)))

# We group families of the count file, by creating a new column, and merge the dfs so we can calculate the % of species with data for a given family
Traitcounts1B <- merge(Traitcounts_genus, N_sp_genus, by = "genus")
Traitcounts1B$percentage <- Traitcounts1B$specie_name.x/Traitcounts1B$specie_name.y*100

# We order Trait count superfamily accordingly to display it in the tree
Traitcounts1BB <- Traitcounts1B[order(Traitcounts1B$trait_category),]
heatxb <-Traitcounts1BB
heatxb$trait_category=as.factor(heatxb$trait_category)
heatxb$genus=as.factor(heatxb$genus)
heatxb$trait_name=as.factor(heatxb$trait_name)

# We also add a summary column stating the number of species with info for a given trait, among the total number of species
summaryb<- heatxb
summaryb$specie_name.x=as.numeric(summaryb$specie_name.x)
summaryb$specie_name.y=as.numeric(summaryb$specie_name.y)
data_sum1b <- aggregate(specie_name.x ~ trait_name+trait_category, summaryb, sum)       
data_sum1b$specie_name.y=sum(N_sp_genus$specie_name)

# Since we also want to display the number of species per genus in the labels of the figure, 
# we create a new column that display both name of the genus and number of species
heatxb$genus_combined <-paste(heatxb$genus, " (",heatxb$specie_name.y, ") ", sep = "")

# We remove the genus column in the df, since we wont need it anymore 
heatxb <- subset(heatxb, select = -genus )

# Now we need to add missing traits
# Get unique trait names from df1
all_traits_df1 <- unique(heatx1[c("trait_name", "trait_category")])

# Get unique trait names from df2
all_traits_df2 <- unique(heatxb[c("trait_name", "trait_category")])

# Identify missing trait names
missing_traits <- setdiff(all_traits_df1$trait_name, all_traits_df2$trait_name)

# Create a data frame of missing traits with their associated category
missing_rows <- data.frame(
  trait_name = missing_traits,
  trait_category = all_traits_df1$trait_category[match(missing_traits, all_traits_df1$trait_name)],
  stringsAsFactors = FALSE
)

# Get unique specie names from df2
all_species <- unique(heatxb$genus)

# Repeat each row of missing_rows for each specie name
missing_rows <- do.call(rbind, lapply(all_species, function(specie) {
  missing_rows$genus_combined <- specie
  missing_rows
}))

# Add n_observations as 0
missing_rows$n_observations <- 0
names(missing_rows)[names(missing_rows) == "n_observations"] <- "percentage"

# remove species_x and y from heatxb
heatxc <- heatxb[, !(names(heatxb) %in% c("specie_name.x","specie_name.y"))]

# we bind both df
heatxd<-rbind(heatxc,missing_rows)
heatxe<-heatxd[order(heatxd$trait_name), ]
heatxd <- heatxd[order(heatxd$trait_name),]
heatxd$trait_name <- as.character(heatxd$trait_name)
heatxd <- heatxd[order(heatxd$trait_category, heatxd$trait_name), , drop = FALSE]
heatxd$trait_name<-as.factor(heatxd$trait_name)

# Now we plot the figure

p2b_supp<-ggplot(heatxd, aes(genus_combined, reorder(trait_name, desc(trait_name)))) +
  geom_tile(aes(fill = percentage), colour = "white") +
  labs(x = "Genus", y = "Trait name")+
  #scale_y_discrete(limits = unique(rev(heatx1$trait_name)))+
  facet_grid(vars(trait_category), space="free", scales="free", switch="y") +
  #scale_x_discrete(limits = unique(heatx$family))+
  scale_fill_gradient(low = "white", high = "darkorchid4")+
  theme(axis.text.x = element_text(face = "italic"))+
  theme_classic()

p3b_supp<-p2b_supp+ theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))

# We save it
ggsave("Figure_S2.pdf", width = 25,  height =35,  units = "cm", limitsize = FALSE)

##END OF THE CODE

