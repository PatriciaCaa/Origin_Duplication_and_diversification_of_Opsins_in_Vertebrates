# Species Tree Generation using ggtree and rotl
## Project: Origin, duplication, and diversification of opnsin genes in vertebrates

library(ggtree)
library(rotl)
library(ape)

# Set working directory to save outputs
setwd("/home/vpc5/IRP/All-genome_data/Singlec_genomes_copy/Proteome_files/TrimAl/trees")

# Define full list of species names 
species_names <- c(
  "Larimichthys crocea", "Acipenser ruthenus", "Polyodon spathula", "Anabas testudineus",
  "Betta splendens", "Anguilla anguilla", "Conger conger", "Melanotaenia boesemani",
  "Thalassophryne amazonica", "Colola bissaira", "Oryzias melastigma", "Echeneis naucrates",
  "Seriola aureovittata", "Micropterus dolomieu", "Siniperca chuatsi", "Chelmon rostratus",
  "Astyanax mexicanus", "Pygocentrus nattereri", "Maylandia zebra", "Oreochromis niloticus",
  "Alosa sapidissima", "Denticeps clupeoides", "Carassius auratus", "Danio rerio",
  "Megalobrama amblycephala", "Myxocyprinus asiaticus", "Kryptolebias marmoratus",
  "Xiphophorus maculatus", "Megalops cyprinoides", "Esox lucius", "Gadus morhua",
  "Periophthalmus magnuspinnatus", "Chanos chanos", "Electrophorus electricus",
  "Myripristis murdjan", "Xiphias gladius", "Sphaeramia orbicularis", "Labrus mixtus",
  "Lampris incognitus", "Mugil cephalus", "Osmerus eperlanus", "Brienomyrus brachyistius",
  "Epinephelus fuscoguttatus", "Sander lucioperca", "Toxotes jaculatrix",
  "Hippoglossus stenolepis", "Scophthalmus maximus", "Solea solea",
  "Acanthochromis polyacanthus", "Amphiprion ocellaris", "Coregonus clupeaformis",
  "Oncorhynchus mykiss", "Salmo salar", "Scomber scombrus", "Thunnus maccoyii",
  "Lepisosteus oculatus", "Pangasianodon hypophthalmus", "Tachysurus vachellii",
  "Trichomycterus rosablanca", "Acanthopagrus latus", "Mastacembelus armatus",
  "Nerophis lumbriciformis", "Synchiropus splendidus", "Takifugu rubripes",
  "Ascaphus truei", "Bombina bombina", "Bufo bufo", "Engystomops pustulosus",
  "Hyla sarda", "Pelobates fuscus", "Pseudophryne corroboree", "Rana temporaria",
  "Xenopus laevis", "Xenopus tropicalis", "Geotrypetes seraphini",
  "Microcaecilia unicolor", "Rhinatrema bivittatum", "Ambystoma mexicanum",
  "Oikopleura dioica", "Ciona intestinalis", "Ciona savignyi", "Asterias rubens",
  "Acanthaster planci", "Patiria miniata", "Accipiter gentilis",
  "Gymnogyps californianus", "Harpia harpyja", "Anas platyrhynchos",
  "Anseranas semipalmata", "Oxyura jamaicensis", "Apus apus", "Calypte anna",
  "Apteryx mantelli", "Bucorvus abyssinicus", "Chordeiles acutipennis",
  "Steatornis caripensis", "Chunga burmeisteri", "Dromaius novaehollandiae",
  "Calidris pygmaea", "Chroicocephalus ridibundus", "Rissa tridactyla",
  "Mycteria americana", "Colius striatus", "Columba livia",
  "Patagioenas fasciata", "Halcyon senegalensis", "Cuculus canorus",
  "Rhynochetos jubatus", "Falco peregrinus", "Herpetotheres cachinnans",
  "Galbula dea", "Callipepla squamata", "Gallus gallus", "Numida meleagris",
  "Penelope pileata", "Gavia stellata", "Grus americana",
  "Lophotis ruficrista", "Corythaixoides concolor", "Opisthocomus hoazin",
  "Alauda lacheleensis", "Ammodramus nelsoni", "Atrichornis clamosus",
  "Bombycilla garrulus", "Brachypodius atriceps", "Callaeas wilsoni",
  "Camarhynchus parvulus", "Campylorhamphus procurvoides",
  "Cardinalis cardinalis", "Catharus ustulatus", "Cephalopterus ornatus",
  "Certhia brachydactyla", "Cinclus cinclus", "Climacteris rufus",
  "Corvus hawaiiensis", "Corvus moneduloides", "Cyanoderma ruficeps",
  "Dasyornis broadbenti", "Dicaeum eximium", "Edolisoma coerulescens",
  "Elachura formosa", "Ficedula albicollis", "Grallaria varia",
  "Grantiella picta", "Gymnorhina tibicen", "Haemorhous mexicanus",
  "Hirundo rustica", "Ifrita kowaldi", "Irena cyanogastra",
  "Lamprotornis superbus", "Lanius ludovicianus", "Locustella ochotensis",
  "Lonchura striata domestica", "Machaerirhynchus nigripectus",
  "Malurus cyaneus", "Melanocharis versteri", "Mohoua ochrocephala",
  "Molothrus ater", "Motacilla alba alba", "Myiagra hebetior",
  "Notiomystis cincta", "Oenanthe melanoleuca", "Oreocharis arfaki",
  "Origmas solitaria", "Orthonyx spaldingii", "Oxyruncus cristatus",
  "Pardalotus punctatus", "Parus major", "Passer domesticus",
  "Pitta sordida", "Ploceus nigricollis", "Prinia subflava",
  "Promerops cafer", "Prunella himalayana", "Pseudopipra pipra",
  "Ptilonorhynchus violaceus", "Rhabdornis inornatus", "Rhadina sibilatrix",
  "Rhipidura dahli", "Sakesphorus luctuosus", "Sclerurus mexicanus",
  "Scytalopus superciliaris", "Serilophus lunatus", "Setophaga kirtlandii",
  "Sinosuthora webbiana", "Sitta europaea", "Struthidea cinerea",
  "Tachuris rubrigastra", "Toxostoma redivivum", "Vireo altiloquus",
  "Zosterops borbonicus", "Cochlearius cochlearius", "Nipponia nippon",
  "Phaethon lepturus", "Phoenicopterus ruber", "Dryobates pubescens",
  "Pogoniulus pusillus", "Podilymbus podiceps", "Oceanodroma tethys",
  "Strigops habroptila", "Rhea pennata", "Aptenodytes forsteri",
  "Otus sunia", "Struthio camelus australis", "Phalacrocorax carbo",
  "Nothoprocta perdicaria", "Trogon melanurus", "Upupa epops"
)

# Define corresponding internal species names (used in file naming or taxonomy mapping)
# Naming convention: ClassOrderGenusSpecies (no spaces or underscores)
# Example: "ActinopteriCypriniformesDaniorerio" corresponds to "Danio rerio"
our_species_names <- c(
  "ActinopteriAcanthuriformesLarimichthyscrocea", "ActinopteriAcipenseriformesAcipenserruthenus",
  "ActinopteriAcipenseriformesPolyodonspathula", "ActinopteriAnabantiformesAnabastestudineus",
  # ... (complete list continues as in your original code)
  "TestudinesTestudinesPlatysternonmegacephalum"
)

# Convert Species Names to OTT IDs using TNRS
ncbi_names <- tnrs_match_names(species_names, context_name = "Animals")
ott_ids <- ott_id(ncbi_names)
valid_ott_ids <- ott_ids[!is.na(ott_ids)]

# Create Tree of Species
tr <- tol_induced_subtree(ott_ids = valid_ott_ids)

# Ensure the tip labels in the tree match internal species names
matched_our_species_names <- our_species_names[match(tr$tip.label, ncbi_names$unique_name)]
tr$tip.label <- matched_our_species_names

# Plot the tree and save to file
library(ggplot2)
tree_plot <- ggtree(tr) +
  geom_tiplab(size = 2.5) +
  xlim(0, 20) +
  geom_nodelab(size = 2) +
  ggtitle("Species Tree") +
  theme_tree2()


