##############################################
############# PREPARATION  ###################
##############################################

#Choix du dossier de travail, devant intégrer les fichiers de données, et où les fichiers de sortie seront créés
setwd("C:/Users/a.queffelec/Documents/Analyses/Networks/Aurignacien")

# Lancement des packages nécessaires
library(statnet)
library(tnet)

#######################################################
############# IMPORTATION DES DONNEES  ################
#######################################################


# Importer les données du fichier Data.csv (séparation point-virgule) dans un objet nommé Aurig.
Aurig <- read.csv(file = "Data_quanti.csv", sep=";", row.names = 1)
# Import des attributs des ensembles archéologiques à partir du fichier csv Attributes.csv
Aurig.attr <- read.csv("Attributes.csv", sep=";", row.names = 1)

#####################################################
######    CALCUL DE LA MATRICE DE SIMILARITE   ######
#####################################################

### Creation de la matrice de similarité de co-présence (1 ou 0 si les ensembles archéologiques définis partagent le critère)
co.p <- function(x) {
  x <- as.matrix(x)
  # redéfinit comme présent (1) quelque soit la valeur dans le tableau
  x[x >= 1] <- 1
  # redéfinit tout le reste comme absent (0)
  x[x < 1] <- 0
  # Calcul matriciel pour trouver les co-presence (%*% indique une multiplication matricielle), t(x) est la transposée de de la matrice x
  out <- x %*% t(x)
  return(out)
}
# Appliquer cette fonction co.p à nos données
AurigP <- co.p(Aurig)

####################################################################################
#############      CREATION DES RESEAUX PONDERES EN FONCTION          ##############
##########   DU SEUIL ET CREATION DES FIGURES REGROUPEES DANS UN PDF    ############
####################################################################################

### Création du fichier pdf
namepdf=paste0("Réseau pondéré de co-présence avec seuil variable et centralité degré.pdf")
pdf(namepdf,height=40,width=20)

### Création de la grille permettant la mémorisation de plusieurs graphiques les uns à côté des autres
par(mfrow = c(7, 2))

### Création nécessaire d'une variable en dehors de la boucle
Temp_AurigP = AurigP

### Boucle de réalisation des réseaux et graphiques en fonction du seuil variable T 
### représentant la force des liens entre les noeuds du réseau
for (T in seq(0,6)) {

    ### Calcul du réseau en ne gardant que les liens dont le poids est supérieur ou égal au seuil  
    Pnet2 <- network(event2dichot(AurigP, method = "absolute", thresh = 2), 
                 directed = F,ignore.eval = F, names.eval = "weight") # Ceci ne conserve que les liens ayant une force supérieur ou égale au seuil T mais en leur donnant une valeur de 1
    set.edge.value(Pnet2, "weight", AurigP) # Ceci redonne leur valeur initiale aux liens, mais seulement à ceux qui avaient été conservés
    weight = get.edge.value(Pnet2, "weight") # Ceci crée une variable qui sera utilisée pour donner l'épaisseur et la couleur des liens dans le graphique
    Pnet2 %v% "vertex.names" <- row.names(AurigP) # Ceci donne les noms des ensembles archéologiques aux noeuds du réseau

    ### Calcul de la centralité des noeuds
    Temp_AurigP[Temp_AurigP <= T] <- 0 # on ne garde que les liens dont le poids est supérieurs ou égal au seuil
    dg.wt <- as.matrix(rowSums(Temp_AurigP) - 1) # la centralité d'un noeud est la somme des poids des liens qui lui sont connectés (degree centrality)
    dg.wt[dg.wt == -1] <- 1 # valeur toute petite de 1 même quand il n'y a pas de lien pour que le noeud soit tout de même visible dans le graphique

    # Création des variables pour l'attribution des couleurs des noeuds et des liens
    reg.col <- as.factor(Aurig.attr$Groupe)
    edge.cols <- colorRampPalette(c("blue", "yellow"))(7)

    # Réaliser le graphique en utilisant la disposition générique
    plot(Pnet2, 
     #edge.col = "gray",
     edge.lwd = 2*weight,
     edge.col = edge.cols[get.edge.value(Pnet2, "weight")],
     vertex.cex = 0.05*dg.wt,
     vertex.col = reg.col,
     displaylabels = TRUE,
     main = paste("Réseau pondéré de co-présence, seuil =",T))
    legend("left", title = "Force du lien",legend = seq(0, 6), lwd = seq(0, 6), col = edge.cols)
    legend("bottomleft", title = "Groupe", legend = levels(reg.col), pch=16, pt.cex = 2 ,col = as.factor(levels(reg.col)))
    legend("topleft", legend = "Taille du noeuds = centralité degré")

    # Réaliser le graphique en utilisant les coordonnées géographiques
    plot(Pnet2, 
     edge.lwd = 2*weight,
     edge.col = edge.cols[get.edge.value(Pnet2, "weight")], 
     vertex.cex = 0.05*dg.wt,
     coord = cbind(jitter(Aurig.attr[,2],factor=5),jitter(Aurig.attr[,3],factor=40)),
     displaylabels = TRUE,
     vertex.col = reg.col,
     main = paste("Réseau pondéré de co-présence, seuil =",T))
    legend("left", title = "Force du lien",legend = seq(0, 6), lwd = seq(0, 6), col = edge.cols)
    legend("bottomleft", title = "Groupe", legend = levels(reg.col), pch=16, pt.cex = 2 ,col = as.factor(levels(reg.col)))
    legend("topleft", legend = "Taille du noeuds = centralité degré")
}
dev.off()

