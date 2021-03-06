##############################################
############# PREPARATION  ###################
##############################################

#Choix du dossier de travail, devant int�grer les fichiers de donn�es, et o� les fichiers de sortie seront cr��s
setwd("C:/Users/a.queffelec/Documents/Analyses/Networks/Aurignacien")

# Lancement des packages n�cessaires
library(statnet)
library(tnet)

#######################################################
############# IMPORTATION DES DONNEES  ################
#######################################################


# Importer les donn�es du fichier Data.csv (s�paration point-virgule) dans un objet nomm� Aurig.
Aurig <- read.csv(file = "Data_quanti.csv", sep=";", row.names = 1)
# Import des attributs des ensembles arch�ologiques � partir du fichier csv Attributes.csv
Aurig.attr <- read.csv("Attributes.csv", sep=";", row.names = 1)

#####################################################
######    CALCUL DE LA MATRICE DE SIMILARITE   ######
#####################################################

### Creation de la matrice de similarit� de co-pr�sence (1 ou 0 si les ensembles arch�ologiques d�finis partagent le crit�re)
co.p <- function(x) {
  x <- as.matrix(x)
  # red�finit comme pr�sent (1) quelque soit la valeur dans le tableau
  x[x >= 1] <- 1
  # red�finit tout le reste comme absent (0)
  x[x < 1] <- 0
  # Calcul matriciel pour trouver les co-presence (%*% indique une multiplication matricielle), t(x) est la transpos�e de de la matrice x
  out <- x %*% t(x)
  return(out)
}
# Appliquer cette fonction co.p � nos donn�es
AurigP <- co.p(Aurig)

####################################################################################
#############      CREATION DES RESEAUX PONDERES EN FONCTION          ##############
##########   DU SEUIL ET CREATION DES FIGURES REGROUPEES DANS UN PDF    ############
####################################################################################

### Cr�ation du fichier pdf
namepdf=paste0("R�seau pond�r� de co-pr�sence avec seuil variable et centralit� degr�.pdf")
pdf(namepdf,height=40,width=20)

### Cr�ation de la grille permettant la m�morisation de plusieurs graphiques les uns � c�t� des autres
par(mfrow = c(7, 2))

### Cr�ation n�cessaire d'une variable en dehors de la boucle
Temp_AurigP = AurigP

### Boucle de r�alisation des r�seaux et graphiques en fonction du seuil variable T 
### repr�sentant la force des liens entre les noeuds du r�seau
for (T in seq(0,6)) {

    ### Calcul du r�seau en ne gardant que les liens dont le poids est sup�rieur ou �gal au seuil  
    Pnet2 <- network(event2dichot(AurigP, method = "absolute", thresh = 2), 
                 directed = F,ignore.eval = F, names.eval = "weight") # Ceci ne conserve que les liens ayant une force sup�rieur ou �gale au seuil T mais en leur donnant une valeur de 1
    set.edge.value(Pnet2, "weight", AurigP) # Ceci redonne leur valeur initiale aux liens, mais seulement � ceux qui avaient �t� conserv�s
    weight = get.edge.value(Pnet2, "weight") # Ceci cr�e une variable qui sera utilis�e pour donner l'�paisseur et la couleur des liens dans le graphique
    Pnet2 %v% "vertex.names" <- row.names(AurigP) # Ceci donne les noms des ensembles arch�ologiques aux noeuds du r�seau

    ### Calcul de la centralit� des noeuds
    Temp_AurigP[Temp_AurigP <= T] <- 0 # on ne garde que les liens dont le poids est sup�rieurs ou �gal au seuil
    dg.wt <- as.matrix(rowSums(Temp_AurigP) - 1) # la centralit� d'un noeud est la somme des poids des liens qui lui sont connect�s (degree centrality)
    dg.wt[dg.wt == -1] <- 1 # valeur toute petite de 1 m�me quand il n'y a pas de lien pour que le noeud soit tout de m�me visible dans le graphique

    # Cr�ation des variables pour l'attribution des couleurs des noeuds et des liens
    reg.col <- as.factor(Aurig.attr$Groupe)
    edge.cols <- colorRampPalette(c("blue", "yellow"))(7)

    # R�aliser le graphique en utilisant la disposition g�n�rique
    plot(Pnet2, 
     #edge.col = "gray",
     edge.lwd = 2*weight,
     edge.col = edge.cols[get.edge.value(Pnet2, "weight")],
     vertex.cex = 0.05*dg.wt,
     vertex.col = reg.col,
     displaylabels = TRUE,
     main = paste("R�seau pond�r� de co-pr�sence, seuil =",T))
    legend("left", title = "Force du lien",legend = seq(0, 6), lwd = seq(0, 6), col = edge.cols)
    legend("bottomleft", title = "Groupe", legend = levels(reg.col), pch=16, pt.cex = 2 ,col = as.factor(levels(reg.col)))
    legend("topleft", legend = "Taille du noeuds = centralit� degr�")

    # R�aliser le graphique en utilisant les coordonn�es g�ographiques
    plot(Pnet2, 
     edge.lwd = 2*weight,
     edge.col = edge.cols[get.edge.value(Pnet2, "weight")], 
     vertex.cex = 0.05*dg.wt,
     coord = cbind(jitter(Aurig.attr[,2],factor=5),jitter(Aurig.attr[,3],factor=40)),
     displaylabels = TRUE,
     vertex.col = reg.col,
     main = paste("R�seau pond�r� de co-pr�sence, seuil =",T))
    legend("left", title = "Force du lien",legend = seq(0, 6), lwd = seq(0, 6), col = edge.cols)
    legend("bottomleft", title = "Groupe", legend = levels(reg.col), pch=16, pt.cex = 2 ,col = as.factor(levels(reg.col)))
    legend("topleft", legend = "Taille du noeuds = centralit� degr�")
}
dev.off()

