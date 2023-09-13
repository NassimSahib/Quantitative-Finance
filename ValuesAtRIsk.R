#Code modélisation de la VaR par la Théorie des Valeurs Extrêmes (TVE) sur l'indice FTSE MIB


#On importe les librairies nécessaires
library(quantmod)
library(DistributionUtils)
library(ggplot2)
library(moments)
library(rugarch)
library(LSTS)
library(forecast)
library(aTSA)



# Téléchargement des données
getSymbols("FTSEMIB.MI",src="yahoo",auto.assign=TRUE)

#Sauvegarde
save(FTSEMIB.MI,file = "FTSEMIB.rdata")

# Traitement des variables manquantes
sum(is.na(FTSEMIB.MI))
FTSEMIB.MI<-na.fill(FTSEMIB.MI,fill = "extend")
sum(is.na(FTSEMIB.MI))

#Graphique de l'indice et création des rendements du FTSE MIB à partir de l'indice ajusté

FTSEMIB.Index<-FTSEMIB.MI$FTSEMIB.MI.Adjusted
X11()
plot(as.zoo(FTSEMIB.Index))

rend_FTSEMIB<-na.omit(diff(log(FTSEMIB.Index)))
length(rend_FTSEMIB)
summary(rend_FTSEMIB)
X11()
plot(as.zoo(rend_FTSEMIB))

#calcul des moments statistiques de distributions de probabilité

kurtosis(as.numeric(rend_FTSEMIB), na.rm = TRUE)
skewness(as.numeric(rend_FTSEMIB), na.rm = TRUE)

#Affichage de l'histogramme des rendements / densité de répartition comparé à la loi normale
X11()
hist(rend_FTSEMIB, breaks=100,xlim = c(-0.14,0.10))
hist(rend_FTSEMIB, breaks=100,xlim = c(-0.14,0.10),prob="TRUE")

lines(density(rend_FTSEMIB, na.rm = "TRUE"),lwd=2,col="orange")
lines(density(rnorm(nrow(rend_FTSEMIB),mean=mean(rend_FTSEMIB,na.rm = TRUE),sd=sd(rend_FTSEMIB,na.rm = TRUE))),lwd=2,col="red")


# QQ Plot pour vérifier la normalité des rendements
X11()
qqnorm((rend_FTSEMIB-mean(rend_FTSEMIB,na.rm = TRUE))/sd(rend_FTSEMIB,na.rm = TRUE));qqline((rend_FTSEMIB-mean(rend_FTSEMIB,na.rm = TRUE))/sd(rend_FTSEMIB,na.rm = TRUE))

#Tests sur la normalité des rendements du FTSE MIB
agostino.test(as.numeric(rend_FTSEMIB)) #H0 : Skewness = 0
anscombe.test(as.numeric(rend_FTSEMIB)) #H0 : Kurtosis = 3
jarque.test(as.numeric(rend_FTSEMIB)) #H0 : la série suit une loi normale


# On montre que la volatilité n'est pas constante:
X11()
vol_glissante<-rollapply(rend_FTSEMIB,width=20,FUN = sd)
plot(as.zoo(vol_glissante),main="volatilité glissante sur une fenetre de 20 jours")
abline(h=as.numeric(sd(rend_FTSEMIB,na.rm = TRUE)),col="red")




## VAR NON CONDITIONNELLE (Processus i.i.d)


#VaR historique

#Représentation graphique à partir d'un histogramme du seuil de 1%
alpha <- 0.01
q <- quantile(as.numeric(rend_FTSEMIB), probs=alpha)
hist(as.numeric(rend_FTSEMIB), breaks=50, freq=FALSE, 
     xlab = "Rendements du FTSE MIB",
     main = "Histogramme du FTSE MIB au seuil de 1%")
abline(v=q, col="red")

#Valeur de la VaR historique à 1 jour au seuil de 1%
VaR_HS<-(-q)
VaR_HS


#Représentation des rendements et de la Var historique avec ses dépassements
X11()
plot(as.zoo(rend_FTSEMIB),
     main = "Rendements et Value at Risk historique",
     xlab = "Temps",
     ylab = "Rendements")
abline(h=-VaR_HS, col="red")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_HS]), pch=20, col="red")
legend("top", legend = "VaR Historique", lty = 1, col = "red")


#VaR variance-covariance / hyp : les rendements suivent une loi normale

#Calcul de la VaR variance-covariance : -moyenne + déviation standard * quantile
VaR_VCov<-(-(mean(rend_FTSEMIB)+sd(rend_FTSEMIB)*qnorm(0.01)))
VaR_VCov

#Représentation des rendements et de la Var variance-covariance avec ses dépassements
X11()
plot(as.zoo(rend_FTSEMIB),
     main = "Rendements et Value at Risk variance-covariance",
     xlab = "Temps",
     ylab = "Rendements")
abline(h= -VaR_VCov, col="blue")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_VCov]), pch=20, col="blue")
legend("top", legend = "VaR Vcov", lty = 1, col = "blue")


#VaR Cornish-Fisher

#Calcul de la VaR Cornish-Fisher : méthode dans le but d'améliorer la précision de la VaR Vcov
K<-kurtosis(rend_FTSEMIB)
S<-skewness(rend_FTSEMIB)

(k<-as.numeric(qnorm(0.01)+(qnorm(0.01)^2-1)*(S/6)+(qnorm(0.01)^3-3*qnorm(0.01))*(K/24)+(2*qnorm(0.01)^3-5*qnorm(0.01))*(S^2/36)))

VaR_CFisher<-as.numeric((-sd(rend_FTSEMIB)*k))
VaR_CFisher


#Représentation des rendements et de la VaR Cornish-Fisher avec ses dépassements
X11()
plot(as.zoo(rend_FTSEMIB),
     main = "Rendements et Value at Risk Cornish-Fisher",
     xlab = "Temps",
     ylab = "Rendements")
abline(h= -VaR_CFisher, col="green")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_CFisher]), pch=20, col="green")
legend("top", legend = "VaR CFisher", lty = 1, col = "green")


#Représentation des VaR : HS, VC et Cornish Fisher et de leurs dépassements. 
X11()
plot(as.zoo(rend_FTSEMIB), 
     main = "Rendements et VaR non-conditionnelles",
     xlab = "Temps",
     ylab = "Rendements")
abline(h =c(-VaR_HS, -VaR_VCov, -VaR_CFisher), 
       col = c("red", "blue", "green"))
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_VCov]), pch=20, col="blue")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_HS]), pch=20, col="red")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_CFisher]), pch=20, col="green")
legend("top", 
       legend = c("VaR HS", "VaR Vcov", "VaR CFisher"),
       lty = 1, 
       col = c("red", "blue", "green"))

# Graphique
X11()
plot(as.zoo(rend_FTSEMIB),main="Les différentes mesures de la VaR",ylim=c(-0.12,0.12))
abline(h=c(-VaR_HS,-VaR_VCov,-VaR_CFisher), col=c("red","blue","green"))
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_VCov]), pch=20, col="blue")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_HS]), pch=20, col="red")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_CFisher]), pch=20, col="green")
legend("top",legend = c("VaR_HS","VaR_VCov","VaR_CFisher"), col=c("red","blue","green"),pch = "_")
#### Test sur la correlation des rendements


#Autocorrélogrammes des séries
X11()
plot(acf(rend_FTSEMIB,lag=100)[1:100,],lwd=3,main="Corrélation des rendements de l'indice FTSE MIB")
X11()
plot(acf(rend_FTSEMIB^2,lag=100)[1:100,],lwd=3,main="Corrélation des rendements au carré de l'indice FTSE MIB")

#Présence de corrélation dans les rendements et rendements au carré, ainsi modèle GARCH avec un processus ARMA pour la moyenne conditionnelle à envisager. 


#Tests de Ljung Box 
Box.Ljung.Test(rend_FTSEMIB,20)
Box.Ljung.Test(rend_FTSEMIB^2,40)

Box.test(rend_FTSEMIB,lag=1,type = "Ljung-Box")  
Box.test(rend_FTSEMIB,lag=5,type = "Ljung-Box")  
Box.test(rend_FTSEMIB,lag=20,type = "Ljung-Box") 

#On constate une autocorrélation significative au niveau des rendements du FTSE MIB à 1,5 et 20 retards

Box.test(rend_FTSEMIB^2,lag=1,type = "Ljung-Box") 
Box.test(rend_FTSEMIB^2,lag=5,type = "Ljung-Box")  
Box.test(rend_FTSEMIB^2,lag=20,type = "Ljung-Box") 

#On constate également une autocorrélation significative des rendements du FTSE MIB à 1,5 et 20 retards

# Modélisation ARMA de REND FTSE MIB

X11()
plot(acf(rend_FTSEMIB,lag=100)[1:100,],lwd=3,main="Corrélation totale de la série REND FTSEMIB")
X11()
pacf(rend_FTSEMIB, main="Corrélation partielle REND FTSEMIB",lag.max=100,na.action = na.pass)

# Meilleur modèle serait donc un processus ARMA
#
#### Sélection automatique du "meilleur" modèle ARMA

choix_mod<-auto.arima(rend_FTSEMIB,max.p=10,max.q = 10,d=0,D=0,max.P=2,max.Q=2,ic = "aic",stepwise=TRUE,stationary = TRUE,seasonal = TRUE,approximation = TRUE,allowmean = TRUE)
choix_mod


### Estimation du modèle

X11()
estimate(rend_FTSEMIB,p=0,d=0,q=1)

#Un modèle ARMA avec p =0, d=0, q=1 semble meilleur car minimisant les critères d'information AIC...

X11()
estimate(rend_FTSEMIB,p=1,d=0,q=1)

# Test de surdimensionnement
X11()
estimate(rend_FTSEMIB,p=1,d=0,q=1)

#
# Analyse des résidus
X11()
residu<-estimate(rend_FTSEMIB,p=1,d=0,q=1)$residuals
X11()
plot(as.zoo(residu))
#
# Test du portemanteau
Box.test(residu,lag=1,type = "Ljung-Box")
Box.test(residu,lag=5,type = "Ljung-Box")
Box.test(residu,lag=20,type = "Ljung-Box") # L'autocorrélation est significative à 20 retards

#on affiche les 3ème et 4ème moments
kurtosis(residu,na.rm = TRUE)
skewness(residu,na.rm = TRUE)
#
### Conclusion: on retient le ARMA(1,1) mais il peut rester de la correlation sur les residus



# Test Effet ARCH
# Test LM d'Engle (1982)  

arma11<-arima(rend_FTSEMIB,order = c(1,0,1),include.mean = TRUE)
X11()
arch.test(arma11,output = TRUE)

#Estimations de l'équation de variance conditionnelle

### Estimation d'un ARMA(1,1)-GARCH(1,1) pour rend_FTSEMIB avec z suivant N(0,1)
#
## Pour spécifier le modèle: ugarchspec
#
arma11garch11<-ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1),variance.targeting=FALSE),mean.model = list(armaOrder=c(1,1),include.mean=TRUE),distribution.model = "norm")
arma11garch11
#
## Pour éstimer le modèle spécifié: ugarchfit
#
arma11garch11_fit<-ugarchfit(arma11garch11,rend_FTSEMIB)
arma11garch11_fit
# Le biais de signe calcule le test de biais de signe d'Engle et Ng (1993) et est également affiché dans le résumé.
# Cela teste la présence d'effets de levier dans les résidus standardisés (pour capturer une éventuelle # erreur de spécification du modèle GARCH),
# en régressant les carrés des résidus standardisés sur les chocs négatifs et positifs décalés
#
persistence(arma11garch11_fit)
halflife(arma11garch11_fit)
# demi-vie = log(0.5).log(alpha1+beta1)
#
sd(as.numeric(rend_FTSEMIB)) # vol non conditionnelle calculée  de la série
sqrt(uncvariance(arma11garch11_fit)) # vol non conditionnelle estimée de la série

#volatilité non cond de la série et volatilité estimée similaires car ciblage de variance 

#
X11()
plot(arma11garch11_fit,which="all")
X11()
plot(arma11garch11_fit,which=8)
# Les résidus ne suivent pas exactement une loi de Student mais ils sen rapprochent plus que pour une loi normale

# graphe de la volatilité conditionnelle estimée

X11()
plot(as.zoo(sigma(arma11garch11_fit)),
     main="volatilité du FTSE MIB", 
     xlab = "Temps", 
     ylab = "Volatilité")




#Il est possible d’effectuer des améliorations, en effet il reste de la kurtosis passé , ainsi un GJR-GARCH permettrait de modéliser des effet asymétriques et un ciblage de variance pour tenir compte de la non significativité de la constante.

# GARCH Asymétrique: ARMA(1,1)-GJR-GARCH(1,1) et ciblage de variance pour rend_FTSE MIB avec z suivant t(v)
arma11gjrgarch11<-ugarchspec(variance.model = list(model="gjrGARCH",garchOrder=c(1,1),variance.targeting=TRUE),mean.model = list(armaOrder=c(1,1),include.mean=TRUE),distribution.model = "std")
arma11gjrgarch11_fit<-ugarchfit(arma11gjrgarch11,rend_FTSEMIB)
arma11gjrgarch11_fit
#
persistence(arma11gjrgarch11_fit)
halflife(arma11gjrgarch11_fit)
#
X11()
plot(arma11gjrgarch11_fit,which="all")

#
sd(as.numeric(rend_FTSEMIB)) # vol non conditionnelle calculée  de la série
sqrt(uncvariance(arma11gjrgarch11_fit)) # vol non conditionnelle estimée de la série
#
# Graphique de comparaison des vol estimées
X11()
plot(as.zoo(cbind(sigma(arma11garch11_fit),sigma(arma11gjrgarch11_fit))),screens=1,col = c("blue","red"),main="Volatilité conditionnelle estimée")
legend("topleft",c("GARCH(1,1)","GJR-GARCH(1,1)"),col = c("blue","red"),lty = c(1,1),lwd = c(2,2))


#############################  VAR CONDITIONNELLES


# Calculde la VaR Risk-metrics
#Calcul de la VaR
VaR_Risk <- (-as.numeric(quantile(arma11gjrgarch11_fit, probs = alpha)))
VaR_Risk

#Représentation Grapique de la VaR et de ses dépassements
X11()
display = {plot(as.zoo(cbind(rend_FTSEMIB,-VaR_Risk)),
                screens=1,
                col = c("black","orange"),
                main="VaR ARMA(1,1)-RiskMetrics", 
                xlab = "Temps", 
                ylab = "Rendements")
  
  points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_Risk]), pch=20, col="blue")
  
  legend("bottom",c("Rendements","VaR RiskMetrics", "Dépassements de la VaR"),
         col = c("black","orange", 'blue'),lty = c(1,1,NA),lwd = c(2,2,NA),pch = c(NA,NA,20))}

### Calcul de la VaR GARCH
alpha<- 0.01
X11()
mu<-fitted(arma11gjrgarch11_fit)
X11()
plot(as.zoo(cbind(rend_FTSEMIB,mu)),screens=1,col = c("black","red"),main="Rend observé vs Rend estimé")
legend("topleft",c("Observé","Estimé"),col = c("black","red"),lty = c(1,1),lwd = c(2,2))
#
## Calcul "manuel" de la VaR GARCH 
(v <- arma11gjrgarch11_fit@fit$coef["shape"]) 
VaR_Garch_calcul <- as.numeric(mu + sigma(arma11gjrgarch11_fit)*sqrt((v-2)/v) * qt(alpha, df = v)) 
X11()
plot(as.zoo(VaR_Garch_calcul),main="VaR GARCH estimée")
#
# Calcul automatique
VaR_Garch <- (-as.numeric(quantile(arma11gjrgarch11_fit, probs = alpha)))
VaR_Garch
# Comparaison
all.equal(-VaR_Garch_calcul, VaR_Garch)

#
# Graphique
X11()
plot(as.zoo(cbind(rend_FTSEMIB,-VaR_Garch)),screens=1,col = c("black","pink"),main="VaR Garch")
points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_Garch]), pch=20, col="blue")
legend("bottom",c("Rend","VaR99%"),col = c("black","pink"),lty = c(1,1),lwd = c(2,2))
#
# Graphique de comparaison des différentes VaR
q <- quantile(as.numeric(rend_FTSEMIB), probs=alpha)
(VaR_HS<-(-q))
#
(VaR_VCov<-(-(mean(rend_FTSEMIB,na.rm=TRUE)+sd(rend_FTSEMIB,na.rm=TRUE)*qnorm(0.01))))
#
library(DistributionUtils)
K<-kurtosis(rend_FTSEMIB, na.rm = TRUE)
S<-skewness(rend_FTSEMIB, na.rm = TRUE)
(k<-as.numeric(qnorm(0.01)+(qnorm(0.01)^2-1)*(S/6)+(qnorm(0.01)^3-3*qnorm(0.01))*(K/24)+(2*qnorm(0.01)^3-5*qnorm(0.01))*(S^2/36)))
(VaR_CFisher<-as.numeric((-sd(rend_FTSEMIB)*k)))
#
(fitdist(distribution = 'std' , x = rend_FTSEMIB)$pars)
(qt(0.01,df=fitdist(distribution = 'std' , x = rend_FTSEMIB)$pars["shape"]))
(VaR_student<-(-sd(rend_FTSEMIB)*qt(0.01,df=fitdist(distribution = 'std' , x = rend_FTSEMIB)$pars["shape"])*sqrt((qt(0.01,df=fitdist(distribution = 'std' , x = rend_FTSEMIB)$pars["shape"])-2)/qt(0.01,df=fitdist(distribution = 'std' , x = rend_FTSEMIB)$pars["shape"]))))
#

#VaR TVE
install.packages("extRemes")
library(extRemes)
install.packages("moments")

phi <- 80 #taille des blocs 
num_blocs <- ceiling(length(rend_FTSEMIB) / phi) #Nombre de blocs

valeurs_min_blocs <- c()

for (i in 1:num_blocs) {
  debut <- (i - 1) * phi + 1
  fin <- min(i * phi, length(rend_FTSEMIB))
  bloc <- rend_FTSEMIB[debut:fin]
  valeur_min <- min(bloc)
  valeurs_min_blocs <- c(valeurs_min_blocs, valeur_min)
}

valeurs_min_blocs

paramètres.GEV <- fevd(valeurs_min_blocs, type="GEV") #Estimation des paramètres
mu <- paramètres.GEV$results$par[1] #Paramètre Mu
sigma <- paramètres.GEV$results$par[2] #Paramètre sigma
epsilon <- paramètres.GEV$results$par[3] #Paramètre epsilon

qtve <- 1-phi*(1-0.99) #quantile TVE
(VaR_TVE <- as.numeric(mu-(sigma*epsilon**(-1))*(1-(-log(qtve))**(-epsilon))))
VaR_TVE

X11()
{plot(as.zoo(rend_FTSEMIB),
      main = "VaR TVE et dépassements",
      xlab = "Temps",
      ylab = "Rendements")
  abline(h=c(VaR_TVE, -VaR_VCov), col=c('red', 'blue'))
  points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_VCov]), pch=20, col="blue")
  points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < VaR_TVE]), pch=20, col="red")
  legend("bottom", legend = c("VaR gaussienne", "VaR TVE Phi = 99"), lty = 1, col=c("blue", "red"))}

# La VaR GARCh TVE est estimée à partir du modèle ARMA(1,1)-GJR-GARCH(1,1)


arma11gjrgarch11<-ugarchspec(variance.model = list(model="gjrGARCH",garchOrder=c(1,1),variance.targeting=TRUE),mean.model = list(armaOrder=c(1,1),include.mean=TRUE),distribution.model = "std")
# 
arma11gjrgarch11_fit<-ugarchfit(arma11gjrgarch11,rend_FTSEMIB)

residus <- arma11gjrgarch11_fit@fit$z #récupèration des résidus du modèle GJRGARCH


#On estime la TVE avec un nombre de 80 valeurs par blocs
phi=80 #taille des blocs
num_blocs <- ceiling(length(residus) / phi) #Calcul du nombre de blocs en fonction de phi et des données
valeurs_min_blocs <- c() #On initie un vecteur vide qui contiendra nos valeurs extrêmes de chaque bloc

for (i in 1:num_blocs) {  
  debut <- (i - 1) * phi + 1
  fin <- min(i * phi, length(residus))
  bloc <- residus[debut:fin]
  valeur_min <- min(bloc)
  valeurs_min_blocs <- c(valeurs_min_blocs, valeur_min) 
}

paramètres.GEV <- fevd(valeurs_min_blocs, type="GEV") #Estimation des paramètres de la GEV sur l'échantillon de valeurs maximales
mu <- paramètres.GEV$results$par[1] #Paramètre Mu
sigma <- paramètres.GEV$results$par[2] #Paramètre sigma
epsilon <- paramètres.GEV$results$par[3] #Paramètre epsilon

xq <- mu - (sigma*epsilon**(-1))*(1-(-log(0.99))**(-epsilon))


#Calcul de la VaR TVE
VaR_GJRGARCH_TVE_80 <- (-arma11gjrgarch11_fit@fit$fitted.values -arma11gjrgarch11_fit@fit$sigma*xq)


#Représentation graphique de la VaR et de ses dépassements
X11()
{plot(as.zoo(cbind(rend_FTSEMIB,-VaR_GJRGARCH_TVE_80)),
      screens=1,
      col = c("black","violet"),
      main="VaR ARMA(1,1)-GJRGARCH(1,1) TVE-80", 
      xlab = "Temps", 
      ylab = "Rendements")
  
  points(as.zoo(rend_FTSEMIB[as.numeric(rend_FTSEMIB) < -VaR_GJRGARCH_TVE_80]), pch=20, col="blue")
  
  legend("bottom",c("Rendements","VaR GJRGARCH(1,1) TVE-80", "Dépassements de la VaR"),
         col = c("black","violet", 'blue'),lty = c(1,1,NA),lwd = c(2,2,NA),pch = c(NA,NA,20))}


#Backtesting des VaR 


# fonction HIT avec probabilités de dépassements sachant un dépassement la veille ou non. 
fonc_N_ij <- function(HIT,i,j){
  compteur <- 0
  for(k in 2:length(HIT)){
    if(HIT[k]==j & HIT[k-1]==i){ 
      compteur <- compteur+1 
    }
  }
  return(compteur)
}



# test de Kupiec et de Christoffersen / teste l'indépendance et le taux de couverture non conditionnelle. 

# On obtiens la statistique LM du taux de couverture conditionnelle, la stat LM du taux de couverture non conditionnelle et la statistique LM de l'indépendance. 
Backtest_function <- function(HIT, sample, alpha=0.01){
  
  N <- sum(HIT) #Nombre total de dépassement
  Total <- length(sample) #Nombre d'observations 
  failure.rate <- N/Total #Ratio de dépassement
  
  n_00 <- fonc_N_ij(as.numeric(HIT),0,0) #Le nombre de fois que HIT_t=0 sachant que HIT_t-1 = 0
  n_01 <- fonc_N_ij(as.numeric(HIT),0,1) #Le nombre de fois que HIT_t=1 sachant que HIT_t-1 = 0
  n_10 <- fonc_N_ij(as.numeric(HIT),1,0) #Le nombre de fois que HIT_t=0 sachant que HIT_t-1 = 1
  n_11 <- fonc_N_ij(as.numeric(HIT),1,1) #Le nombre de fois que HIT_t=1 sachant que HIT_t-1 = 1
  
  pi_01 <- n_01/(Total-N) #Nombre de dépassement sachant que la veille n'avait pas dépassée sur le nombre de jours sans dépassements
  pi_11 <- n_11/N #Nombre de dépassements sachant que la veille était un dépassement sur le nombre de dépassements
  
  LR_cc <- -2*log(((1-alpha)**(Total-N))*alpha**N)+
    2*log(((1-pi_01)**n_00)*(pi_01**n_01)*((1-pi_11)**n_10)*pi_11**n_11) #Statistique LM de couverture conditionnelle : CHi-deux(2)
  
  LR_uc <- 2*log((((1-failure.rate)/(1-alpha))**(Total-N))*(failure.rate/alpha)**N)#statistique LM de couverture non conditionnelle : Chi-deux(1)
  
  LR_ind <- -2*log(((1-(N/Total))**(Total-N))*(N/Total)**N)+
    2*log(((1-pi_01)**n_00)*(pi_01**n_01)*((1-pi_11)**n_10)*pi_11**n_11) #statistique LM déindépendance : Chi-deux(1)
  
  return(c(LR_cc, LR_uc, LR_ind))
  
}





#Création des sample pour estimer et tester les VaR
taille.echantillon <- length(rend_FTSEMIB)-1000
sample_test <- rend_FTSEMIB[1:taille.echantillon] #sample sur lequel on estime la VaR
sample_backtest <- rend_FTSEMIB[(taille.echantillon+1):length(rend_FTSEMIB)] #1000 jours de test de la VaR



#Estimation de la VaR historique sur l'échantillon de test
{alpha <- 0.01
  q <- quantile(as.numeric(sample_test), probs=alpha)
  (VaR_HS_test<-(-q))}

#Création de la fonction HIT qui prend la valeur 1 en cas de dépassements de la VaR
HIT_HS <- ifelse(sample_backtest < -VaR_HS_test, 1, 0)

#Résultats du backtesting
(res_backtest_HS <- Backtest_function(HIT_HS, sample_backtest))



#Estimation de la VaR variance-covariance sur l'échantillon de test
(VaR_VCov_test<-(-(mean(sample_test)+sd(sample_test)*qnorm(0.01))))

#Création de la fonction HIT qui prend la valeur 1 en cas de dépassements de la VaR
HIT_Vcov <- ifelse(sample_backtest < -VaR_VCov_test, 1, 0)

#Résultats du backtesting
(res_backtest_VCov <- Backtest_function(HIT_Vcov, sample_backtest))




#Estimation de la VaR Cornish-Fisher sur l'échantillon de test
{K<-kurtosis(sample_test)
  S<-skewness(sample_test)
  
  (k<-as.numeric(qnorm(0.01)+(qnorm(0.01)^2-1)*(S/6)+(qnorm(0.01)^3-3*qnorm(0.01))*(K/24)+(2*qnorm(0.01)^3-5*qnorm(0.01))*(S^2/36)))
  
  (VaR_CFisher_test<-as.numeric((-sd(sample_test)*k)))}

#Création de la fonction HIT qui prend la valeur 1 en cas de dépassements de la VaR
HIT_Cfisher <- ifelse(sample_backtest < -VaR_CFisher_test, 1, 0)


#Résultats du backtesting
(res_backtest_Cfisher <- Backtest_function(HIT_Cfisher, sample_backtest))




#Test sur la VaR ARMA(1,1)-GJRGARCH(1,1) 

#Spécification du modèle 
spec <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), variance.targeting=TRUE),
                   mean.model = list(armaOrder = c(1, 1), include.mean=TRUE),
                   distribution.model = "std")

#J'estime le modèle sur toutes les observations sauf les 1000 dernières
fit <- ugarchfit(spec, data = rend_FTSEMIB, out.sample = 1000)
coefficients <- coef(fit)
shape <- coefficients["shape"]

#Je forecast en fenêtre glissante à 1 pas de prévision pour le reste des valeurs en prenant les paramètres estimés sur les 3000 premières valeurs. 
forecast <- ugarchforecast(fit, data = rend_FTSEMIB, n.roll = 999, n.ahead = 1)

#J'enregistre les valeurs de mu et sigma estimées
sigma_gjrgarch <- forecast@forecast$sigmaFor
mu_gjrgarch <- forecast@forecast$seriesFor



#Calcul de la VaR GARCH suivant une loi de student
VaR_GJRGarch11_test <- as.numeric(mu_gjrgarch)+as.numeric(sigma_gjrgarch)*sqrt((shape-2)/shape)* qt(alpha, df = shape)


#Création de la variable HIT
HIT_GJRgarch <- ifelse(sample_backtest < VaR_GJRGarch11_test, 1, 0)

#Test de backtesting
(res_backtest_GJRgarch <- Backtest_function(HIT_GJRgarch, sample_backtest))

#Test sur la VaR TVE

# Spécification du modèle GJR-GARCH(1,1) avec une moyenne ARMA(1,1)
spec_tve <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), variance.targeting=TRUE),
                       mean.model = list(armaOrder = c(1, 1), include.mean=TRUE),
                       distribution.model = "std")

# Estimation du modèle sur toutes les observations sauf les 1000 dernières
fit_tve <- ugarchfit(spec_tve, data = rend_FTSEMIB, out.sample = 1000)
coefficients <- coef(fit_tve)
shape <- coefficients["shape"]

# Forecast en fenêtre glissante à 1 pas de prévision pour le reste des valeurs en prenant les paramètres estimés sur les 3000 premières valeurs
forecast_tve <- ugarchforecast(fit_tve, data = rend_FTSEMIB, n.roll = 999, n.ahead = 1)

# Enregistrement des valeurs de mu et sigma estimées
sigma_tve <- forecast_tve@forecast$sigmaFor
mu_tve <- forecast_tve@forecast$seriesFor

# Calcul de la VaR TVE suivant une loi de Student
VaR_tve_test <- as.numeric(mu_tve) + as.numeric(sigma_tve) * sqrt((shape - 2) / shape) * qt(alpha, df = shape)

# Création de la variable HIT
HIT_tve <- ifelse(sample_backtest < VaR_tve_test, 1, 0)

# Test de backtesting
(res_backtest_tve <- Backtest_function(HIT_tve, sample_backtest))
