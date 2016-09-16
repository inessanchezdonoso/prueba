install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")
install.packages("bbmle")
library(xlsx)
library(Deducer)
library(ggplot2)
library(boot)
library(car)
library(MASS)
library(survival)
library(RVAideMemoire)
library(doBy)
library(AICcmodavg)
library(ggplot2) 
library(betareg)
library(glmmADMB)
library(bbmle) 
library(lme4) 
library(lattice)
citation("glmmADMB")


setwd("C:/Users/Usuario/Dropbox/Research/5 Varios/Carcass degradation_Jen")


basefinal<-read.table("basefinal.txt",header=T, dec=",")
d=subset(basefinal, basefinal$Time!=0)  #sin casos en que el tiempo es 0.

str(d)
  d$Treecoverage<-factor(d$Treecoverage)
  d$PhaseInes<-factor(d$PhaseInes)
  d$InflectionPoint1<-factor(d$InflectionPoint1)
  d$InflectionPoint2<-factor(d$InflectionPoint2)

  
  
###SCI over time, según Class
#############################
  
#Con Time ajuste malo:
  
ggplot(data=d, aes(x=Time, y=SCIprop, group=Class, colour=Class)) + geom_point(size = 3, shape = 1) + theme_bw() 
  
  m1<-glmmadmb(SCIprop~Class+Time+Treecoverage+Locality+Sps+Start+Age+(1|ID), family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m1), residuals(m1), main="m1")
    abline(0,0)
    qqnorm(residuals(m1), main="qqnorm m1")
    qqline(residuals(m1)) #Ajuste malo

  m2<-glmmadmb(SCIprop~Class+Time+Treecoverage+Locality+(1|ID), family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m2), residuals(m2), main="m2")
    abline(0,0)
    qqnorm(residuals(m2), main="qqnorm m2")
    qqline(residuals(m2)) #Ajuste malo
    summary(m2)

  m3<-glmmadmb(SCIprop~Class*Time+Treecoverage+Locality+(1|ID), family="beta", data=d)  
      par(mfrow=c(1, 2))
      plot(fitted(m3), residuals(m3), main="m3")
      abline(0,0)
      qqnorm(residuals(m3), main="qqnorm m3")
      qqline(residuals(m3)) #Ajuste malo
      summary(m3)

#Además, con Time yo veo que los datos no siguen una distribución lineal, que es indispensable para que pueda analizarlos
#con modelos lineales. Tengo que transformar alguna variable para que mis datos sigan una distrib lineal y entonces 
#poderlos analizar con GLM o GLMM. Si no, debería meterme en modelos no lineales.


#Transformo Time
d$eTime<-exp(-d$Time)

ggplot(data=d, aes(x=eTime, y=SCIprop, group=Class, colour=Class)) + geom_point(size = 3, shape = 1) + theme_bw() 

  m1<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+Sps+Start+Age+(1|ID), family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m3), residuals(m3), main="m3")
    abline(0,0)
    qqnorm(residuals(m3), main="qqnorm m3")
    qqline(residuals(m3))
    #Ajuste ok
    summary(m3) #no sig: Age, Start, Locality aquí sale no sig, pero en el siguiente (sin Age ni Start) sale sig!
                #Sig: Class III, eTime, Treecoverage
    

  m2<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+Start+Age+(1|ID), family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m2), residuals(m2), main="m2")
    abline(0,0)
    qqnorm(residuals(m2), main="qqnorm m2")
    qqline(residuals(m2))
    #Ajuste ok
    summary(m2) 
    
    
    #locality me sale sig en cuanto saco Sps. Éso significa que locality tiene algo que ver con sps. 
    #Efectivamente, hay un sesgo claro entre localidades: las especies. En SNS sólo hay C.elaphus.
    
    #Para comparar localidades debería comparar sólo Cervus elaphus, ya que es la única especie que tiene SNS.
    #No es la localidad lo interesante, si no la especie.
    
    Cervus=subset(d, d$Sps=="Celaphus")
    str(Cervus)
    
    m<-glmmadmb(SCIprop~eTime+Treecoverage+Locality+(1|ID), family="beta", data=Cervus) #he tenido que eliminar Class, claro
    par(mfrow=c(1, 2))
    plot(fitted(m), residuals(m), main="m")
    abline(0,0)
    qqnorm(residuals(m), main="qqnorm m")
    qqline(residuals(m))
    summary(m)
    
 
    
  m3<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+Age+(1|ID), family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m3), residuals(m3), main="m3")
    abline(0,0)
    qqnorm(residuals(m3), main="qqnorm m3")
    qqline(residuals(m3))
    #Ajuste ok
    summary(m3)
    
    
  m4<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+(1|ID), family="beta", data=d)  #ESTE EN ARTÍCULO
                                                    #Indicando con o sin zeroInflation me sale lo mismo
    par(mfrow=c(1, 2))
    plot(fitted(m4), residuals(m4), col=as.numeric(d$Class), main="m4")
    abline(0,0)
    qqnorm(residuals(m4), main="qqnorm m4")
    qqline(residuals(m4))
    #Ajuste ok
    summary(m4) #todo sig. Classes: ClassI=ClassII>ClassIII, eTime, Locality
    d$Class=relevel(d$Class,ref="III")
    #Por default usa el link logit. 
    

ggplot(data=d, aes(x=Time, y=SCIprop, group=Treecoverage, colour=Treecoverage)) + geom_point(size = 3, shape = 1) + theme_bw()     

ggplot(data=d, aes(x=Time, y=SCIprop, group=Locality, colour=Locality)) + geom_point(size = 3, shape = 1) + theme_bw() 


Cervus=subset(d, d$Sps=="Celaphus")
RBD=subset(d, d$Locality=="RBD")

m<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Sps+(1|ID), family="beta", data=RBD)
par(mfrow=c(1, 2))
plot(fitted(m), residuals(m), main="m")
abline(0,0)
qqnorm(residuals(m), main="qqnorm m")
qqline(residuals(m))
summary(m)
Anova(m, test.statistic = "Chisq")
#Cuando trabajamos con sólo Doñana, se vé que sps no tiene efecto. La variabilidad de especie está recogida de alguna
#manera en size class.

SNS=subset(d, d$Locality=="SNS")

m<-glmmadmb(SCIprop~eTime+Treecoverage+(1|ID), family="beta", data=SNS)
par(mfrow=c(1, 2))
plot(fitted(m), residuals(m), main="m")
abline(0,0)
qqnorm(residuals(m), main="qqnorm m")
qqline(residuals(m))
#Ajuste ok
summary(m) 

#Los datos de SNS son los que están dando la significación a Tree coverage


##Pruebas:
  #Prueba Ccon interacción clase y tiempo:
  m5<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+(1|ID), 
             zeroInflation=TRUE, family="beta", data=d)
      par(mfrow=c(1, 2))
      plot(fitted(m5), residuals(m5), main="m5")
      abline(0,0)
      qqnorm(residuals(m5), main="qqnorm m5")
      qqline(residuals(m5))    #Ajuste ok
      summary(m5)
      #La interacción Class y eTime no es significativa.
      d$Class=relevel(d$Class,ref="II")

  #Otras pruebas:
  m5<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+Sps+(1|ID), 
             zeroInflation=TRUE, family="beta", data=d)
    par(mfrow=c(1, 2))
    plot(fitted(m5), residuals(m5), main="m5")
    abline(0,0)
    qqnorm(residuals(m5), main="qqnorm m5")
    qqline(residuals(m5))     #Ajuste ok
    summary(m5)

  m6<-glmmadmb(SCIprop~Class+eTime+Treecoverage+(1|ID), 
             zeroInflation=TRUE, family="beta", data=d) #No me lo hace, no se por qué
    par(mfrow=c(1, 2))
    plot(fitted(m6), residuals(m6), main="m6")
    abline(0,0)
    qqnorm(residuals(m6), main="qqnorm m6")
    qqline(residuals(m6)) 
    summary(m6)

    

    
#DIFERENCIAS DE SCI EN INFLECTION POINTS ENTRE CLASES
#####################################################

    
#InflectionPoint1:
###################
#Hago una selección sólo de las observaciones que son punto de inflexión entre Phase1-2
#OJO!! uso la base de datos en la que si hay Time=0, porque hay algunos casos en que 
    #el InflectionPoint1 es en T=0.
    
  InfPoint1=subset(basefinal, basefinal$InflectionPoint1==1)
  str(InfPoint1)

  InfPoint1RBD=subset(InfPoint1, InfPoint1$Locality=="RBD")
  InfPoint1SNS=subset(InfPoint1, InfPoint1$Locality=="SNS")


  #VR SCIprop:
    m=betareg(SCIprop~Class+Locality+Treecoverage+Sps+Age+Start, family="beta", data=InfPoint1) #ERROR varios
    #quizá si le reduzco VE
    
    m=betareg(SCIprop~Class+Locality+Treecoverage,  data=InfPoint1)  #ERROR. Pero sale nada sig.
    summary(m)

    m4=betareg(SCIprop~Class+Locality, data=InfPoint1) #ESTE Y...
        par(mfrow=c(1, 2))
        plot(fitted(m4), residuals(m4), main="m4")
        abline(0,0)
        qqnorm(residuals(m4), main="qqnorm m4")
        qqline(residuals(m4)) #Ajuste muy raro.
      summary(m4) #Nada es sig. 


    m6=betareg(SCIprop~Class+Treecoverage, data=InfPoint1)  #... ESTE
        par(mfrow=c(1, 2))
        plot(fitted(m6), residuals(m6), main="mm64")
        abline(0,0)
        qqnorm(residuals(m6), main="qqnorm m6")
        qqline(residuals(m6)) #Ajuste muy raro.
      summary(m6) #Nada es sig.

    m5=betareg(SCIprop~Class, data=InfPoint1)
        par(mfrow=c(1, 2))
        plot(fitted(m5), residuals(m5), main="m5")
        abline(0,0)
        qqnorm(residuals(m5), main="qqnorm m5")
        qqline(residuals(m5)) #Ajuste muy raro.
      summary(m5) #No sig. 
    

    kruskal.test(SCIprop ~ Class, data = InfPoint1) #No sig
  #No diferencias entre clases ni con betareg (ni controlado por Locality ni 
  #controlado por treecoverage), ni con Kruskal-Wallis (prueba no paramétrica 
  #para varias muestras independientes).

  #Ploteo:
      par(mfrow=c(1, 3))
      plot(InfPoint1$Class, InfPoint1$SCIprop, ylab="SCIprop", xlab="Class",  main="both loc")
      plot(InfPoint1RBD$Class, InfPoint1RBD$SCIprop, ylab="SCIprop", xlab="Class", main="RBD")
      plot(InfPoint1SNS$Class, InfPoint1SNS$SCIprop, ylab="SCIprop", xlab="Class", main="SNS")
    #Poner el plot the both loc, porque no hay diferencias entre locs.

  
  #Creo que si lo hiciera sin los datos T=0, quizá el modelo se ajustaría mejor.
  #Ahora creo que mejor incluirle los casos en que T=0.
    InfPoint1.sinT0=subset(d, d$InflectionPoint1==1)
    
      m2=betareg(SCIprop~Class+Locality+Start+Treecoverage, data=InfPoint1.sinT0) #ERROR

      m3=betareg(SCIprop~Class+Locality+Treecoverage, data=InfPoint1.sinT0) #funciona
        par(mfrow=c(1, 2))
        plot(fitted(m3), residuals(m3), main="m3")
        abline(0,0)
        qqnorm(residuals(m3), main="qqnorm m3")
        qqline(residuals(m3)) #Ajuste muy raro.
        summary(m3) #no lo calcula. Errores con NA ????

      m4=betareg(SCIprop~Class+Locality, data=InfPoint1.sinT0) 
        par(mfrow=c(1, 2))
        plot(fitted(m4), residuals(m4), main="m4")
        abline(0,0)
        qqnorm(residuals(m4), main="qqnorm m4")
        qqline(residuals(m4)) #Ajuste muy raro.
        summary(m4) #Class II sig.

      m5=betareg(SCIprop~Class, data=InfPoint1.sinT0) 
        par(mfrow=c(1, 2))
        plot(fitted(m5), residuals(m5), main="m5")
        abline(0,0)
        qqnorm(residuals(m5), main="qqnorm m5") #este no lo calcula
        qqline(residuals(m5)) #Ajuste muy raro.
        summary(m5) #Class II sig.
    #No me gustan estos modelos. No son buenos.


  #VR Time:

  mt1=glm(Days~Class+Locality+Start+Treecoverage, data=InfPoint1, family="poisson")
    par(mfrow=c(1, 2))
    plot(fitted(mt1), residuals(mt1), main="mt1")
    abline(0,0)
    qqnorm(residuals(mt1), main="qqnorm mt1")
    qqline(residuals(mt1))  #Ajuste aceptable.
    summary(mt1) 
    #Cambio el orden de los niveles del factor Class, para testar diferencias entre Class II y Class III:
    InfPoint1$Class=relevel(InfPoint1$Class,ref="II")
    #Treecoverage sig.
    #ClassII>ClassI=ClassIII

  mt2=glm(Days~Class+Treecoverage, data=InfPoint1, family="poisson") #ESTE
    par(mfrow=c(1, 2))
    plot(fitted(mt2), residuals(mt2), main="mt2")
    abline(0,0)
    qqnorm(residuals(mt2), main="qqnorm mt2")
    qqline(residuals(mt2))  #Ajuste aceptable.
    summary(mt2)
    #Cambio el orden de los niveles del factor Class, para testar diferencias entre Class II y Class III:
    InfPoint1$Class=relevel(InfPoint1$Class,ref="I")
    #Treecoverage sig.
    #ClassII>ClassI=ClassIII

    ggplot(data=d, aes(x=Time, y=SCIprop, shape=InflectionPoint1, colour=Class)) + geom_point() + theme_bw() + geom_point(size = 3)

  #Ploteo:
    par(mfrow=c(1, 3))
    plot(InfPoint1$Class, InfPoint1$Time, ylab="Time", xlab="Class",  main="both loc")
    plot(InfPoint1RBD$Class, InfPoint1RBD$Time, ylab="Time", xlab="Class", main="RBD")
    plot(InfPoint1SNS$Class, InfPoint1SNS$Time, ylab="Time", xlab="Class", main="SNS")
    #Poner el plot de both loc porque no hay dif entre locs.


#InflectionPoint2:
#################    
#Hago una selección sólo de las observaciones que son punto de inflexión entre Phase2 y Phase3
    
  InfPoint2=subset(basefinal, basefinal$InflectionPoint2==1)
  str(InfPoint2)
  InfPoint2RBD=subset(InfPoint2, InfPoint2$Locality=="RBD")
  InfPoint2SNS=subset(InfPoint2, InfPoint2$Locality=="SNS")

  #VR SCIprop:
  
    m3=betareg(SCIprop~Class+Locality+Treecoverage, data=InfPoint2)
        par(mfrow=c(1, 2))
        plot(fitted(m3), residuals(m3), main="m3")
        abline(0,0)
        qqnorm(residuals(m3), main="qqnorm m3")
        qqline(residuals(m3)) #Buen ajuste!!
        summary(m3)
        #Locality tiene efecto. En SNS pasa a PhaseIII con un SCI más bajo. Treecoverage no tiene efecto.

    m4=betareg(SCIprop~Class+Locality, data=InfPoint2) #ESTE
        par(mfrow=c(1, 2))    
        plot(fitted(m4), residuals(m4), main="m4")
        abline(0,0)
        qqnorm(residuals(m4), main="qqnorm m4")
        qqline(residuals(m4)) #Buen ajuste!!
        summary(m4)
        InfPoint2$Class=relevel(InfPoint2$Class,ref="III")
        #Locality sig. En SNS pasa a PhaseIII con un SCI más bajo.
        #Class: ClassI~>ClassII>ClassIII

    m5=betareg(SCIprop~Class, data=InfPoint2RBD)  #Y ESTE
        par(mfrow=c(1, 2))
        plot(fitted(m5), residuals(m5), main="m5")
        abline(0,0)
        qqnorm(residuals(m5), main="qqnorm m5")
        qqline(residuals(m5)) #Buen ajuste!!
        summary(m5)
        InfPoint2RBD$Class=relevel(InfPoint2RBD$Class,ref="III")
        
        #m=glmmadmb(SCIprop~Class+Locality+Start+Treecoverage, data=InfPoint, family="beta")
        #summary(m)
        #glm con family beta no me deja. No sé si glmmadmb sería correcto. 
        #No puedo testar todas las VE que querría, porque hay más variables que casos. Así que no pruebo ni Age, ni Sps, porque 
        #ya sé de los otros modelos que no me afectan la SCIprop.
        #Start y Treecoverage no tiene efecto. Locality marginalm, Class III dif, Class II marg

  #Ploteo:
      par(mfrow=c(1, 3))
      plot(InfPoint2$Class, InfPoint2$SCIprop, ylab="SCIprop", xlab="Class",  main="both loc")
      plot(InfPoint2RBD$Class, InfPoint2RBD$SCIprop, ylab="SCIprop", xlab="Class", main="RBD")
      plot(InfPoint2SNS$Class, InfPoint2SNS$SCIprop, ylab="SCIprop", xlab="Class", main="SNS")
      #Cuando compara las clases, lo hace sólo con los datos de RBD. Por lo que debería poner el boxplot por clases sólo de RBD.

    ggplot(data=d, aes(x=Time, y=SCIprop, shape=InflectionPoint, colour=Class)) + geom_point() + 
      theme_bw() + geom_point(size = 3)


  #VR Time:
    plot(d2InfPoint$Class, d2InfPoint$Days, ylab="Days", xlab="Class")

    mt1=glm(Days~Class+Locality+Start+Treecoverage, data=InfPoint2, family="poisson")
        summary(mt1)
        #Cambio el orden de los niveles del factor Class, para testar diferencias entre Class II y Class III:
        InfPoint$Class=relevel(InfPoint$Class,ref="II")

        #ClassII>ClassI>ClassIII
          ggplot(data=d, aes(x=Time, y=SCIprop, shape=InflectionPoint, colour=Class)) + geom_point() + theme_bw() + geom_point(size = 3)
        #LocalitySNS-
          ggplot(data=d, aes(x=Time, y=SCIprop, shape=InflectionPoint, colour=Locality)) + geom_point() + theme_bw() + geom_point(size = 3)
        #StartWST+
          ggplot(data=d, aes(x=Time, y=SCIprop, shape=InflectionPoint, colour=Start)) + geom_point() + theme_bw() + geom_point(size = 3)

    #creo que no tiene sentido meter Start, porque sólo tengo 1 caso de WST
    mt2=glm(Days~Class+Locality, data=InfPoint2, family="poisson")
        summary(mt2)
        InfPoint2$Class=relevel(InfPoint2$Class,ref="I")
        #LocalitySNS+
          ggplot(data=d, aes(x=Days, y=SCIprop, shape=InflectionPoint, colour=Locality)) + geom_point() + theme_bw() + geom_point(size = 3)
        #ClassII>ClassI>ClassIII
          ggplot(data=d, aes(x=Days, y=SCIprop, shape=InflectionPoint, colour=Class)) + geom_point() +theme_bw() + geom_point(size = 3)

        ggplot(data=d, aes(x=Days, y=SCIprop, shape=InflectionPoint, colour=Class)) + geom_point() + theme_bw() + geom_point(size = 3)

        #Cómo puede cambiar el signo de la significación de una VE, cuando tiene o no otra VE no sig en el modelo? MANUELA!!!

  #Ploteo:
    par(mfrow=c(1, 3))
    plot(InfPoint2$Class, InfPoint2$Days, ylab="Days", xlab="Class",  main="both loc")
    plot(InfPoint2RBD$Class, InfPoint2RBD$Days, ylab="Days", xlab="Class",  main="RBD")
    plot(InfPoint2SNS$Class, InfPoint2SNS$Days, ylab="Days", xlab="Class",  main="SNS")

  #Cuando compara las clases, lo hace sólo con los datos de RBD. Por lo que debería poner el boxplot por clases sólo de RBD.
  #Si comparo sólo las clases en RBD, porque es donde tengo de las 3:
    ggplot(data=d, aes(x=Days, y=SCIprop, shape=InflectionPoint, colour=Locality)) + geom_point() + theme_bw() + geom_point(size = 3)
    ggplot(data=d, aes(x=Days, y=SCIprop, shape=InflectionPoint, colour=Class)) + geom_point() + theme_bw() + geom_point(size = 3)

    mt3=glm(Days~Class, data=InfPoint2RBD, family="poisson")
        summary(mt3)
        #ClassII>ClassI>ClassIII
        #No le veo sentido, creo que es por la N tan pequeña.
        #No creo que pueda hacer este test, con esta N tan pequeña.

      #OJO!!! De estos modelos también debería testar su ajuste a las asunciones, no?




#Ahora trabajo con PhaseInes para cada uno de los puntos
  dPhaseInes=subset(d, d$PhaseInes!="NA")

str(d)
m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+Sps+Locality+Start+Age+Treecoverage+(1|ID), family="poisson", data=d)
#negative values not allowed for the 'Poisson' family"

m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+Sps+Locality+Start+Age+Treecoverage+(1|ID), 
               family="poisson", link="probit", data=d) 
#negative values not allowed for the 'Poisson' family"

m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+Sps+Locality+Start+Age+Treecoverage+(1|ID), 
               family="poisson", link="logit", data=d) 
#negative values not allowed for the 'Poisson' family"

m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+Sps+Locality+Start+Age+Treecoverage+(1|ID), 
               family="binomial", link="logit", data=d) 
#pwrssUpdate did not converge in (maxit) iterations

m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+Sps+Locality+Start+Age+Treecoverage+(1|ID), 
               family="binomial", link="probit", data=d) 
#pwrssUpdate did not converge in (maxit) iterations


m1.phase=glmer(PhaseInes~Class+eTime+SCIprop+(1|ID), 
               family="poisson", data=d) 

m1.phase=nlmer(PhaseInes~Class+eTime+SCIprop+ID, 
               family="poisson", data=d) 

m1.phase=glmer(PhaseInes~Class+(1|ID), 
               family="", data=d) 

kruskal.test(PhaseInes~Class, data=d)


Probit/Logit Regression 


anova(PhaseInes)






#Anterior____________________________
d<-read.table("soon1.txt",header=T)
d$Treecoverage<-factor(d$Treecoverage)
d$Phase<-factor(d$Phase)
str(d)
d3=subset(d, d$Time!=0)
str(d3)

par(mfrow=c(1,1))

with(d, plot(Time, SCIprop, col=as.numeric(Treecoverage), pch=19))

with(d3, plot(Days, SCIprop, col=as.numeric(Class), pch=19, main="d3"))

ggplot(data=d, aes(x=Time, y=SCIprop, group=Treecoverage, colour=Treecoverage)) + geom_point() + 
  theme_bw() + geom_point(size = 3)

ggplot(data=d, aes(x=Time, y=SCIprop, shape=Phase, colour=Class)) + geom_point() + 
  theme_bw() + geom_point(size = 3)


#TRANSFORMO TIME
  d3$eTime<-exp(-d3$Time)
  d3$eTime0.1<-exp(-0.1*d3$Time)
  d3$eDays<-exp(-d3$Days)
  d3$eDays0.1<-exp(-0.1*d3$Days)


m.d3.Time<-glmmadmb(SCIprop~Class*Time+Treecoverage+Locality+Start+(1|ID), 
                    zeroInflation=TRUE, family="beta", data=d3)

m.d3.eTime<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+Start+(1|ID), 
                     zeroInflation=TRUE, family="beta", data=d3)

#m.d3.eTime0.1<-glmmadmb(SCIprop~Class*eTime0.1+Treecoverage+Locality+Start+(1|ID), 
#                        zeroInflation=TRUE, family="beta", data=d3)                #error
#m.d3.eTime0.1<-glmmadmb(SCIprop~Class*eTime0.1+Treecoverage+Locality+Start+(1|ID), 
#                        zeroInflation=TRUE, family="gaussian", data=d3)  #no me calcula los residuales

m.d3.eTime0.1<-lmer(SCIprop~Class*eTime0.1+Treecoverage+Locality+Start+(1|ID), 
                      data=d3)  #no me calcula los residuales

m.d3.Days<-glmmadmb(SCIprop~Class*Days+Treecoverage+Locality+Start+(1|ID), 
                    zeroInflation=TRUE, family="beta", data=d3)

m.d3.eDays<-glmmadmb(SCIprop~Class*eDays+Treecoverage+Locality+Start+(1|ID), 
                    zeroInflation=TRUE, family="beta", data=d3)

m.d3.eDays0.1<-lmer(SCIprop~Class*eDays0.1+Treecoverage+Locality+Start+(1|ID), 
                      data=d3)

    par(mfrow=c(2, 6))
      plot(fitted(m.d3.Time), residuals(m.d3.Time), main="m.d3.Time")
        abline(0,0)
        qqnorm(residuals(m.d3.Time), main="qqnorm m.d3.Time")
        qqline(residuals(m.d3.Time))
      plot(fitted(m.d3.eTime), residuals(m.d3.eTime), main="m.d3.eTime")
        abline(0,0)
        qqnorm(residuals(m.d3.eTime), main="qqnorm m.d3.eTime")
        qqline(residuals(m.d3.eTime))
      plot(fitted(m.d3.eTime0.1), residuals(m.d3.eTime0.1), main="m.d3.eTime0.1")
        abline(0,0)
        qqnorm(residuals(m.d3.eTime0.1), main="qqnorm m.d3.eTime0.1")
        qqline(residuals(m.d3.eTime0.1))
      plot(fitted(m.d3.Days), residuals(m.d3.Days), main="m.d3.Days")
        abline(0,0)
        qqnorm(residuals(m.d3.Days), main="qqnorm m.d3.Days")
        qqline(residuals(m.d3.Days))
      plot(fitted(m.d3.eDays), residuals(m.d3.eDays), main="m.d3.eDays")
        abline(0,0)
        qqnorm(residuals(m.d3.eDays), main="qqnorm m.d3.eDays")
        qqline(residuals(m.d3.eDays))
      plot(fitted(m.d3.eDays0.1), residuals(m.d3.eDays0.1), main="m.d3.eDays0.1")
        abline(0,0)
        qqnorm(residuals(m.d3.eDays0.1), main="qqnorm m.d3.eDays0.1")
        qqline(residuals(m.d3.eDays0.1))

#El modelo que mejor se ajusta es m.d3.eTime

m.d3.eTime<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+Start+(1|ID), 
                     zeroInflation=TRUE, family="beta", data=d3)
      summary(m.d3.eTime) #eTime +, ClassIII -, Treecoverage1 -, LocalitySNS -

m.d3.eTime.2<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+(1|ID), 
                     zeroInflation=TRUE, family="beta", data=d3) #the best
      summary(m.d3.eTime.2) #eTime +, ClassIII -, Treecoverage1 -, LocalitySNS -

  AICtab(m.d3.eTime, m.d3.eTime.2) #mejor el m.d3.eTime.2

    #m.d3.eTime.2 también se ajusta bien a las asunciones.

    #Cambio el orden de los niveles del factor Class, para ver si encuentro 
      #diferencias entre Class II y Class III:
        d3$Class=relevel(d3$Class,ref="I")
    #ClassIII es distinta de ClassI y ClassII, pero ClasI y ClassII son iguales

  par(mfrow=c(2, 2))
    plot(fitted(m.d3.eTime), residuals(m.d3.eTime), main="m.d3.eTime")
        abline(0,0)
        qqnorm(residuals(m.d3.eTime), main="qqnorm m.d3.eTime")
        qqline(residuals(m.d3.eTime))
    plot(fitted(m.d3.eTime.2), residuals(m.d3.eTime.2), main="m.d3.eTime.2")
        abline(0,0)
        qqnorm(residuals(m.d3.eTime.2), main="qqnorm m.d3.eTime.2")
        qqline(residuals(m.d3.eTime.2))


#Cambio a soon2, en el que no hay nonatos y sí Vaca3-4 y los que tienen sólo 1 loc:

d2<-read.table("soon2.txt",header=T)
  d2$Treecoverage<-factor(d2$Treecoverage)
  d2$Phase<-factor(d2$Phase)
  d2$InflectionPoint<-factor(d2$InflectionPoint)
  str(d2)
  d2sinT0=subset(d2, d2$Time!=0)
  d2sinT0.sinVacas34=subset(d2sinT0, d2sinT0$ID!="Vacas3_4")
  str(d2sinT0.sinVacas34)

par(mfrow=c(1,1))


ggplot(data=d2sinT0, aes(x=Time, y=SCIprop, shape=Phase, colour=Class)) + geom_point() + 
  theme_bw() + geom_point(size = 3)

#Transformo Time
d2sinT0$eTime<-exp(-d2sinT0$Time)
d2sinT0.sinVacas34$eTime<-exp(-d2sinT0.sinVacas34$Time)

#el mejor modelo anteriormente
m.d2sinT0.eTime.2<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+(1|ID), 
                       zeroInflation=TRUE, family="beta", data=d2sinT0)

m.d2sinT0.eTime.2.sinVacas34<-glmmadmb(SCIprop~Class+eTime+Treecoverage+Locality+(1|ID), 
                            zeroInflation=TRUE, family="beta", data=d2sinT0.sinVacas34)

summary(m.d2sinT0.eTime.2) #eTime +, ClassIII -, Treecoverage1 deja de ser significativo si meto la Vacas 3_4, LocalitySNS -
summary(m.d2sinT0.eTime.2.sinVacas34) #eTime +, ClassIII -, Treecoverage1 -, LocalitySNS -

#Cambio el orden de los niveles del factor Class, para ver si encuentro diferencias entre Class II y Class III:
d2sinT0$Class=relevel(d2sinT0$Class,ref="II")
#ClassIII es distinta de ClassI y ClassII, pero ClasI y ClassII son iguales

par(mfrow=c(2,1))
ggplot(data=d2sinT0, aes(x=Time, y=SCIprop, colour=Treecoverage)) + geom_point() + 
  theme_bw() + geom_point(size = 3)
ggplot(data=d2sinT0.sinVacas34, aes(x=Time, y=SCIprop, colour=Treecoverage)) + geom_point() + 
  theme_bw() + geom_point(size = 3)


#m.d2sinT0.eTime.2 se ajusta bien a las asunciones.
par(mfrow=c(2, 2))
plot(fitted(m.d2sinT0.eTime.2), residuals(m.d2sinT0.eTime.2), main="m.d2sinT0.eTime.2")
abline(0,0)
qqnorm(residuals(m.d2sinT0.eTime.2), main="qqnorm m.d2sinT0.eTime.2")
qqline(residuals(m.d2sinT0.eTime.2))

#m.d2sinT0.eTime.2.sinVacas34 se ajusta mejor sin Vacas
plot(fitted(m.d2sinT0.eTime.2.sinVacas34), residuals(m.d2sinT0.eTime.2.sinVacas34), main="m.d2sinT0.eTime.2.sinVacas34")
abline(0,0)
qqnorm(residuals(m.d2sinT0.eTime.2.sinVacas34), main="qqnorm m.d2sinT0.eTime.2.sinVacas34")
qqline(residuals(m.d2sinT0.eTime.2.sinVacas34))


#Con Vacas34 el modelo se me ajusta peor, y además, treecoverge pierde significación.
#Además, las Vacas34 no me añaden info en cuanto al análisis del cambio de Phase2 a 3.
#Jennifer dice que también vé mejor no incluir Vaca3-4

#IMPORTANTE que en los primeros modelos de prueba meta Age, para ver que no tiene efecto, y poderlo discutir.




#PHASES

phase2<-read.table("StartPhase2.txt",header=T)
str(phase2)

#VR Time ó Days:
plot(phase2$Class, phase2$Days, ylab="Days", xlab="Class", main="All")
m=glm(Days~Class+Locality, data=phase2, family="poisson")
summary(m)
  #ClassI diferente a ClassIII, pero igual a ClassII.
    #Locality tiene efecto, pero creo que es porque Class I y Class III sólo están
    #presentes en RBD. 
    #No puedo obviar Locality. Si lo hago sin las de SNS, me sale lo mismo. Si 
    #obvío locality me salen cosas distintas (las tres clases diferentes, que es lo 
    #que se vé en el boxplot Days~Class):

phase2RBD=subset(phase2, phase2$Locality=="RBD")
plot(phase2RBD$Class, phase2RBD$Days, ylab="Days", xlab="Class", main="Only RBD")
m3=glm(Days~Class, data=phase2RBD, family="poisson")
summary(m3)
    #Sin locality me saldrían diferencias entre las tres clases, pero porque las de SNS 
    #(todas ClassII) tienen valores inferiores de tiempo. Los tejidos blandos se degradan
    #antes en SNS.

#m2=glm(Days~Class, data=phase2, family="poisson") #no es correcto
#summary(m2)


#VR SCI:
plot(phase2$Class, phase2$SCIprop, ylab="SCIprop", xlab="Class", main="all")
x=glmmadmb(SCIprop~Class+Locality, data=phase2, family="beta")
summary(x)
#glm con family beta no me deja. No sé si glmmadmb sería correcto. 
#No dif entre clases. Sí debido a Locality.

#Analizando sólo los datos de RBD:
plot(phase2RBD$Class, phase2RBD$SCIprop, ylab="SCIprop", xlab="Class", main="Only RBD")
x2=glmmadmb(SCIprop~Class, data=phase2RBD, family="beta")
summary(x2)
#No dif en SCIprop segun Class.

#Seguramente debería eliminar la ClassIII, porque es N=1.

#Está bien. Puedo definir bien las clases según el SCI, y veo que no hay diferencias
#entre clases en el momento en que se degradan los tejidos blandos (cambio a Phase2).
  
#y si trabajo con un rango (último punto phase anterior y primer punto phase
#posterior) en vez de con un dato (inicio de phase)??


phase3<-read.table("StartPhase3.txt",header=T)
str(phase3)

#No tiene sentido hacerlo ahora mismo, ya que no hay classI, y ClassIII N=1.
plot(phase3$Class, phase3$Days, ylab="Days", xlab="Class", main="All")







plot(d3$Phase, d3$Time, ylab="Time", xlab="Phase")
plot(d3$Phase, d3$Days, ylab="Days", xlab="Phase")
plot(d3$Time, d3$Phase, ylab="Phase", xlab="Time")
d3$Phase=relevel(d3$Phase,ref="1")

m=glmmadmb(Days~Phase+Class+(1|ID), data=d3, family="Poisson") 
summary(m)
d3$Phase=relevel(d3$Phase,ref="1")

m2=glmmadmb(Days~Phase+(1|ID), data=d3, family="Poisson") 
summary(m2)

d3$Phase=relevel(d3$Phase,ref="3") #Cuando hago relevel, creoq ue no me calcula bien el summary.

http://en.wikipedia.org/wiki/List_of_probability_distributions

dPhase1=subset(d3, d3$Phase==1)



str(d3)

boxplot( ~ Time | Phase, data = d3)


coplot(Time~Phase | Class, data = d3, pch =19)





#________________

  AICtab(m.d.Time, m.d.eTime, m.d.eTime0.1, m.d3.Time, m.d3.eTime) 
  #De mejor a peor AIC: m.d.eTime0.1, m.d.eTime, m.d.Time, m.d3.Time, m.d3.eTime
       #Mejor AIC con los time=0. Por qué? porque tiene más casos, quizás.
  #Pero el mejor ajuste es el de m.d3.eTime, el que tiene peor AIC.

  par(mfrow=c(3, 4))
    plot(fitted(m.d3.Time), residuals(m.d3.Time), main="m.d3.Time")
      abline(0,0)
      qqnorm(residuals(m.d3.Time), main="qqnorm m.d3.Time")
      qqline(residuals(m.d3.Time))
    plot(fitted(m.d3.eTime), residuals(m.d3.eTime), main="m.d3.eTime")
      abline(0,0)
      qqnorm(residuals(m.d3.eTime), main="qqnorm m.d3.eTime")
      qqline(residuals(m.d3.eTime))


summary(m.d.Time) #Time -, ClassIII*Time -, SNS - 
summary(m.d.eTime) #eTime +, ClassIII -, SNS -, Treecov1 marg-, WST marg-
summary(m.d.eTime0.1) #eTime0.1 +, ClassII -, ClassIII -, ClassII*eTime0.1 +, ClassIII*eTime0.1 +, SNS -, Treecov1 -, WST +
summary(m.d3.Time) #Time -, SNS marg-, Treecov1 -
summary(m.d3.eTime) #eTime +, ClassIII -, SNS -, Treecov1 -




  #por la forma de los modelos, creo que es mejor eliminando los casos en que time=0.
  #Hoy 18/3/15, un día más tarde, no me deja calcular el modelo si lo trabajo 
    #con eTime en la base de datos sin Time=0. Ni siquiera sacándole el 
    #zeroInflation. Desesperante. No se por qué.
    #Igualmente, las significaciones son iguales con o sin Time=0.

  #En m5 (era con eTime0.1 y sin Time=0) me salían diferencias significativas 
  #en todas las variables:
      #Class II -
      #Class III -
      #eTime +
      #Treecoverage1 -
      #LocalitySNS -  
      #StartWST +
      #ClassII:eTime +
      #ClassIII:eTime +
  
  #Cambio el orden de los niveles del factor Class, para ver si encuentro diferencias entre Class II y Class III:
    d2tree.sintime0$Class=relevel(d2tree.sintime0$Class,ref="II")

    m5b<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+Start+(1|ID), 
                      zeroInflation=TRUE, family="beta", data=d2tree.sintime0) # ESTE
      summary(m5b)
  #Si, hay diferencias entre las tres clases, Class I, II y III. Aunque la interacción Tiempo con Class no es distinto
  #para Class II y III. O sea, el punto de origen es distinto para Class II y III, pero no la pendiente que siguen.


par(mfrow=c(2,1))
with(d2tree.sintime0, plot(eTime, SCIprop, col=as.numeric(Class), pch=19, main="d2tree.sintime0"))
with(d2tree, plot(eTime,SCIprop, col=as.numeric(Class), pch=19, main="d2tree"))


#Cuanto más tiempo pasa, más degradación.
#Class I tiene el mayor intercepto con SCI
#Class II tiene un intercepto intermedio
#Class III tiene el menor intercepto
#Las tres son distintas. Cuanta menor masa inicial (mayor Class), más degradación hay.
#La pendiente de degradación, es menos pronunciada para la Class I (mayor masa) que para las Class II y III.
  ggplot(data=d2tree.sintime0, aes(x=Time, y=SCIprop, group=Class, colour=Class)) + geom_point() + 
    theme_bw() + geom_point(size = 3)

#Con covertura arbórea el SCI es menor, o sea que hay más degradación.
  ggplot(data=d2tree.sintime0, aes(x=Time, y=SCIprop, group=Treecoverage, colour=Treecoverage)) + geom_point() + 
    theme_bw() + geom_point(size = 3)
  plot(d2tree.sintime0$Treecoverage, d2tree.sintime0$SCIprop, xlab="Treecoverage", ylab="SCIprop")

#En SNS el SCI es menor, o sea hay más degradación.
  ggplot(data=d2tree.sintime0, aes(x=Time, y=SCIprop, group=Locality, colour=Locality)) + geom_point() + 
    theme_bw() + geom_point(size = 2)
  plot(d2tree.sintime0$Locality, d2tree.sintime0$SCIprop, xlab="Locality", ylab="SCIprop")

#Los cuerpos que no tienen tejidos blandos presentan menores mayor proporción de cuerpo, o sea menor degradación.
#Podría ser porque empiezan más tarde están en una fase posterior (en Phase II) del proceso de degradación, 
#por lo que el descenso del SCI se produce más lentamente.
  ggplot(data=d2tree.sintime0, aes(x=Time, y=SCIprop, group=Start, colour=Start)) + geom_point() + 
    theme_bw() + geom_point(size = 3)
  plot(d2tree.sintime0$Start, d2tree.sintime0$SCIprop, xlab="Start", ylab="SCIprop")

ggplot(data=d2tree.sintime0, aes(x=Time, y=Lossrate, group=ID, colour=Class)) + geom_point() + geom_line() +
  theme_bw() + geom_point(size = 1) + geom_line(size = 1)

ggplot(data=d2tree.sintime0, aes(x=Time, y=SCIprop, group=ID, colour=Class)) + geom_point() + geom_line() +
  theme_bw() + geom_point(size = 1) + geom_line(size = 1)



str(d2tree.sintime0)

#Creo que no estoy interpretando bien los resultados. Creo que no es degradación. 
#Podría ser tasa de degradación?? No, porque el tiempo no está metido en la VR.

#Me parece que no esta ploteando bien los daots sin time=0. Crear yo una nueva base de datos.
#Probándolo con una base de datos nueva, en realidad me sale igual. Creo que sí me lo estaba cogiendo bien, pero 
#como hay tantos puntos en Tiempo=0 y SCIprop=1, se superponen y en realidad sólo se vé uno al final.


m4.conTime<-glmmadmb(formula = SCIprop ~ Class * Time + Treecoverage + Locality + 
           Start + (1 | ID), data = dtree, family = "beta", zeroInflation = TRUE)
summary(m4.conTime)
#+ tiempo, - SCI, + degradación
#sí Treecoverage, - SCI, + degradación
#SNS, - SCI, + degradación
#WST, + SCI, - degradación
#Interacción tiempo classII, más pronunciada a menos SCI, + tasa degradación
#Interacción tiempo classIII, la más pronunciada a menos SCI, ++ tasa degradación




#Sin los nonatos (Vaca5 y Vaca6). No me mejora.
dtree.sinVaca5=subset(dtree,dtree$ID!="Vaca5")
dtree.sinnonatos=subset(dtree.sinVaca5, dtree.sinVaca5$ID!="Vaca6")

m3.sinnonatos<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+(1|ID), zeroInflation=TRUE, family="beta", data=dtree.sinnonatos)
summary(m3.sinnonatos)
AICtab(m3, m3.sinnonatos) 

par(mfrow=c(2,2))
plot(fitted(m3), residuals(m3))
abline(0,0)
qqnorm(residuals(m3), main="qqnorm m3")
qqline(residuals(m3))
plot(fitted(m3.sinnonatos), residuals(m3.sinnonatos))
abline(0,0)
qqnorm(residuals(m3.sinnonatos), main="qqnorm m3.sinnonatos")
qqline(residuals(m3.sinnonatos))
par(mfrow=c(1,2))
with(dtree, plot(Time,SCIprop, col=as.numeric(Class), pch=19))
with(dtree.sinnonatos, plot(Time,SCIprop, col=as.numeric(Class), pch=19))
with(d2tree, plot(Time,SCIprop, col=as.numeric(Class), pch=19))


#Comparo el efecto de los buitres (cubierto o no) en SNS: sí hay efecto. Si está cubierto, la degradación es menor. 
#Sin diferencias entre la clase II y III!!

dtreeSNS=subset(dtree,dtree$Locality=="SNS")

mt<-glmmadmb(SCIprop~Class*eTime+Treecoverage+(1|ID), zeroInflation=TRUE, family="beta", data=dtreeSNS)
summary(mt)
mr<-glmmadmb(SCIprop~Class+eTime+Treecoverage+(1|ID), zeroInflation=TRUE, family="beta", data=dtreeSNS) #este es un 
#buen modelo (no tiene RBD ni clase I, la interacción entre clase y tiempo no existe)
summary(mr)
plot(fitted(mr), residuals(mr))
abline(0,0)
qqnorm(residuals(mr), main="qqnorm mr")
qqline(residuals(mr))


#Sólo con los datos de RBD, en que tengo Date y Location.
d2tree.RBD=subset(d2tree,d2tree$Locality=="RBD")
d2tree.RBD.Date=subset(d2tree.RBD,d2tree.RBD$Date!="NA")
d2tree.RBD.Date.Loc=subset(d2tree.RBD.Date,d2tree.RBD.Date$Location!="NA")
str(d2tree.RBD.Date.Loc)
d2tree.RBD.Date.Loc$eTime<-exp(-0.1*d2tree.RBD.Date.Loc$Time)

mRBD<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Date+Location+(1|ID), 
               zeroInflation=TRUE, family="beta", data=d2tree.RBD.Date.Loc) #no sé por qué, pero no me lo calcula si 
#le meto Start
summary(mRBD)
#No hay diferencias entre fechas, ni entre localizaciones.
#No vale la pena incluir estas VE.




#Anterior_______________


#Manuela:
#Mirar de obtener el tiempo en días, no en fracción. Éso quizá te cambie la distribución de los datos.
#Si uso eTime0.1, los datos (las observaciones reales) ya no siguen una distribución beta (en la que 
#hay puntos entre 0 y 1 y que siguen una distribución curvada), y sí lineal, por lo que debo usar
#una gaussian en vez de una beta. Lo puedo hacer con el glmmadmb o con el lmer.
#Vé mejor que haya eliminado los casos en que Tiempo=0.

#PLOTS
ggplot(data=d, aes(x=Time, y=SCI, group=Treecoverage, colour=Treecoverage)) + geom_line() + geom_point() + 
  theme_bw() + geom_point(size = 2) +  geom_line(size = 0.2)  #por clase

with(d, plot(Time,SCIprop, col=as.numeric(Treecoverage), pch=19))
#Las que no están tapadas se degradan más lentamente?? Y las que están tapadas más rápido?

plot(d$Class, d$Treecoverage, xlab="Class", ylab="Treecoverage", main="Treecoverage")
#Creo que tiene que ver con que las clase I y clase II estan sobrerepresentadas entre las no tapadas, y las clase III
#lo están entre las tapadas.

par(mfrow=c(2,2))
with(d, plot(Time,SCIprop, col=as.numeric(Class), pch=19))

plot(d$Locality, d$Treecoverage, xlab="Locality", ylab="Treecoverage", main="Treecoverage")
#No parece que haya mucha correlación entre el área y el treecoverage.

par(mfrow=c(1,1))
with(d3, plot(Time,SCIprop, col=as.numeric(Class), pch=19, main="d3"))
with(d3, plot(eTime0.1, SCIprop, col=as.numeric(Class), pch=19, main="d3"))
with(d3, plot(Days, SCIprop, col=as.numeric(Class), pch=19, main="d3"))

#para meter Treecoverage tengo que sacar manualmente los casos en que tengo NA (sólo es 1: 5NAVAS)
#  dtree=subset(d,d$Treecoverage!="NA")
#  d2tree=subset(d2,d2$Treecoverage!="NA")


#m1<-glmmadmb(SCIprop~Class*Time+Treecoverage+Locality+(1|ID), 
#             zeroInflation=TRUE, family="beta", data=d)
#      summary(m1)
#Hay diferencias entre Locality!!! Hay más degradación en SNS (donde hay más buitres!)
#Diferencias entre Locality: en SNS hay buitres, en RBD hay jabalíes.
#Y parece que Treecoverage también tiene algo de significación! Hay menos SCI, más degradación cuando sí hay 
#cobertura arbórea. Según Eloísa: los buitres no pueden degradarlo.

#        plot(fitted(m1), residuals(m1))
#          abline(0,0)
#        qqnorm(residuals(m1), main="qqnorm m1")
#          qqline(residuals(m1))

#Sin Locality
#m2<-glmmadmb(SCIprop~Class*Time+Treecoverage+(1|ID), zeroInflation=TRUE, family="beta", data=dtree)
#  summary(m2)
#  plot(fitted(m2), residuals(m2))
#    abline(0,0)
#  qqnorm(residuals(m2), main="qqnorm m2")
#    qqline(residuals(m2))




#Hay 3 puntos outliers: Gamo1, 12NAVAS (tejón) y 9 NAVAS (meloncillo).
#El tejón y el meloncillo Eloísa dice que son casos extremos, se momificaron. 
#Intuyen que puede ser porque murieran por envenenamiento, y por éso no se los comieran.
#Con el Gamo1 no sé qué pudo pasar.
#Lo mejor será quitar estos casos de la base de datos.
#Así, añado 4 casos más a la base inicial de Eloísa (son todos ClassII, 3 de ellos de 
#RBD)
#Creo que lo mejor va a ser que deje sólo los casos que ella tenía, sin añadir 
#estos 4. Porque se me añaden dos puntos a los 63 meses con un SCI más alto, y porque 
#así es más fácil de explicar en el texto, diciendo que se trabaja solo con aquellas
#carcasas que se muestrearon más de una vez, aunque luego no use el dato de Tiempo=0.


#m3<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+(1|ID), 
#             zeroInflation=TRUE, family="beta", data=d)
#  summary(m3)
#    par(mfrow=c(1,2))
#    plot(fitted(m3), residuals(m3))
#      abline(0,0)
#    qqnorm(residuals(m3), main="qqnorm m3")
#      qqline(residuals(m3))
#  AICtab(m1, m3) 
#m3 (con exp(-0.1*Time)) mejor modelo que m1 (con Tiempo)


#m.d.Time<-glmmadmb(SCIprop~Class*Time+Treecoverage+Locality+Start+(1|ID), 
#             zeroInflation=TRUE, family="beta", data=d)

#m.d.eTime<-glmmadmb(SCIprop~Class*eTime+Treecoverage+Locality+Start+(1|ID), 
#                    zeroInflation=TRUE, family="beta", data=d)

#m.d.eTime0.1<-glmmadmb(SCIprop~Class*eTime0.1+Treecoverage+Locality+Start+(1|ID), 
#             zeroInflation=TRUE, family="beta", data=d)


#El Start tiene efecto! Los que empiezan sin tejido blando (WST) se degradan más rápido.

#  AICtab(m.d.Time, m.d.eTime0.1, m.d.eTime)
#De mejor a peor AIC: m.d.eTime0.1, m.d.eTime, m.d.Time
#De mejor a peor por la forma del ajuste: m.d.eTime

#  par(mfrow=c(3,2))
#    plot(fitted(m.d.Time), residuals(m.d.Time), main="m.d.Time")
#      abline(0,0)
#    qqnorm(residuals(m.d.Time), main="qqnorm m.d.Time")
#      qqline(residuals(m.d.Time))
#    plot(fitted(m.d.eTime), residuals(m.d.eTime), main="m.d.eTime")
#      abline(0,0)
#    qqnorm(residuals(m.d.eTime), main="qqnorm m.d.eTime")
#      qqline(residuals(m.d.eTime))
#    plot(fitted(m.d.eTime0.1), residuals(m.d.eTime0.1), main="m.d.eTime0.1")
#      abline(0,0)
#    qqnorm(residuals(m.d.eTime0.1), main="qqnorm m.d.eTime0.1")
#      qqline(residuals(m.d.eTime0.1))


#Pruebo trabajar sin las primeras localizaciones, en que tiempo=0 y SCI=100, para tener una varianza de los datos
#más homogénea, y que pueda ajustar mejor los modelos.
#Todas las carcasas tienen SCI=100 a Time=0.

#d.sintime0=subset(d,d$Time!="0")
#d2tree.sintime0=subset(d2tree,d2tree$Time!="0")
#d2tree.sintime0$eTime<-exp(-0.1*d2tree.sintime0$Time)


m3x<-glmmadmb(SCIprop~Class*eTime+Treecoverage+(1|Locality)+(1|ID), zeroInflation=TRUE, family="beta", data=dtree)
  summary(m3x)
  plot(fitted(m3x), residuals(m3x))
    abline(0,0)
  qqnorm(residuals(m3x), main="qqnorm m3x")
    qqline(residuals(m3x))
#Se ajusta igual de mal que m3, pero tiene sig marginal Class II, mayor Log-likelihood que m3 y menor AIC que m3.
#Quizá este sea mejor??
#Creo que las localidades en realidad me interesan más como en el m3, porque quiero ver si hay diferencias entre ésas 
#localidades en concreto, no cualquier localidad.
AICtab(m1, m3) 

#Según http://glmmadmb.r-forge.r-project.org/glmmADMB.html, me es mejor que m1

m4<-glmmadmb(SCIprop~Class*eTime+Locality+(1|ID), zeroInflation=TRUE, family="beta", data=dtree) #error. No me lo calcula.
d$eTime<-exp(-0.1*d$Time)
m4<-glmmadmb(SCIprop~Class*eTime+Locality+(1|ID), zeroInflation=TRUE, family="beta", data=d) #error. No me lo calcula.
#No me importa mucho, porque treecoverage si dá sig. marginal.



m1<-glmmadmb(SCIprop~Class*Time+(1|ID)+(1|Sps), zeroInflation=TRUE, family="beta", data=d)
m2<-glmmadmb(SCIprop~Class+Time+Sps+Start+(1|ID), zeroInflation=TRUE, family="beta", data=d)
mx<-glmmadmb(SCIprop~Class*Time+Sp+(1|ID), zeroInflation=TRUE, family="beta", data=d)
my<-glmmadmb(SCIprop~Class*Time+(1|ID), zeroInflation=TRUE, family="beta", data=d)
#Creo que este sería el más acertado, en cuanto a qué VE meto y como
#El modelo no me mejora controlando o no por Sps
#por lo que me contaron, ni Sp ni edad no tiene por qué influir en la degradación. 


m3<-glmmadmb(SCIprop~Class+Time+Start+(1|ID), zeroInflation=TRUE, family="beta", data=d)

ma<-glmmadmb(SCIprop~Class+Time+(1|ID), zeroInflation=TRUE, family="beta", data=d)
summary(ma)


m4<-glmmadmb(SCIprop~Class*Time+(Sps|ID), zeroInflation=TRUE, family="beta", data=d) #dá error
m5<-glmmadmb(SCIprop~Class*Time+(1|ID)+(Sps|ID), zeroInflation=TRUE, family="beta", data=d) #dá error



#Mirar correlación de las VE

m<-glmmadmb(SCIprop~Class*Time+Treecoverage+(1|ID), zeroInflation=TRUE, family="beta", data=d)



#_____________________________________________________________________________________________________________
#Cuatro modelos probando con Mass en vez de Class. No quieren hacerlo con Mass, quieren hacerlo con class, porque 
#son clases predefinidas en otros trabajos.

  m6<-glmmadmb(SCIprop~Mass*Time+(1|ID)+(1|Sps), zeroInflation=TRUE, family="beta", data=d) 
  
  m7<-glmmadmb(SCIprop~Mass*Time+Sps+(1|ID), zeroInflation=TRUE, family="beta", data=d)

  m8<-glmmadmb(SCIprop~Mass+Time+Sps+(1|ID), zeroInflation=TRUE, family="beta", data=d)

  m9<-glmmadmb(SCIprop~Mass*Time+(1|ID), zeroInflation=TRUE, family="beta", data=d)

  m10<-glmmadmb(SCIprop~Mass*Time+(1|Sps), zeroInflation=TRUE, family="beta", data=d)

      summary(m6)  

  #m6 mejor que m7 (por AIC)
  #Controlando por Sps o por ID me dá lo mismo! (m9vs m10)

  #Los 5 modelos estan igual de mal ajustados.
  par(mfrow=c(3,2))
    plot(fitted(m6), residuals(m6))
      abline(0,0)
    plot(fitted(m7), residuals(m7))
      abline(0,0)
    plot(fitted(m8), residuals(m8))
      abline(0,0)
    plot(fitted(m9), residuals(m9))
      abline(0,0)
  plot(fitted(m10), residuals(m10))
      abline(0,0)

  par(mfrow=c(3,2))
    qqnorm(residuals(m6))
      qqline(residuals(m6))
    qqnorm(residuals(m7))
      qqline(residuals(m7))
    qqnorm(residuals(m8))
      qqline(residuals(m8))
    qqnorm(residuals(m9))
      qqline(residuals(m9))
    qqnorm(residuals(m10))
      qqline(residuals(m10))
#_____________________________________________________________________________________________________________


#m9 vs my: con Class el AIC es menor
  summary(m9)
  summary(my)  


#Como veo que los gráficos me muestran que los datos siguen una distribución exponencial negativa,
#voy a intentar transformar la VE Tiempo, para linealizar los datos y para que los modelos 
#se ajusten más

par(mfrow=c(2,2))
with(d, plot(Time,SCIprop))
with(d, plot(exp(-0.2*Time),SCIprop))
qqnorm(exp(-0.2*Time))
qqline(with(d, plot(exp(-0.2*Time),SCIprop)))

with(d, plot(exp(-0.1*Time),SCIprop))
with(d, plot(exp(-Time),SCIprop))

with(d, plot(Time,SCIprop))
with(d, plot(Time,sqrt(SCIprop)))
with(d, plot(Time,(SCIprop^(0.25))))
with(d, plot(Time,exp(-SCIprop)))
with(d, plot(Time,exp(-0.1*SCIprop)))
with(d, plot(Time,exp(-5*SCIprop)))


qqnorm(d$SCIprop)
qqline(d$SCIprop)

qqnorm(exp(-d$SCIprop))
qqline(exp(-d$SCIprop))

qqnorm(exp(-0.1*d$SCIprop))
qqline(exp(-0.1*d$SCIprop))


#Pruebo a cambiar la VE Time -> Me dá error. Parece que no puede trabajar con ésa transformación.
  with(d, plot(exp(-0.1*Time),SCIprop))
  with(d, plot(Mass,SCIprop, col=as.numeric(ID), pch=19))
 
  d$eTime<-exp(-0.1*d$Time)

  m<-glmmadmb(SCIprop~Class*eTime+(1|ID)+(1|Sps), zeroInflation=TRUE, family="beta", data=d) #error
  m<-glmmadmb(SCIprop~Class+eTime+(1|ID), zeroInflation=TRUE, family="beta", data=d) #error
summary(m)
  m<-glmmadmb(SCIprop~Class*eTime+(1|ID), zeroInflation=TRUE, family="beta", data=d) #error

  m<-glmmadmb(SCIprop~Mass*eTime+(1|ID)+(1|Sps), zeroInflation=TRUE, family="beta", data=d) #error
  m<-glmmadmb(SCIprop~Mass+eTime+(1|ID), zeroInflation=TRUE, family="beta", data=d) #error
  m<-glmmadmb(SCIprop~Mass*eTime+(1|ID), zeroInflation=TRUE, family="beta", data=d) #error
  
summary(m)

d$eTime<-exp(-d$Time)


#Manuela dice:
#No se le ocurre cómo mejorar el modelo.
#Puedo intentar algún otro link.
#Puedo modificar alguna variable, como tiempo. Según lo que veo, parecería que los datos no siguen
#un modelo lineal, sino exponencial negativo. Si uso e^-tiempo puedo hacerlo exponencial negativo.
#Lo he provado, pero no mejora. También lo hemos probado multiplicando por un coeficiente, y
#al * por 0.1 el tiempo sale algo mejor el modelo. Pero Manuela dice que es un poco cutresalchichero.
#Dice que use masa, en vez de clase.
#No se le ocurre cómo mejorarlo. Dice que se debe mucho a que mis datos son un poco chungos. Tengo 
#muchos para una clase y muy pocos para las otras.
#Además, según el individuo, también tengo más o menos datos, y muchos más datos tomados al poco
#tiempo de morir y pocos datos después de mucho tiempo de morir.

#Dice que quizá también pueda intentarlo cambiando la función link.


m11<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="log", zeroInflation=TRUE,  data=d) #error

m12<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="probit", zeroInflation=TRUE,  data=d)
#mal ajuste
    summary(m12)
    par(mfrow=c(1,1))
    plot(fitted(m12), residuals(m12))
    abline(0,0)
    qqnorm(residuals(m12))
    qqline(residuals(m12))

m13<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="identity", zeroInflation=TRUE,  data=d) #error
m13<-glmmadmb(SCIprop~Mass+Time+(1|ID), family="beta", link="identity", zeroInflation=TRUE,  data=d) #error

m14<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="cloglog", zeroInflation=TRUE,  data=d) #error
m14b<-glmmadmb(SCIprop~Mass+Time+(1|ID), family="beta", link="cloglog", zeroInflation=TRUE,  data=d) 
m14c<-glmmadmb(SCIprop~Class+Time+(1|ID), family="beta", link="cloglog", zeroInflation=TRUE,  data=d) 
summary(m14c)
par(mfrow=c(1,1))
plot(fitted(m14c), residuals(m14c))
abline(0,0)
qqnorm(residuals(m14c))
qqline(residuals(m14c))
qqnorm(residuals(m14b))
qqline(residuals(m14b))


m15<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="sqrt", zeroInflation=TRUE,  data=d) #error
#mal escrita la link

m16<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="inverse", zeroInflation=TRUE,  data=d) #error
m16<-glmmadmb(SCIprop~Mass+Time+(1|ID), family="beta", link="inverse", zeroInflation=TRUE,  data=d) #error

m17<-glmmadmb(SCIprop~Mass*Time+(1|ID), family="beta", link="loglog", zeroInflation=TRUE,  data=d) #error

dC$Class=relevel(dC$Class, ref="2005")  #Cambias el orden de los factores

?family
??link




##INNECESARIO
##SCI en las N=30 carcasses recién encontradas (C en Start)________________________________________

dC=subset(d,d$Start=="C")
dC
str(dC)

ggplot(data=dC, aes(x=Time, y=SCI, group=Class, colour=Class)) + geom_point() + 
  theme_bw() + geom_point(size = 2)   #por clase

m<-glmmadmb(SCIprop~Class+Time+Sps+(1|ID), family="beta", data=dC)
summary(m)

m<-glmmadmb(SCIprop~Class+Time+(1|ID), family="beta", data=dC) #por defecto en R, la distribución 
#beta usa el link logit. logit: g(x) = log [(x/(1-x)]
#Me parece que puede usar otras (ver pdf "beta regresion"), pero esta es la que usa por defecto. 
summary(m)


m<-betareg(SCIprop~Class, data=dC)
summary(m)
plot(m)




#There are some kinds of proportion data, such as percentage cover, which are best analysed using conventional 
#models (normal errors and constant variance) following arcsine transformation. The response variable, y, 
#measured in radians,is sin???1???001×p, where p is percentage cover.
#asin(sqrt(d$SCI/100))


Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  logBroodSize=log(BroodSize),
                  NCalls=SiblingNegotiation)

fit_zipoiss <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+
                          offset(logBroodSize)+(1|Nest),
                        data=Owls,
                        zeroInflation=TRUE,
                        family="poisson")
summary(fit_zipoiss)
