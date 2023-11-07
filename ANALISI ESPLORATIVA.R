
sani = read.csv('HOP_controlND_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
sani = sani[, -1]
chirurgia = read.csv('HOP_RI_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
chirurgia = chirurgia[, -1]
fisioterapia = read.csv('HOP_TPPI_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
fisioterapia = fisioterapia[,-1]

#### PROCESSING ####

# nei file qual salvo le variaibli qualitative tipo codice identificativo paziente

#PROCESSING SANI
df = NULL
nomi = colnames(sani)

for (i in 1:dim(sani)[2]) {
  if (sani[1, i] == 'Knee' && sani[4, i]=='X' && grepl('R_01', nomi[i], fixed=T)) {
    df = cbind(df, sani[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.sani = as.data.frame(cbind(colnames(df), t(df)))

# PROCESSING CHIRURGIA
df = NULL
nomi = colnames(chirurgia)

for (i in 1:dim(chirurgia)[2]) {
  if (chirurgia[1, i] == 'Knee' && chirurgia[4, i]=='X' && grepl('R_01', nomi[i], fixed=T)) {
    df = cbind(df, chirurgia[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.chir = as.data.frame(cbind(colnames(df), t(df)))


# PROCESSING FISIOTERAPIA
df = NULL
nomi = colnames(fisioterapia)

for (i in 1:dim(fisioterapia)[2]) {
  if (fisioterapia[1, i] == 'Knee' && fisioterapia[4, i]=='X' && grepl('R_01', nomi[i], fixed=T)) {
    df = cbind(df, fisioterapia[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.fisio = as.data.frame(cbind(colnames(df), t(df)))

sani=dati.sani[, 6:1682]
chir=dati.chir[, 6:2229]
fisio=dati.fisio[, 6:2010]

qual.sani = dati.sani[ , 1:5]
qual.chir = dati.chir[ , 1:5]
qual.fisio = dati.fisio[ , 1:5]


rm(nomi, df, i, fisioterapia, chirurgia, dati.chir, dati.fisio, dati.sani)


#### ANALISI ESPLORATIVA SANI ####

sani = as.matrix(sani)
sani = t(sani)
sani=='NaN' # per vedere quali mancano
sani = sani[1:673]
sani = as.numeric(sani)

x11()
matplot(sani, main='Individui sani', type='l', col='darkgreen', lwd=2)

x11()
boxplot(sani, col='darkgreen')

graphics.off()



#### ANALISI ESPLORATIVA CHIRURGIA ####

chir = t(chir)
chir = as.matrix(chir)

tc = NULL
for (i in 1:23) {
tc = c(tc, 2224 - sum(chir[, i]==NaN)) }

# tutti quelli oltre 1729 compreso sono NaN
chir = chir[1:1728, ]
chir = matrix(as.numeric(chir), ncol = ncol(chir))


x11()
matplot(chir, main='Individui sottoposti a chirurgia', type='l', lwd=2)

# vedere se qualcuno ha dati mancanti nel mezzo
for ( i in 1:23) {
print(sum(chir[1:tc[i], i]=='NaN')) }
# solo il paziente 5 ne ha 79
  


#### ANALISI ESPLORATIVA FISIOTERAPIA ####

fisio = t(fisio)
fisio = as.matrix(fisio)

tf = NULL
for (i in 1:23) {
  tf = c(tf, 2005 - sum(fisio[, i]==NaN)) }


# tutti quelli oltre 1246 compreso sono NaN
fisio = fisio[1:1245, ]
fisio = matrix(as.numeric(fisio), ncol = ncol(fisio))

x11()
matplot(fisio, main='Individui sottoposti a fisioterapia', type='l', lwd=2)

# vedere se qualcuno ha dati mancanti nel mezzo
for ( i in 1:23) {
  print(sum(fisio[1:tf[i], i]=='NaN')) }
# solo il paziente 16 ne ha 86


#### INDICI DEI VARI GRUPPI ####

mean(sani)
median(sani)
quantile(sani, 0.25)
quantile(sani, 0.75)

# PER FARE GLI INDICI, DEVO TOGLIERE I NAN -> salvo in aux i dati senza NaN e poi la elimino
aux = chir[!is.nan(chir)]
mean(aux)
median(aux)
quantile(aux, 0.25)
quantile(aux, 0.75)

aux = fisio[!is.nan(fisio)]
mean(aux)
median(aux)
quantile(aux, 0.25)
quantile(aux, 0.75)

#### BOXPLOT COMPARATI ####

# SANI VS CHIRURGIA (diviso in due perché troppo grosso se no)
x11()
par(mfrow=c(2,6))
boxplot(sani, col='darkgreen')
for ( i in 1:11) {
boxplot(chir[,i], col='firebrick2') 
}

x11()
par(mfrow=c(2,6))
for ( i in 12:23) {
  boxplot(chir[,i], col='firebrick2') 
}

# SANI VS FISIOTERAPIA (diviso in due perché troppo grosso se no)
x11()
par(mfrow=c(2,6))
boxplot(sani, col='darkgreen')
for ( i in 1:11) {
  boxplot(fisio[,i], col='dodgerblue') 
}

x11()
par(mfrow=c(2,6))
for ( i in 12:23) {
  boxplot(fisio[,i], col='dodgerblue') 
}

