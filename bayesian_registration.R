
sani = read.csv('HOP_controlND_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
sani = sani[, -1]
chirurgia = read.csv('HOP_RI_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
chirurgia = chirurgia[, -1]
fisioterapia = read.csv('HOP_TPPI_ThoraxPelvisHipKneeFoot_all.csv', sep=';')  
fisioterapia = fisioterapia[,-1]

#### PROCESSING ####

# PROCESSING SANI
df = NULL
nomi = colnames(sani)

for (i in 1:dim(sani)[2]) {
  if (sani[1, i] == 'Knee' && sani[4, i]=='X' && grepl('_01', nomi[i], fixed=T)) {
    df = cbind(df, sani[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.sani = as.data.frame(cbind(colnames(df), t(df)))

sani = as.data.frame(t(dati.sani[, 6:dim(dati.sani)[2]]))
qual.sani = as.data.frame(t(dati.sani[, 1:5])) # salvo le variabili qualitative che non si sa mai


# salvo il numero di osservazioni per ogni paziente sano
ts = 1:dim(sani)[2]

for (i in 1:length(ts)) {
  ts[i] = sum(sani[,i]!='NaN') + sum(sani[1:sum(sani[,i]!='NaN'), i]=='NaN')
}
sani = sani[1:max(ts), ]

# vedere se qualcuno ha dati mancanti nel mezzo
for ( i in 1:dim(sani)[2]) {
  print(sum(sani[1:ts[i], i]=='NaN')) }
# i sani non hanno dati mancanti nel mezzo


# PROCESSING CHIRURGIA
df = NULL
nomi = colnames(chirurgia)

for (i in 1:dim(chirurgia)[2]) {
  if (chirurgia[1, i] == 'Knee' && chirurgia[4, i]=='X' && grepl('_01', nomi[i], fixed=T)) {
    df = cbind(df, chirurgia[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.chir = as.data.frame(cbind(colnames(df), t(df)))

chirurgia = as.data.frame(t(dati.chir[, 6:dim(dati.chir)[2]]))
qual.chir = as.data.frame(t(dati.chir[, 1:5])) # salvo le variabili qualitative che non si sa mai


# salvo il numero di osservazioni per ogni paziente di chirurgia
tc = 1:dim(chirurgia)[2]

for (i in 1:length(tc)) {
  tc[i] = sum(chirurgia[,i]!='NaN') + sum(chirurgia[1:sum(chirurgia[,i]!='NaN'), i]=='NaN')
}
chirurgia = chirurgia[1:max(tc), ]

# vedere se qualcuno ha dati mancanti nel mezzo
for ( i in 1:dim(chirurgia)[2]) {
  print(sum(chirurgia[1:tc[i], i]=='NaN')) }
# in chirurgia il paziente 6 ha dati mancanti nel mezzo



# PROCESSING FISIOTERAPIA
df = NULL
nomi = colnames(fisioterapia)

for (i in 1:dim(fisioterapia)[2]) {
  if (fisioterapia[1, i] == 'Knee' && fisioterapia[4, i]=='X' && grepl('_01', nomi[i], fixed=T)) {
    df = cbind(df, fisioterapia[, i])
    colnames(df)[dim(df)[2]]=nomi[i]
  }
}

df = as.data.frame(df)

dati.fisio = as.data.frame(cbind(colnames(df), t(df)))

fisioterapia = as.data.frame(t(dati.fisio[, 6:dim(dati.fisio)[2]]))
qual.fisio = as.data.frame(t(dati.fisio[, 1:5])) # salvo le variabili qualitative che non si sa mai


# salvo il numero di osservazioni per ogni paziente di fisioterapia
tf = 1:dim(fisioterapia)[2]

for (i in 1:length(tf)) {
  tf[i] = sum(fisioterapia[,i]!='NaN') + sum(fisioterapia[1:sum(fisioterapia[,i]!='NaN'), i]=='NaN')
}
fisioterapia = fisioterapia[1:max(tf), ]

# vedere se qualcuno ha dati mancanti nel mezzo
for ( i in 1:dim(fisioterapia)[2]) {
  print(sum(fisioterapia[1:tf[i], i]=='NaN')) }
# in fisioterapia il paziente 27 ha dati mancanti nel mezzo



# trasformo tutti i dataframe in numerici
for (i in 1:dim(sani)[2])
sani[, i]= as.numeric(sani[,i])
for (i in 1:dim(chirurgia)[2])
chirurgia[, i]= as.numeric(chirurgia[,i])
for (i in 1:dim(fisioterapia)[2])
fisioterapia[, i]= as.numeric(fisioterapia[,i])

# pulizia di cose che non servono
rm(nomi, df, i, dati.chir, dati.fisio, dati.sani)



#### PLOT BASIC E INDICI DESCRITTIVI ####

# BASIC WARPING PER PLOT PRELIMINARI
minimo = min(c(min(sani, na.rm=T), min(chirurgia, na.rm=T), min(fisioterapia, na.rm=T)))
massimo = max(c(max(sani, na.rm=T), max(chirurgia, na.rm=T), max(fisioterapia, na.rm=T)))

x11()
par(mfrow=c(1,3))
plot( seq(0,1,length.out = ts[1]), sani[ 1:ts[1], 1], main='Control group', type='l', lwd=2, ylim=c(minimo, massimo), xlab='Time', ylab = 'Knee flexion angle')
for ( i in 1:dim(sani)[2]) {
  lines(seq(0,1,length.out = ts[i]), sani[ 1:ts[i], i], col=i)
}
plot( seq(0,1,length.out = tc[1]), chirurgia[ 1:tc[1], 1], main='Undergone surgery', type='l', lwd=2, ylim=c(minimo, massimo),  xlab='Time', ylab = 'Knee flexion angle')
for ( i in c(1:5, 7:dim(chirurgia)[2])) {
  lines(seq(0,1,length.out = tc[i]), chirurgia[ 1:tc[i], i], col=i)
}
plot( seq(0,1,length.out = tf[1]), fisioterapia[ 1:tf[1], 1], main='Undergone physioterapy', type='l', lwd=2, ylim=c(minimo, massimo),  xlab='Time', ylab = 'Knee flexion angle')
for ( i in c(1:26, 28:dim(fisioterapia)[2])) {
  lines(seq(0,1,length.out = tf[i]), fisioterapia[ 1:tf[i], i], col=i)
}

graphics.off()



# INDICI DEI VARI GRUPPI

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



#### BOXPLOT E FUNCTIONAL BOXPLOT ####

# SANI (primi 12 pazienti)
x11()
par(mfrow=c(2,6))
boxplot(sani[,1], col='darkgreen', ylim=c(min(sani[, 1:12], na.rm=T), max(sani[, 1:12], na.rm=T)))
for ( i in 1:11) {
boxplot(sani[,i], col='darkgreen', ylim=c(min(sani[, 1:12], na.rm=T), max(sani[, 1:12], na.rm=T))) }

# CHIRURGIA (primi 12 pazienti)
x11()
par(mfrow=c(2,6))
boxplot(chirurgia[,1], col='firebrick2', ylim=c(min(chirurgia[, 1:12], na.rm=T), max(chirurgia[, 1:12], na.rm=T)))
for ( i in 1:11) {
  boxplot(chirurgia[,i], col='firebrick2', ylim=c(min(chirurgia[, 1:12], na.rm=T), max(chirurgia[, 1:12], na.rm=T))) }


# FISIOTERAPIA (primi 12 pazienti)
x11()
par(mfrow=c(2,6))
boxplot(fisioterapia[,1], col='dodgerblue', ylim=c(min(fisioterapia[, 1:12], na.rm=T), max(fisioterapia[, 1:12], na.rm=T)))
for ( i in 1:11) {
  boxplot(fisioterapia[,i], col='dodgerblue', ylim=c(min(fisioterapia[, 1:12], na.rm=T), max(fisioterapia[, 1:12], na.rm=T))) }



library(fda)

x11()
fbplot(sani[1:min(ts), ])

x11()
fbplot(chir[1:min(tc), c(1:5, 7:dim(chirurgia)[2])])  # escludo paziente 6 perché fbplot non regge gli Na nel mezzo

x11()
fbplot(fisio[1:min(tf), c(1:26, 28:dim(fisioterapia)[2])]) # escludo paziente 28 perché fbplot non regge gli Na nel mezzo







