dados = read.table(file.choose(), header =TRUE)
dados = read.csv(file.choose(), header= TRUE, sep=",")



# Comandos para instalação dos pacotes
install.packages(sp, dependences=T)
install.packages(gstat, dependences =T)
install.packages(geoR, dependences = T)

# Comandos para carregar os pacotes 
library(sp)
library(geoR)
library(gstat)

data(walker)            #Baixa o conjunto de dados do Walker Lake
head(walker)            # Observa os primeiros valores do Walker Lake
tail(walker)            # Observa os últimos valores do Walker Lake
summary(walker)         # Realiza um sumário estatístico do Walker Lake


for (i in 1:10){
  print(i)
}

i = 1
while(i <= 10){
  print(i)
  i = i + 1
}

?plot
plot(walker, cex=TRUE )

head(walker)

# HISTOGRAMS 
?par
dados2 = dados[dados$Au > 0,]
par(mfrow = c(1,2))
hist(walker$U, main= "histograma da variável U ", breaks =15, xlab= "U", ylab= "Frequência")
hist(walker$V, main= "histograma da variável V", breaks = 15,  xlab= "V", ylab = "Frequência")
dev.off()

par(mfrow = c(1,2))
boxplot(walker$U, main= "Boxplot da variável U ", breaks =15, ylab= "U")
boxplot(walker$V, main= "Boxplot da variável V", breaks = 15,  ylab= "V")

#Disperion 
linear = lm(walker$U~walker$V)
summary(linear)
plot(walker$V, walker$U, xlab="V", ylab="U")
abline(linear)

Walker_out = walker[!is.na(walker$U),]
walker_out =walker_out[walker_out$U < 5000,]

# PREPARING GRID 

# ASSOCIAR AS COORDENADAS NO ESPAÇO
coordinates(walker) = c("X", "Y")

# SSPLOT DA VARIÁVEL v
spplot(walker, c("V"),scales = list(draw =T))

# ssPLOT DA VARIÁVEL V E U 
spplot(walker, c("V","U"),scales = list(draw =T))


var(dados$Cu)
points(as.geodata(walker$V))


# Criar um grid 
grid_stat = makegrid(walker,cellsize = 5)
grid_stat = SpatialPixels(SpatialPoints(grid_stat))


# Calcular vizinho mais próximo 

ca = raster(walker, res= 5)
gs = gstat(NULL, "V", V~1, walker, nmax=1)
nn = interpolate(ca, gs)
print(nn)
summary(nn$V.pred)
cellStats(nn$V.pred,stat='mean')
plot(nn, axes=T)


#CALCULATE IDW
idw_result = idw(V~1, walker, grid_stat, idp =1.0)
summary(idw_result)
spplot(idw_result["var1.pred"])

#CACULATE VARIOGRAMS
hscat(V ~1 , walker, (0:9)*10)


# Criar um grid 
grid_stat = makegrid(walker,cellsize = 5)
grid_stat = SpatialPixels(SpatialPoints(grid_stat))



#Variância da variável V
var(walker$V)

#Variogramas experimentais

v.dir = variogram(V~1, walker, width=10, cutoff= 200, tol.hor=45, alpha = (0:7)*22.5 )

# Modelagem da primeira estrutura 

v.anis1 = vgm(59929,"Sph",40,20000, anis=c(157,0.6))

# Modelagem da segunda estrutura 

v.anis2 = vgm(20000,"Exp",100, 0, anis=c(157,0.6),add.to =v.anis1 )

# Plotagem do gráfico
plot(v.dir, v.anis2)


# CROSS VALIDATION KRIGING
?krige.cv
cv = krige.cv(V~1, walker, v.anis2 , nmin= 3, nmax=10, maxdist=100, nfold=20)
summary(cv)
bubble(cv["residual"])


print(v.dir)
plot(v.dir)

vgm()



library(geoR)
library(gstat)

# Baixar os dados do Walker Lake
data(walker)

# Associar coordenadas ao banco de dados 
coordinates(walker) = c('X', 'Y')


# H-scatterplot da variável V
hscat(V~1, walker, (0:9)*5)

#Criar variograma experimental utilizando a biblioteca geoR
var = variog(as.geodata(walker["V"]), uvec = seq(from=0,to=500,by=10),max.dist = 500, direction= 157.5*pi/180, tolerance = pi/4)

#Fitar manualmente o variograma  
ve.eye = eyefit(var)

#Transformar o modelo fitado em um modelo do gstat
ve.fit = as.vgm.variomodel(v.eye[[1]])



v.anis = vgm(69929,"Exp",20,20000, anis=c(157,0.6))
plot(v.dir, v.anis)



v.dir = variogram(V~1, walker, width=10, cutoff= 200, tol.hor=45, alpha = (0:7)*22.5, covariogram = T )
v.anis = vgm(69929,"Exp",20,20000, anis=c(157,0.6))
plot(v.dir, v.anis)




# Criar um grid 
grid_stat = makegrid(walker,cellsize = 5)
grid_stat = SpatialPixels(SpatialPoints(grid_stat))



# KRIGING VALUES 
kriged = krige(V~1, walker, grid_stat,nmin=1, nmax=3, maxdist=100, v.anis2)
spplot(kriged['var1.pred'], scales = list(draw =T))
spplot(kriged['var1.var'] , scales = list(draw =T))
summary(kriged)



# PERFORM COKRIGING
dados = dados[dados$Au >= 0,]
g = gstat(NULL, "cu", Cu~1, dados)
g = gstat(g , "Au", Au~1, dados)
print(g)
v.dir = variogram(g)
vm.fit = fit.lmc(v.dir, g, vgm(7,"Sph",30,1, anis=c(157.5,0.5)))
plot(v.dir, vm.fit)

kriged2 = predict(vm.fit, grid_stat)
spplot(kriged2["Au.pred"])
