# simulation of data 
#
# 
# TABLE 1
names = expand.grid(sample(LETTERS,10),sample(letters,10),sample(LETTERS,10))
names = apply(names,1,paste,collapse="" )
names = sample(names,200)

#50 lineas serán con la misma distribución normal
#5 lineas con una distribución exponencial 
#5 lineas con una distribución 

matrixids<-1:200
norm_list1=sample(matrixids,130) #equal norm
matrixids<-matrixids[!(matrixids %in% norm_list1)]
norm_list2=sample(matrixids, 10) #dif norm
matrixids<-matrixids[!(matrixids %in% norm_list2)]
nonorm_list1=sample(matrixids,50) #equal no_norm
matrixids<-matrixids[!(matrixids %in% nonorm_list1)]
nonorm_list2=matrixids #last 10 numbers not iequal no normal
rm(matrixids)

mus=sample(1:1000,size = 130)
sds=runif(1)*sample(1:100,130,replace=T)

Matriz1<-matrix(data = 0,nrow = 200,ncol = 1000)
Matriz2<-matrix(data = 0,nrow = 200,ncol = 1000)

Matriz1[norm_list1,] = do.call(rbind,lapply(seq_len(130),FUN= function(x) rnorm(n=1000, mean=mus[x], sd=sds[x] ) ))
Matriz2[norm_list1,] = do.call(rbind,lapply(seq_len(130),FUN= function(x) rnorm(n=1000, mean=mus[x], sd=sds[x] ) ))

Matriz1[norm_list2,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rnorm(n=1000, mean=sample(1:1000,size = 1), sd=runif(1)*sample(1:100,1) ) ))
Matriz2[norm_list2,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rnorm(n=1000, mean=sample(1:1000,size = 1), sd=runif(1)*sample(1:100,1) ) ))


shpe=runif(10,1,100)

Matriz1[nonorm_list1,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rgamma(1000, shape=shpe )*runif(n = 1000,min=1, max=100)))
Matriz2[nonorm_list1,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rgamma(1000, shape=shpe )*runif(n = 1000,min=1, max=100)))

Matriz1[nonorm_list2,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rgamma(1000, shape=runif(1,1,100) )*runif(n = 1000,min=1, max=100)))
Matriz2[nonorm_list2,] = do.call(rbind,lapply(seq_len(10),FUN= function(x) rgamma(1000, shape=runif(1,1,100) )*runif(n = 1000,min=1, max=100)))


save(Matriz1,Matriz2, file = "sampler_tables.RData")
