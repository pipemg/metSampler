setwd("~/Dropbox/Compu_Proteinas/Doctorado/Cursos/TIB2018-Bioconductor/metabolome/metaboloma/metSampler")

dim(Matriz1)
dim(Matriz2)

p <- ks.test(Matriz1, Matriz2)
class(p)
p$p.value

matrix <- c(Matriz1, Matriz2)


for (i in 1:2) {
  p <- ks.test(Matriz1[i,], Matriz2[i,])
  p[i]$p.value
}

length(p[3])

result <- ks.test(Matriz1[1,], Matriz2[1,])
result$p

test <- function(x,y){
  result[i] <- ks.test(Matriz1[i,], Matriz2[i,])
  result[i]$p.value
}

apply(matrix, MARGIN = 1, test)
