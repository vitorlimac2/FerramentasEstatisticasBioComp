## Autor: Vitor Lima Coelho
## Aula: Ferramentas estatísticas para Biologia Computacional e Bioinformatica
## Data: 11/03/2021

## importacoes

library(caret)

#################################################################
## Normalizacao

#### dados de exemplo

b1 = c(56,75,45,71,62,64,58,80,76,61)
b2 = c(66,70,40,60,65,56,59,77,67,63)
genes = c("Gene1",
          "Gene2",
          "Gene3",
          "Gene4",
          "Gene5",
          "Gene6",
          "Gene7",
          "Gene8",
          "Gene9",
          "Gene10")

df = data.frame(row.names = genes, b1,b2)

df

rownames(df) <- genes

## Padronizacao


padronizacao <- as.data.frame(scale(df))

## escalonamento

preproc2 <- preProcess(df, method=c("range"))
escalonamento <- predict(preproc2, df)

## Transformacao log

norm_log = log(df)

##################################
### gráfico boxplot


boxplot(df, 
        horizontal = TRUE, 
        axes = TRUE, 
        staplewex = 1, 
        las=1, 
        col=c("cornflowerblue", "deepskyblue4"),
        xlab= "Expressão")
text(x=fivenum(df$b1), labels=fivenum(df$b1), y=1.5, col="cornflowerblue")
text(x=fivenum(df$b2), labels=fivenum(df$b2), y=2.5, col="deepskyblue4")




##############################################
### sintese e exploracao dos dados


mean(df$b1)
mean(df$b2)

median(df$b1)
median(df$b2)

min(df$b1)
min(df$b2)

max(df$b1)
max(df$b2)

quantile(df$b1)
quantile(df$b2)


IQR(df$b1)
IQR(df$b2)
################################################
##### Agrupamento

library(NOISeq)

mycounts <- read.table("C:/Users/vitor/Downloads/quantificacao.tsv", row.names = 1, header = TRUE)

head(mycounts)

myfactors <- data.frame(Stage = c("Epi", "Epi","Epi","Trypo","Trypo","Trypo"),
                        StageRun = c("Epi_1", "Epi_2","Epi_3","Trypo_1","Trypo_2","Trypo_3"),
                        Run = c(rep("R1", 6)))
myfactors

mydata <- readData(data = mycounts, factors = myfactors)
par(mfrow = c(1, 1))
myPCA = dat(mydata, type="PCA")
explo.plot(myPCA, factor = "Stage")

mydata_corr1 = ARSyNseq(mydata, factor="Stage", batch = FALSE, norm = "tmm",logtransf = FALSE)
myPCA_corr1 = dat(mydata_corr1, type="PCA")
explo.plot(myPCA_corr1, factor = "Stage")

