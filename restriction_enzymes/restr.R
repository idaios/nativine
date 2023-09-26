library(DECIPHER)
library(spgs)


data(RESTRICTION_ENZYMES)
ren = RESTRICTION_ENZYMES

tmp1 = gsub(pattern="/", "", x=RESTRICTION_ENZYMES)
tmp2 = grep(pattern="[^ACGT]", x=tmp1, perl=TRUE, invert=TRUE)

infofile = read.table("ddrad RE double digest.xlsx - RE characteristics.tsv", header=TRUE, sep="\t")
palindromic.indexes = infofile$Palindromic.enzymes == "Yes"

restr = infofile$RE[palindromic.indexes]
names(infofile)
restrsites = infofile[palindromic.indexes,"recognition.site..cut"]
restrsites = gsub(pattern="/", "", x=restrsites)
names(restrsites) = infofile$RE[palindromic.indexes]
restrsites

revComplement=reverseComplement(restrsites, case="as is")
restrSym = restrsites[revComplement == restrsites]

frequent = names(restrSym[sapply(restrSym, nchar) == 4])
rare = names(restrSym[sapply(restrSym, nchar) > 4])

frequent
rare

grep(restrSym["ApaI"], x=restrSym, ignore.case=TRUE)
restrSym[53]


file = "Vitis_vinifera.fa.gz"
seqs = readDNAStringSet(file)
seqs

resMat = matrix(0, nrow=length(frequent), ncol=length(rare))
colnames(resMat) = rare
rownames(resMat) = frequent
resMatConf = resMat
resMatLength =resMatMean = resMat
resMatDistr = matrix(list(), nrow=length(frequent), ncol=length(rare))
colnames(resMatDistr) = rare
rownames(resMatDistr) = frequent

dim(resMat)

## find the cases without a conflict
for(i in 1:length(frequent)){
    for(j in 1:length(rare)){
        resMatConf[i,j] = length(grep(restrSym[frequent[i]], restrSym[rare[j]])) == 0
    }
}

## analyze the results
for(chr in 1:length(seqs)){
    print(paste("CHR: ", chr))
    for(i in 1:nrow(resMat)){
        print(i)
        if ( is.na(ren[frequent[i]]) ){
            print( frequent[i] )
            next
        }
        for(j in 1:ncol(resMat)){
            if ( is.na(ren[rare[j]]) ){
                print(rare[j])
                next
            }
            print(j)
            digest1 = DigestDNA(ren[frequent[i]], seqs[chr], type="positions", strand="top")
            digest2 = DigestDNA(ren[rare[j]], seqs[chr], type="positions", strand="top")
            digest1Pos = digest1[[1]]$top
            digest2Pos = digest2[[1]]$top
            df = data.frame(x = c (rep(1, length(digest1Pos)), rep(2, length(digest2Pos) ) ), y = c(digest1Pos, digest2Pos))
            df = df[order(df[,2]),]
            pos1 = which(diff(df[,1], 1) != 0)
            pos2 = pos1+1
            lengths = df[pos2, 2] - df[pos1,2]
            resMat[i,j] = resMat[i,j] + length(lengths)
            resMatLength[i,j] = resMatLength[i,j] + sum(lengths)
            resMatMean[i,j] = resMatMean[i,j] + sum(lengths)
            resMatDistr[i,j][[1]] = c(resMatDistr[i,j][[1]], lengths)
        }
    }
}




for(i in 1:nrow(resMat)){
    print(i)
    if ( is.na(ren[frequent[i]]) ){
        print( frequent[i] )
        next
    }
    for(j in 1:ncol(resMat)){
        if ( is.na(ren[rare[j]]) ){
            print(rare[j])
            next
        }
        print(j)
        resMatMean[i,j] = resMatMean[i,j]/resMat[i,j]
    }
}


library(gplots)

resMatquantMin = resMat
resMatquantMax = resMat
resMatPerc = resMat
for(i in 1:nrow(resMat)){
    for(j in 1:ncol(resMat)){
        resMatquantMin[i,j] = quantile(resMatDistr[i,j][[1]], probs=0.05)
        resMatquantMax[i,j] = quantile(resMatDistr[i,j][[1]], probs=0.95)
        resMatPerc[i,j] = sum(resMatDistr[i,j][[1]] > 200 & resMatDistr[i,j][[1]] < 600)
    }
}




mypallette = colorRampPalette(c("blue", "white", "red"))(n=256)
mypallette
pdf("qualitativeFeaturesChr1VitisAll.pdf")
heatmap.2(x=resMatPerc, Rowv=F, Colv=F, dendrogram='none', trace='none', cexRow=1, main="# > 200 & < 600")
heatmap.2(x=resMat, Rowv=F, Colv=F, dendrogram='none', trace='none', cexRow=1, main="# fragments")
heatmap.2(x=resMatLength, Rowv=F, Colv=F, dendrogram='none', trace='none', cexRow=1, main="Total Length Fragments")
heatmap.2(x=log10(resMatMean/400), Rowv=F, Colv=F, dendrogram='none', trace='none', cexRow=1, col=mypallette, main="mean length dev from 400 (white is 400)")
dev.off()

save.image("back.RData")


