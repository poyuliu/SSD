# Taxa lineage contribution analysis
# http://www.phytools.org/eqg/Exercise_3.2/

PCn <- cmdscale(uk$distMw,k=(ncol(ktu)-1),eig = T)
PCn$eig <- PCn$eig[PCn$eig>0]
prop.table(PCn$eig)*100
kmer <- kluster$kmer.table[,match(rownames(ktu),colnames(kluster$kmer.table))]
ktree <- ape::as.phylo(hclust(as.dist(1-coop::cosine(kmer))))
subktrees <-ape::subtrees(ktree)

subktips <- list()
tipsums <- list()
for(i in 1:length(subktrees)){
  subktips[[i]] <- subktrees[[i]]$tip.label
  tipsums[[i]] <- colSums(ktu[match(subktips[[i]],rownames(ktu)),])
}

tipsums <- do.call(rbind,tipsums)


ori <- as.matrix(dist(rbind(O=0,PCn$points)))
ori <- ori[,1][-1]
#all positive eig
r.ori <- c()
p.ori <- c()
for(j in 1:nrow(tipsums)){
  tests.ori <- cor.test(ori,log1p(t(tipsums))[,j])
  r.ori[j] <- tests.ori$estimate
  p.ori[j] <- tests.ori$p.value
}
pa.ori <- p.adjust(p.ori,"fdr")

#PCo decomposed
r.pc2 <- c()
p.pc2 <- c()
for(j in 1:nrow(tipsums)){
  tests.pc2 <- cor.test(PCn$points[,2],log1p(t(tipsums))[,j])
  r.pc2[j] <- tests.pc2$estimate
  p.pc2[j] <- tests.pc2$p.value
}
pa.pc2 <- p.adjust(p.pc2,"fdr")

ssdiv <- SsDi(ktu,kluster$kmer.table)
#alpha ssdiv
r.alpha <- c()
p.alpha <- c()
for(j in 1:nrow(tipsums)){
  tests.alpha <- cor.test(ssdiv,log1p(t(tipsums))[,j])
  r.alpha[j] <- tests.alpha$estimate
  p.alpha[j] <- tests.alpha$p.value
}
pa.alpha <- p.adjust(p.alpha,"fdr")
R2.alpha <- (r.alpha[which(pa.alpha < 0.001)])^2

which(pa.ori < 0.001)
kaxa[match(subktips[[1009]],rownames(kaxa)),]

ori.sig <- tipsums[which(pa.ori < 0.001),]
rownames(ori.sig) <- which(pa.ori < 0.001)
plot.ordi(pcoa,groups,colset = colset.d.2,ordiellipse = T,ordispider = F,legend = levels(groups),test = T)
fits <- vegan::envfit(pcoa$pcoa,t(ori.sig))
plot(fits,p.max = 0.001)
