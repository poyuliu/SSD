SsDD <- function(KTU,kmer,cores=1){
  require(foreach)
  kmer <- kmer[,match(rownames(KTU),colnames(kmer))]
  #distM <- matrix(data = NA,nrow = ncol(KTU),ncol = ncol(KTU))
  cos <- 1-coop::cosine(kmer)
  
  #for(i in 1:ncol(KTU)){
  #  for(j in 1:ncol(KTU)){
  #    us <- c(1:nrow(KTU))[which((KTU[,i]>0)+(KTU[,j]>0)==1)]
  #    #ss <- intersect(which(KTU[,i]>0),which(KTU[,j]>0))
  #    #distM[i,j] <- sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[ss,ss][lower.tri(cos[ss,ss])])
  #    distM[i,j] <- sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[lower.tri(cos)])
  #  }
  #}
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  distM <- foreach::foreach(i = 1:ncol(KTU), .combine = rbind, .inorder = T) %:% 
    foreach::foreach(j = 1:ncol(KTU), .combine = rbind, .inorder = T) %dopar%{
      us <- c(1:nrow(KTU))[which((KTU[,i]>0)+(KTU[,j]>0)==1)]
      sum(cos[us,us][lower.tri(cos[us,us])])/sum(cos[lower.tri(cos)])
    }
  distM <- Matrix::sparseMatrix(x = distM[,1],i=rep(1:ncol(KTU),each=ncol(KTU)),j=rep(1:ncol(KTU),ncol(KTU)))
  parallel::stopCluster(cl)  
  colnames(distM) <- rownames(distM) <- colnames(KTU)
  
  bc.dist <- as.matrix(vegan::vegdist(t(KTU)))
  
  distMu <- as.dist(distM)
  distMw <- as.dist(distM * bc.dist)
  return(list(distMu=distMu,distMw=distMw))
}


SsDi <- function(KTU,kmer){
  kmer <- kmer[,match(rownames(KTU),colnames(kmer))]
  ssd <- c()
  for(i in 1:ncol(KTU)){
    kmer1 <- kmer[,match(rownames(KTU)[which(KTU[,i]>0)],colnames(kmer))]
    cos1 <- (1-coop::cosine(kmer1))
    diag(cos1) <- NA
    ssd[i] <- mean(colMeans(cos1,na.rm = T)*(sum(KTU[,i]>0)))
  }
  names(ssd) <- names(KTU)
  return(ssd)
}

tetra.freq2 <- function (repseq, pscore = FALSE, file = TRUE, cores = 1) 
{
  require(foreach, quietly = T)
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  rev_comp <- function(x) {
    x.r <- paste0(rev(strsplit(x, "")[[1]]), collapse = "")
    x.r <- gsub("A", "1", x.r)
    x.r <- gsub("T", "2", x.r)
    x.r <- gsub("C", "3", x.r)
    x.r <- gsub("G", "4", x.r)
    x.r <- gsub("1", "T", x.r)
    x.r <- gsub("2", "A", x.r)
    x.r <- gsub("3", "G", x.r)
    x.r <- gsub("4", "C", x.r)
    return(x.r)
  }
  readseq <- function(seq) {
    x <- seq
    x <- paste0(x, collapse = "")
    x.r <- rev_comp(x)
    return(list(x, x.r))
  }
  if (isTRUE(file)) {
    scaffolds = Biostrings::readDNAStringSet(filepath = repseq, 
                                             use.names = T)
  }
  else scaffolds = Biostrings::DNAStringSet(repseq)
  asv.id <- scaffolds@ranges@NAMES
  species <- lapply(scaffolds, function(x) readseq(as.character(x)))
  DNA <- c("A", "T", "C", "G")
  tetra.mer <- expand.grid(DNA, DNA, DNA, DNA)
  tetra.mer <- do.call(paste0, tetra.mer)
  tetra.table <- matrix(nrow = 256, ncol = length(species))
  rownames(tetra.table) <- tetra.mer
  if (isTRUE(pscore)) {
    for (i in 1:256) {
      for (j in 1:length(species)) {
        single.forward <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[1]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[1]]), 0)
        single.reverse <- ifelse(length(grep(tetra.mer[i], 
                                             species[[j]][[2]])) > 0, grep(tetra.mer[i], 
                                                                           species[[j]][[2]]), 0)
        tetra.table[i, j] <- (single.forward + single.reverse)
      }
    }
  }
  else if (!isTRUE(pscore)) {
    tetra.table <- foreach(j = 1:length(species), .combine = cbind) %dopar% 
      {
        single.forward <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
                                                                         Biostrings::DNAString(unlist(species[[j]][1]))))
        single.reverse <- S4Vectors::elementNROWS(Biostrings::matchPDict(Biostrings::PDict(tetra.mer), 
                                                                         Biostrings::DNAString(unlist(species[[j]][2]))))
        tetras <- (single.forward + single.reverse)
        tetras
      }
    rownames(tetra.table) <- tetra.mer
  }
  colnames(tetra.table) <- asv.id
  tetra.table <- prop.table(tetra.table, 2)
  return(tetra.table)
}

ssd.rarecurve <- function(ktu,kmer,steps=20,plot=TRUE){
  ssdi <- list()
  for(i in 1:ncol(ktu)){
    subsampsize <- round(seq(0,colSums(ktu)[i],length.out=(steps+1)))
    subssdi <- c()
    for(j in 1:c(steps+1)){
      subdataset <- as.data.frame(table(sample(rep.int(rownames(ktu),ktu[,i]),size = subsampsize[j])),row.names = 1)
      subssdi[j] <- ifelse(length(subdataset)==0,0,SsDi(subdataset,kmer))
    }
    names(subssdi) <- subsampsize
    ssdi[[i]] <- subssdi
  }
  if(isTRUE(plot)){
    plot(seq(1,max(as.numeric(names(unlist(ssdi)))),length.out=steps),seq(0,max(unlist(ssdi)),length.out=steps),
         xlab="Seq_depth",ylab="SSD",las=1,type="n")
    for(i in 1:ncol(ktu)){
      lo <- loess(ssdi[[i]]~as.numeric(names(ssdi[[i]])))
      lines(as.numeric(names(ssdi[[i]])),predict(lo))
    } 
  }
  invisible(ssdi)
}
