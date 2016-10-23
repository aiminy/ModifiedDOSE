##' interal method for enrichment analysis
##'
##' using the hypergeometric model
##' @title enrich.internal
##' @param gene a vector of entrez gene id.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param qvalueCutoff cutoff of qvalue
##' @param USER_DATA ontology information
##' @param otherGO
##' @param whichway
##' @return  A \code{enrichResult} instance.
##' @importClassesFrom methods data.frame
##' @importFrom qvalue qvalue
##' @importFrom methods new
##' @importFrom stats phyper
##' @importFrom stats p.adjust
##' @keywords manip
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
enricher_internal <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.5,
                              USER_DATA,otherGO,whichway){
  
  
    ##Use otherGO
    cat("Use otherGO\n")
    #print(head(otherGO[[1]]$GO[,1:7]))
    
    pNew<-ReCalculatePvalue(otherGO)
    
    #print(pNew)
    
    cat("done\n")
    
    ## query external ID to Term ID
    gene <- as.character(unique(pNew$gene))
    
    cat("How many genes: ",length(gene),"\n")
    
    #Based pwf DE gene to get all mapped GO terms 
    qExtID2TermID <- EXTID2TERMID(gene, USER_DATA)
    
    #qExtID2TermID<-pNew$gene2cat
    
    #print(qExtID2TermID)
    
    cat("qExtID2TermID\n")
    cat(length(qExtID2TermID),"\n")
    cat("qExtID2TermID_done\n")
    
    qTermID <- unlist(qExtID2TermID)
    if (is.null(qTermID)) {
        message("No gene can be mapped....")
        message("--> return NULL...")
        return(NULL)
    }

    #x=qTermID
    # ## Term ID -- query external ID association list.
    
    #extID TERMID
    #gene1 GO1
    #gene1 GO2
    #gene1 GO3
    
    qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                             times=lapply(qExtID2TermID, length)),
                                   termID=qTermID)

    
    #x=list(qExtID2TermID=qExtID2TermID,qTermID=qTermID)
    
     qExtID2TermID.df <- unique(qExtID2TermID.df)
     
     qTermID2ExtID <- with(qExtID2TermID.df,
                           split(as.character(extID), as.character(termID)))
    
     
     #x=  qTermID2ExtID
     #get background genes
     
     extID <- ALLEXTID(USER_DATA)
     
      if(!missing(universe)) {
          extID <- intersect(extID, universe)
      }
      
     #extID<-unique(unlist(qTermID2ExtID))
     
     cat(dim(otherGO[[2]])[1],"\n")
     
     matchextID<-intersect(extID,rownames(otherGO[[2]]))
     matchPwf<-otherGO[[2]][which(rownames(otherGO[[2]]) %in% matchextID),]
     
     cat(dim(matchPwf)[1],"\n")
     
     extID<-rownames(matchPwf)
     
     qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
    
     
     
     ## Term ID annotate query external ID
     #GO termID
      qTermID <- unique(names(qTermID2ExtID))
      
    #To make a list with GO terms as names,and genes as the elemnets of this list
      #For example:
      #$`GO:0000408`
      #[1] "112858"
      
     termID2ExtID <- TERMID2EXTID(qTermID, USER_DATA)
     
     #termID2ExtID<-pNew$cat2gene
     
     termID2ExtID <- lapply(termID2ExtID, intersect, extID)
    
    
     
     geneSets <- termID2ExtID
     
     #To filters GOs based on size of GO
     idx <- get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
     if (sum(idx) == 0) {
         msg <- paste("No gene set have size >", minGSSize, "...") 
         message(msg)
         message("--> return NULL...")
         return (NULL)
     }
     
     termID2ExtID <- termID2ExtID[idx]
     qTermID2ExtID <- qTermID2ExtID[idx]
     qTermID <- unique(names(qTermID2ExtID))
  
    
     
     
     
     PfromGoseq<-ReCalculatePvalue2(termID2ExtID,matchPwf)
     
     
    
     #x=list(qExtID2TermID=qExtID2TermID,qExtID2TermID.df=qExtID2TermID.df,qTermID2ExtID=qTermID2ExtID,termID2ExtID=termID2ExtID,
     #       PfromGoseq=PfromGoseq)
     
    
     
     pvaluesFromGoSeq<-PfromGoseq$pvals[,2]
     

     
       

    # ## prepare parameter for hypergeometric test

     # get number of genes for each GO term
     k <- sapply(qTermID2ExtID, length)
     k <- k[qTermID]
     M <- sapply(termID2ExtID, length)
     M <- M[qTermID]

     N <- rep(length(extID), length(M))
    # ## n <- rep(length(gene), length(M)) ## those genes that have no annotation should drop.
     n <- rep(length(qExtID2TermID), length(M))

     args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                           numW=M,        ## White balls
                           numB=N-M,      ## Black balls
                           numDrawn=n)    ## balls drawn


     #print(args.df)

     if(whichway == "HT") {
     ## calcute pvalues based on hypergeometric model
      pvalues <- apply(args.df, 1, function(n)
                       phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
                       )
    }else if(whichway == "WHT")
      {pvalues<-pvaluesFromGoSeq}

    # ## gene ratio and background ratio

    # #print(data.frame(a=k, b=n))
    # k: DE genes in this GO
    # n: all DE genes
     
    GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    #
    # #print(data.frame(a=M, b=N))
    #
    # M: all gens
    # N:
    
    BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )
    #
    #
    Over <- data.frame(ID = as.character(qTermID),
                       GeneRatio = GeneRatio,
                       BgRatio = BgRatio,
                       pvalue = pvalues,
                       stringsAsFactors = FALSE)
    
    
  
    #print(head(Over))  
    #print(head(pNew$pvals))
    
    
    
    # # # 
    p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }
    # # 
    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       p.adjust = p.adj,
                       qvalue = qvalues,
                       geneID = geneID,
                       Count = k,
                       stringsAsFactors = FALSE)
    # 
     Description <- TERM2NAME(qTermID, USER_DATA)
    # #  
     if (length(qTermID) != length(Description)) {
         idx <- qTermID %in% names(Description)
         Over <- Over[idx,]
     }
     Over$Description <- Description
     nc <- ncol(Over)
     Over <- Over[, c(1,nc, 2:(nc-1))]
    # 
    # 
     Over <- Over[order(pvalues),]
    #
     Over <- Over[ Over$pvalue <= pvalueCutoff, ]
     Over <- Over[ Over$p.adjust <= pvalueCutoff, ]

     if (! any(is.na(Over$qvalue))) {
         Over <- Over[ Over$qvalue <= qvalueCutoff, ]
     }

     Over$ID <- as.character(Over$ID)
     Over$Description <- as.character(Over$Description)
    # 
     row.names(Over) <- as.character(Over$ID)
    # #  
    # # 
    # #  
    # # # #print(head(Over))
    # #  
    x <- new("enrichResult",
             result         = Over,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             gene           = as.character(gene),
             universe       = extID,
             geneSets       = geneSets,
             organism       = "UNKNOWN",
             keytype        = "UNKNOWN",
             ontology       = "UNKNOWN",
             readable       = FALSE
             )
    return (x)
}

EXTID2TERMID <- function(gene, USER_DATA) {
    EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)

    #cat("OK\n")
    
    
    #like gene2cat
    
    #$ENSG00000000971
    #[1] "GO:0001948" "GO:0002252" "GO:0002253" "GO:0002376" "GO:0002526" "GO:0002673" "GO:0002682" "GO:0002684" "GO:0002697" "GO:0002920" "GO:0003674" "GO:0005488" "GO:0005515" "GO:0005539" "GO:0005575"
    #[16] "GO:0005576" "GO:0005615" "GO:0006508" "GO:0006950" "GO:0006952" "GO:0006954" "GO:0006955" "GO:0006956" "GO:0006957" "GO:0006959" "GO:0008150" "GO:0008152" "GO:0008201" "GO:0009605" "GO:0009611"
    #[31] "GO:0010467" "GO:0010468" "GO:0016485" "GO:0019222" "GO:0019538" "GO:0030162" "GO:0030449" "GO:0031347" "GO:0031982" "GO:0031988" "GO:0032101" "GO:0043167" "GO:0043168" "GO:0043170" "GO:0043226"
    #[46] "GO:0043227" "GO:0043230" "GO:0043394" "GO:0043395" "GO:0044238" "GO:0044421" "GO:0044699" "GO:0044710" "GO:0045087" "GO:0048518" "GO:0048583" "GO:0048584" "GO:0050727" "GO:0050776" "GO:0050778"
    #[61] "GO:0050789" "GO:0050896" "GO:0051246" "GO:0051604" "GO:0060255" "GO:0065007" "GO:0065010" "GO:0070062" "GO:0070613" "GO:0071704" "GO:0072376" "GO:0072562" "GO:0080090" "GO:0080134" "GO:0097367"
    #[76] "GO:1901681" "GO:1903034" "GO:1903317" "GO:1903561" "GO:2000257"
    
    #print(EXTID2PATHID)
    #cat(length(EXTID2PATHID),"\n")
    #cat("OKdone\n")
    
    #remove genes that no GO terms are mapped
    
    qExtID2Path <- EXTID2PATHID[gene]
    len <- sapply(qExtID2Path, length)
    notZero.idx <- len != 0
    qExtID2Path <- qExtID2Path[notZero.idx]

    
    return(qExtID2Path)
}

ALLEXTID <- function(USER_DATA) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- unique(unlist(PATHID2EXTID))
    return(res)
}


TERMID2EXTID <- function(term, USER_DATA) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- PATHID2EXTID[term]
    return(res)
}

TERM2NAME <- function(term, USER_DATA) {
    PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
    if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
        return(as.character(term))
    }
    return(PATHID2NAME[term])
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
    if (is.na(minGSSize) || is.null(minGSSize))
        minGSSize <- 0
    if (is.na(maxGSSize) || is.null(maxGSSize))
        maxGSSize <- .Machine$integer.max

    ## index of geneSets in used.
    ## logical
    idx <- sapply(geneSets, length) > minGSSize & sapply(geneSets, length) < maxGSSize
    return(idx)
}

ReCalculatePvalue <- function(Example.Go.adjusted.by.exon) {
  
  pwf=Example.Go.adjusted.by.exon[[2]]
  pwf=unfactor(pwf)
 
  cat(dim(pwf)[1],"\n")
  
  test.cats=c("GO:CC","GO:BP","GO:MF")
  gene2cat=getgo(rownames(pwf),'hg19','ensGene',fetch.cats=test.cats)
  names(gene2cat)=rownames(pwf)
  
  cat2gene=reversemapping(gene2cat)
  gene2cat=reversemapping(cat2gene)
  nafrac=(sum(is.na(pwf$pwf))/nrow(pwf))*100
  
  if(nafrac>50){
    warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
  }
  pwf$pwf[is.na(pwf$pwf)]=pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)],pwf$bias.data)]
  
  unknown_go_terms=nrow(pwf)-length(gene2cat)
  
  cat(dim(pwf)[1],"\n")
  
  if(unknown_go_terms>0 ){
    message(paste("For",unknown_go_terms,"genes, we could not find any categories. These genes will be excluded."))
    message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf=pwf[rownames(pwf) %in% names(gene2cat),]
  } 
  
  cat(dim(pwf)[1],"\n")
  
  #A few variables are always useful so calculate them
  cats=names(cat2gene)
  DE=rownames(pwf)[pwf$DEgenes==1]
  
  #total number of DE genes
  num_de=length(DE)
  
  #total number of genes
  num_genes=nrow(pwf)
  
  
  pvals=data.frame(category=cats,over_represented_pvalue=NA,under_represented_pvalue=NA,stringsAsFactors=FALSE,numDEInCat=NA,numInCat=NA)
  
  degenesnum=which(pwf$DEgenes==1)
  #Turn all genes into a reference to the pwf object
  cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
  #This value is used in every calculation, by storing it we need only calculate it once
  alpha=sum(pwf$pwf)
  
  pvals[,2:3]=t(sapply(cat2genenum,function(u){
    
    #print(u)
    
    #The number of DE genes in this category
    num_de_incat=sum(degenesnum%in%u)
    #The total number of genes in this category
    num_incat=length(u)
    #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
    avg_weight=mean(pwf$pwf[u])
    weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
    if(num_incat==num_genes){ weight=1 } #case for the root GO terms			
    #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
    
    
    #cat(num_de_incat,"\t",num_incat,"\t",num_genes,"\t",num_incat,"\t",num_de,"\t",weight,"\n")
    
    
    c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
      +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
      pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
  }))
  
  
  
  re<-list(gene=DE,pvals=pvals,gene2cat=gene2cat,cat2gene=cat2gene)
  
  return(re)
  
}


ReCalculatePvalue2 <- function(cat2gene,pwf) {
  
  pwf=pwf
  pwf=unfactor(pwf)
  
  cat(dim(pwf)[1],"\n")
  
 # test.cats=c("GO:CC","GO:BP","GO:MF")
#  gene2cat=getgo(rownames(pwf),'hg19','ensGene',fetch.cats=test.cats)
#  names(gene2cat)=rownames(pwf)
  
  gene2cat=reversemapping(cat2gene)
  cat2gene=reversemapping(gene2cat)
 
  nafrac=(sum(is.na(pwf$pwf))/nrow(pwf))*100
  
  if(nafrac>50){
    warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
  }
  
  pwf$pwf[is.na(pwf$pwf)]=pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)],pwf$bias.data)]
  
  unknown_go_terms=nrow(pwf)-length(gene2cat)
  
  cat(dim(pwf)[1],"\n")
  
  if(unknown_go_terms>0 ){
    message(paste("For",unknown_go_terms,"genes, we could not find any categories. These genes will be excluded."))
    message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf=pwf[rownames(pwf) %in% names(gene2cat),]
  } 
  
  cat(dim(pwf)[1],"\n")
  
  #A few variables are always useful so calculate them
  cats=names(cat2gene)
  DE=rownames(pwf)[pwf$DEgenes==1]
  
  #total number of DE genes
  num_de=length(DE)
  
  #total number of genes
  num_genes=nrow(pwf)
  
  
  pvals=data.frame(category=cats,over_represented_pvalue=NA,under_represented_pvalue=NA,stringsAsFactors=FALSE,numDEInCat=NA,numInCat=NA)
  
  degenesnum=which(pwf$DEgenes==1)
  #Turn all genes into a reference to the pwf object
  cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
  #This value is used in every calculation, by storing it we need only calculate it once
  alpha=sum(pwf$pwf)
  
  pvals[,2:3]=t(sapply(cat2genenum,function(u){
    
    #print(u)
    
    #The number of DE genes in this category
    num_de_incat=sum(degenesnum%in%u)
    #The total number of genes in this category
    num_incat=length(u)
    #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
    avg_weight=mean(pwf$pwf[u])
    weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
    if(num_incat==num_genes){ weight=1 } #case for the root GO terms			
    #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
    
    
    #cat(num_de_incat,"\t",num_incat,"\t",num_genes,"\t",num_incat,"\t",num_de,"\t",weight,"\n")
    
    
    c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
      +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
      pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
  }))
  
  
  
  re<-list(gene=DE,pvals=pvals,gene2cat=gene2cat,cat2gene=cat2gene)
  
  return(re)
  
}