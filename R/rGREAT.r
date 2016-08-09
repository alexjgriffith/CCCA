#!/usr/bin/env R
#
# This file is part of mulcal,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com


#' Gene Ranges
#'
#' Implements the two step Stanford GREAT genomic region generation.
#' First each gene is assigned a basal region defined as the TSS +/-
#' the proxUp and proxDown value (dependent on strand). Second each
#' Gene has its effective region expaned up to the lower of a distal
#' maximum and the basal limits of other genes.
#' The returned values are in the order of the input chrom tss strand
#' variables.
#' note: the indicies of the chrom tss and strand parameters should
#' all relate to the same gene, ie a chrom[100] tss[100] and strand[100]
#' should all be infromation from the same gene
#' @param chrom A list of chromosomes, one for each gene
#' @param tss The list of gene TSS
#' @param strand The stream directions of each of the genes
#' @param proxUp Upstream basal extension (independent of other genes)
#' @param proxDown Downstream basal extension (independent of other genes)
#' @param distal The limit of gene extension (in either direction)
#' @export
geneRanges<-function(chrom,tss,strand,proxUp,proxDown,distal,
                     n=length(chrom)){
    # extendsMax is called once the basal domains are determined.
    # The bounds of each gene are set such that the extended
    # regions do not overlap with the basal regions of other genes
    # and such that the maximum distance is no greater than the
    # distal value from the TSS (transcription start site)
    extendsMax<-function(n,chrom,basalDomains,tss,distal){
        # limits comparision to each chromosome
        withinChrom<-(chrom[n]==chrom)
        # Normalize all basal domains to the TSS of the gene
        # of interest
        # case a sets the apporopriate lower maximum (lm)
        a<-basalDomains[withinChrom,1]-tss[n]
        # Check if there are genes between the basal bounds
        # and the max distance
        ap<-a[a> basalDomains[n,2]-tss[n] & a<distal]
        # Check if there are other overlaping bounds
        ab<-a[a>0 & a< basalDomains[n,2]-tss[n]]        
        if(length(ab)){
            # if there are overlapping bounds set the
            # lower maximum to the lower basal value
            lm<-basalDomains[n,2]
        }
        else if(length(ap)){
            # If there are no overlaping  domains but there
            # are other gene bounds found before the distal
            # limit set the lower maximum to the upper max
            # of the closest gene
            lm<-tss[n]+min(ap)
        }
        else{
            # If no geness are forund between the tss and
            # distal max set the value to the tss + max
            lm<-tss[n]+distal
        }
        # case b sets the apporopriate upper maximum (um)
        # The process is the same for a save for the compaision
        # of the upper max (um) rather than lower.
        b<-tss[n]-basalDomains[withinChrom,2]
        bp<-b[ b> tss[n]-basalDomains[n,1] & b<distal]
        bb<-b[b>0 & b< tss[n]-basalDomains[n,1]]
        if(length(bb))
            um<-basalDomains[n,1]
        else if(length(bp>0))
            um<-tss[n]-min(bp)
        else
            um<-tss[n]-distal
        # cons and return the upper and lower maximums
        geneDomain<-c(um,lm)        
        geneDomain
    }
    # Ensures that the bianary vector always goes from lowest
    # to greatest
    swapif<-function(x){
        if((x[1]<=x[2]))
            x
        else
            cbind(x[2],x[1])
    }
    if( ! all(levels(strand)==c("+","-"))){
        warning("Ensure that : levels(strand) == c(\"+\",\"-\")")
        return(NULL)
    }
    # Make sure that "+" strand -> -1 and "-" strand -> 1
    # Required to acuratly define basal domains
    levels(strand)=c(-1,+1)
    strand<-as.numeric(strand)
    # Extend each gene TSS by the upper proximal distance and 
    # lower proximal distance. (the direction is dependent on the 
    # strandedness of the gene) and then unsure that the upper and
    # lower values are in order    
    basalDomains<-t(apply(
        cbind(tss-strand*proxUp,tss+strand*proxDown),1,swapif))
    # Extend each of the genes basal domains and combine the
    # resulting list into a consL
    genomeDomains<-do.call(
        rbind,
        lapply(seq(n),extendsMax,chrom,basalDomains,tss,distal))
    return(genomeDomains)
}

#' Gene Range to Region
#'
#' Break the genome into regions based on the presence of genes. The
#' function relies on the output of geneRanges. It returnes mutualy
#' exclusive regions of the genome and the genes which relate to each
#' region.
#' @param ranges The genomic ranges for each gene
#' @param chrom A list of chromosomes, sequence should be shared with ranges
#' @param chroms a unique list of all chromomes in the crom list
#' @return A list of conses. <chrom><start><end><gene 1>...<gene N>
#' @export
geneRangeToRegion<-function(ranges,chrom,
                            chroms=unique(chrom)){
    # the output list chromosome Begining End Gene -> beg
    beg<-list()
    n=0
    # Flatens the two column ranges variable and relates each limit
    # to the initial order and wether the region begins or ends
    b<-cbind(c(as.numeric(ranges[,1]),as.numeric(ranges[,2])),
             seq(length(ranges[,1])),
             unlist(lapply(c(0,1),function(x)rep(x,length(ranges[,1])))))
    # each chromosome is assesed individually and then combined in the beg
    for(cro in chroms){    
        region=(chrom==cro)
        # seperate flattened list by chromosome
        k<-b[region,]
        # sort the flatend ranges by genomic coordinate
        # c has three values for each row, the first is the genomic
        # co-ordinate, the second is the gene number and the third is
        # the upper (1) or lowwer (0) bound info
        c<-k[order(k[,1]),]
        # the first region is assosited with the first gene ID
        buffer<-c(c[1,2])
        # reg is the temporary return gene list. It is recylced for each
        # chromosome
        reg<-list()
        # Init the geneomic co-ordinate
        start<-c[1,1]
        # The defined genome regions (defined by the limits of
        # the gene ranges) each have the genes which are pressent
        # assosiated with it. 
        for(i in seq(2,length(which(region))*2)){            
            reg<-append(reg,list(c(cro,start,c[i,1],buffer)))
            # if the region boundary represents an lower bound the gene is
            # added to the buffer list
            if(c[i,3]==0){                
                buffer<-c(buffer,c[i,2])
            }
            # if the boundary is an upper bound the gene is removed
            # from the buffer list
            else{
                buffer<-buffer[c[i,2]!=buffer]}
            # set the genomic co-ordinate
            start<-c[i,1]
        } 
        n<-n+1
        # filter the out the empty regions and add the chromosome regions to
        # the final return list
        beg<-append(beg,reg[unlist(lapply(reg,function(x) length(x)>3))])    
    }
    return(beg)
}

#' Genomic Regions
#' 
#' From a list of genes chromosomes tss and strandedness split the genome
#' into mutualy exclusive regions defined by the gene regions of influence
#' These regions can be used to relate ChIP binding data to genes. This
#' function combines the geneRagnges and geneRangeToRegion function into
#' a single form. 
#' @seealso geneRanges
#' @param chrom A list of chromosomes, one for each gene
#' @param tss The list of gene TSS
#' @param strand The stream directions of each of the genes
#' @param proxUp Upstream basal extension (independent of other genes)
#' @param proxDown Downstream basal extension (independent of other genes)
#' @param distal The limit of gene extension (in either direction)
#' @export
#' @examples
#' heightFile<-"peaks/jan/combined_heights.bed"
#' geneFile<-"~/hg19.RefSeqGenes.csv"
#' geneList<-read.delim(geneFile)
#' chrom<-as.character(geneList$chrom)
#' tss<-as.numeric(geneList$txStart)
#' strand<-geneList$strand
#' # double check to make sure that levels(strand) > c("-","+")
#' levels(strand)<-c(-1,1)
#' strand<-as.numeric(strand)
#' statsAndData<-readSplitTable(heightFile)
#' regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)
genomicRegions<-function(chrom,tss,strand,proxUp,proxDown,distal,
                         n=length(chrom),
                         chroms=unique(chrom)){
    # Find the gene ranges
    a<-geneRanges(chrom,tss,strand,proxUp,proxDown,distal,n)
    # Split the genome into regions based on gene ranges
    genes<-Filter(function(x){(as.numeric(x[2])>=0) || (as.numeric(x[3])>=0)} ,geneRangeToRegion(a,chrom,chroms))
    #set initial -ve values to 0
    genes<-lapply(genes,function(x) {if(as.numeric(x[2])<0){append(c(x[1],0), x[3:length(x)])}else{x}})
    genes
}

#' Great Gene Assoc
#'
#' Compare peak locations to predetermined genomic  regions and return
#' a list of genes that corospond to the overlaped gene regions. This
#' function acts as a filter for a gene list whose indicies are shared
#' with the gene regions. Only the genes within regions that have peaks
#' are returned.
#' @param bedData The binding locations of interest
#' @param genes The genome split into gene regions
#' @param geneList A list with a 1-1 relation to the gene regions
#' @return A filtered geneList
#' @examples
#' bedFile<-"peaks/jan/combined.bed"
#' bed<-read.delim(bedFile)
#' geneFile<-"~/hg19.RefSeqGenes.csv"
#' geneList<-read.delim(geneFile)
#' chrom<-as.character(geneList$chrom)
#' tss<-as.numeric(geneList$txStart)
#' strand<-geneList$strand
#' # double check to make sure that levels(strand) > c("-","+")
#' levels(strand)<-c(-1,1)
#' strand<-as.numeric(strand)
#' statsAndData<-readSplitTable(heightFile)
#' regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)
#' genes<-greatGeneAssoc(bed,regions,geneList)
#' @export
greatGeneAssoc <-function(bedData,genes,geneList){
    # seperate the info contained within the first three values
    # of each list value from the genes (gene regions)
    aChr<-as.character(lapply(genes,"[[",1))
    aStart<-as.numeric(lapply(genes,"[[",2))
    aEnd<-as.numeric(lapply(genes,"[[",3))
    # The summits are simplified to the center of the peaks
    peaks<-apply(bedData[,2:3],1,mean)
    # set the intial pc value, this value will be updated as the
    # bedData is iterated through, so dont count on it remaining
    # constant
    pc<-"chr0"
    # Find the subset of genomic regions which share their chromosome 
    # value with pc
    chroms<-(aChr==pc)
    # Isolate the subset of the starting and ending location for each
    # region on the genome
    cStart<-aStart[chroms]
    cEnd<-aEnd[chroms]
    locs<-c()
    # iterate through the peak data updating the subset of regions of
    # interest every time a new chromosome is encounterd. While the
    # bed data does not have to be sorted it does speed up computation
    # if it is.
    for (i in seq(length(bedData[,1]))){
        # if there is a chromosome change update the regions of interest
        if(pc!=as.character(bedData[i,1])){
            pc<-as.character(bedData[i,1])
            chroms<-which((aChr==pc))
            cStart<-aStart[chroms]
            cEnd<-aEnd[chroms]
        }
        # find the genes which have the peaks fall within them
        gene<-chroms[which(peaks[i]>cStart &peaks[i]<cEnd)]
        locs<-c(locs,gene)
    }
    # select the subset of the geneList which has had peaks bound
    geneList[as.numeric(unique(unlist(
        lapply(genes[unique(locs)],function(x) {x[4:length(x)]})))),]
}
