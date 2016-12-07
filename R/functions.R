#' Function to generate example DEXSeqDataSet
#'
#' Function to generate a DEXSeqDataSet object based on data from
#' the pasilla package
#' 
#' @importClassesFrom DEXSeq DEXSeqDataSet
#' @importMethodsFrom DEXSeq estimateSizeFactors
#' @import DEXSeq
#' @import pasilla
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom statmod glmnb.fit
#' @importFrom stats deviance dnbinom fitted.values model.matrix optimize p.adjust pchisq
#' @importFrom utils data
#' @export

generateExampleDxd <- function()
{
    path <- system.file( package="pasilla", mustWork=TRUE )
    countFiles <- list.files( file.path( path, "extdata" ),
               pattern="fb.txt$",
               full.names=TRUE )
    flattenedfile = file.path( path, "extdata", "Dmel.BDGP5.25.62.DEXSeq.chr.gff")
    data("pasillaDEXSeqDataSet")
    sampleData <- as.data.frame( sampleAnnotation( dxd ) )
    names(countFiles) <- gsub(".txt", "", basename(countFiles))
    DEXSeqDataSetFromHTSeq(
         countfiles=countFiles[as.character(sampleData$sample)],
         sampleData = sampleData,
         design = ~ sample + exon + type:exon + condition:exon,
         flattenedfile=flattenedfile )
}


profileLogLikelihood <- function( disp, mm, y, muhat )
{
   # calculate the log likelihood:
   if(length(disp) != length(y)){
      disp <- rep(disp, length(y))
   }
   ll <- sum( sapply( seq(along=y), function(i)
      dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )
   # transform the residuals, i.e., y - muhat, to the linear
   # predictor scale by multiplying them with the derivative
   # of the link function, i.e., by 1/muhat, and add this to the
   # linear predictors, log(muhat), to get the predictors that
   # are used in the IWLS regression
   z <- log(muhat) + ( y - muhat ) / muhat
   # the variance function of the NB is as follows
   v0 <- muhat + disp * muhat^2
   # transform the variance vector to linear predictor scale by
   # multiplying with the squared derivative of the link to
   # get the (reciprocal) weights for the IWLS
   w <- 1 / ( ( 1 / muhat )^2 * v0 )
   # All we need from the IRLS run is the QR decomposition of
   # its matrix
   qrres <- qr( mm*sqrt(w) )
   # from it, we extract we leverages and calculate the Cox-Reid
   # term:
   cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )
   # return the profile log likelihood:
   ll - cr
}


estimateDispersionForExon <- function( modelMatrix, countVector,
                                      sizeFactors, dispInitialGuess = .5 )
{
   fit1 <- glmnb.fit( modelMatrix, countVector,
                     dispersion=dispInitialGuess,
                     offset=log(sizeFactors) )
   exp( optimize(
      function( logalpha )
         -profileLogLikelihood(
            exp(logalpha),
            modelMatrix,
            countVector,
            fitted.values(fit1) ),
      log( c( 1e-5, 1e3 ) ) )$minimum )
}

#' Function to estimate dispersions for each exon
#'
#' Given a DEXSeqDataSet object, it estimates dispersions
#' for each exon. The functions for shrinkage of dispersions
#' from DESeq2 are not implemented. Thus, it is recomended only
#' for comparisons with large numbers of biological replicates.
#' 
#' @param dxd A \code{DEXSeqDataSet} object.
#' @param BPPARAM A \code{BiocParallel} instance.
#' @param verbose If the function should be verbose or not.
#' @return A \code{DEXSeqDataSet} with dispersion estimations
#' included. Note that no sharing across exons is included, thus it
#' is recommended only for large datasets.
#'
#' @examples
#' data( "pasillaDEXSeqDataSet", package="pasilla" )
#' dxd <- estimateSizeFactors( dxd )
#' dxd <- estimateDispersionsAlt( dxd )
#'
#' @importFrom methods is
#' @export
#'

estimateDispersionsAlt <- function( dxd, BPPARAM = SerialParam(), verbose=TRUE )
{
    stopifnot( is( dxd, "DEXSeqDataSet") )
    if( any( is.na( colData(dxd)$sizeFactor ) ) ){
        stop("No size factors were found. Please run the function estimateSizeFactors first.")
    }
    modelMatrix <- DEXSeq:::rmDepCols(
        model.matrix( design(dxd), colData(dxd) ) )
    thisExons <- colData(dxd)$exon == "this"
    countMatrix <- counts(dxd)
    enoughCounts <- rowSums( countMatrix[,thisExons] ) > 0 & rowSums( countMatrix[,!thisExons] ) > 0 
    multiExonicGenes <- names( which( table(mcols(dxd)$groupID) > 1 ) )
    testableVector <- enoughCounts & mcols(dxd)$groupID %in% multiExonicGenes
    mcols( dxd )$allZero <- !testableVector
    countMatrix <- countMatrix[testableVector,]
    toSplit <- sort( rep( 1:bpnworkers(BPPARAM), length.out=nrow(countMatrix) ) )
    spMat <- split( as.data.frame( countMatrix ), toSplit )
    names( spMat ) <- NULL
    sizeFactors <- colData(dxd)$sizeFactor
    dispsAll <- bplapply( spMat, function(x, testableVectorZ, mmZ, sizeFactorsZ, verboseZ ){
        disps <- sapply( rownames(x), function(ex){
            if( verboseZ ){
                cat(sprintf("fitting %s\n", ex))
            }
            counts <- as.numeric(x[ex,])
            dp <- try( estimateDispersionForExon( mmZ, counts, sizeFactorsZ) )
            if( !inherits( dp, "try-error" ) ){
               dp
            }else{
               NA
            }
        } )
        names( disps ) <- rownames( x )
        disps
    },
    mmZ=modelMatrix, sizeFactorsZ=sizeFactors, verboseZ=verbose, BPPARAM=BPPARAM)
    dispsAll <- unlist(dispsAll)
    mcols( dxd )$dispersion <- NA
    mcols(dxd)$dispersion[match( names( dispsAll), rownames( dxd ) )] <- dispsAll
    dxd
}

#' Function to test for differences in exon usage
#'
#' This function implements testing for differences in exon
#' usage across biological conditions. It is an alternative
#' implementation of the original statistical method implemented
#' in the DEXSeq package. The only difference is that it uses
#' the stamod glmnb.fit instead of the DESeq2 glm fitter.
#'
#' @param dxd A \code{DEXSeqDataSet} object.
#' @param BPPARAM A \code{BiocParallel} instance.
#' @param verbose If the function should be verbose or not.
#' @param reducedModel Formula with the reduced model
#' @param fullModel Formula with the full model
#'
#' @return A \code{DEXSeqDataSet} with p-values and q-values.
#'
#' @examples
#' data( "pasillaDEXSeqDataSet", package="pasilla" )
#' dxd <- estimateSizeFactors( dxd )
#' dxd <- estimateDispersionsAlt( dxd )
#' 
#' @export
#'
#' 

testForDEUAlt <- function( dxd, reducedModel, fullModel, BPPARAM=SerialParam(), verbose=TRUE ){
    modelFrame <- as.data.frame( colData( dxd ) )
    countMatrix <- counts( dxd )
    mmNull <- DEXSeq:::rmDepCols( model.matrix( reducedModel, modelFrame ) )
    mmFull <- DEXSeq:::rmDepCols( model.matrix( fullModel, modelFrame ) )
    countMatrix <- countMatrix[ !is.na( mcols(dxd)$dispersion ),]
    toSplit <- sort( rep( seq_len( bpnworkers( BPPARAM ) ),
                         length.out=nrow( countMatrix ) ) )
    splitMat <- split( as.data.frame( countMatrix ), toSplit )
    names( splitMat ) <- NULL
    rawDisps <- mcols(dxd)$dispersion
    names(rawDisps) <- rownames(dxd)
    cat(sprintf("testing using %s cores\n", bpnworkers(BPPARAM) ) )
    pvalAll <- bplapply( splitMat, function( mat, modelFrameZ, rawDispsZ, mmNullZ, mmFullZ, verboseZ ){
        pvals <- sapply( rownames( mat ), function( ex ){
            if( verboseZ ){
                cat( sprintf( "testing %s\n", ex ) )
            }
            x <- which( rownames( mat ) %in% ex )
            xi <- which( names( rawDispsZ ) %in% ex )
            counts <- as.numeric( mat[x,] )
            disps <- rep( rawDispsZ[xi], length(counts) ) 
            sf <- modelFrameZ$sizeFactor
            fitNull <- try( glmnb.fit( X=mmNullZ, y=counts, dispersion=disps, offset=log( sf ) ), silent=TRUE )
            fitFull <- try( glmnb.fit( X=mmFullZ, y=counts, dispersion=disps, offset=log( sf ) ), silent=TRUE )
            if( inherits( fitNull, "try-error") | inherits( fitFull, "try-error") ){
                return(NA)
            }
            pval <- 1 - pchisq( deviance( fitNull ) - deviance( fitFull ),
                           df=ncol( mmFullZ ) - ncol( mmNullZ ) )
            return( pval )
        } )
        names( pvals ) <- rownames( mat )
        pvals
    }, modelFrameZ = modelFrame, rawDispsZ=rawDisps, verboseZ=verbose,
    mmNullZ = mmNull,
    mmFullZ = mmFull, BPPARAM = BPPARAM )
    pvalAll <- unlist( pvalAll )
    mcols(dxd)$LRTPvalue <- NA
    mcols(dxd)$LRTPvalue[match( names( pvalAll ), rownames( dxd ) )] <- pvalAll
    mcols(dxd)$qvalue <- p.adjust( mcols(dxd)$LRTPvalue, method="BH" )
    mcols(mcols(dxd))[colnames(mcols(dxd)) %in% c("LRTPvalue", "qvalue"), "type"] <- "results"
    attr(dxd, "test") <- "LRT"
    dxd
}
