# function: NORMALIZE
#
row.col.normalize <- function( x, n=1, row.norm=n>1, col.norm=n>1, robust=FALSE, na.rm=TRUE )
{
    ## carrying out row and/or column normalization as many times as
    ## required

    ## handling of ExpressionSet object
    ##
    if ( class(x)=="ExpressionSet" ) {
        exprs(x)  <- row.col.normalize(exprs(x),n=n,row.norm=row.norm,col.norm=col.norm,robust=robust,na.rm=na.rm)
        return(x)
    }
    ## handling of matrix object
    ##
    if (  !row.norm && !col.norm )
        stop( "row.norm and col.norm can't be both false" )

    if ( is.vector(x) ) {
        x <- matrix( x, 1, length(x) )
    }
    if ( !row.norm || !col.norm ) n <- 1

    ctr <- if (robust) median else mean
    loc <- if (robust) mad else sd

    for ( i in (1:n) )
    {
        if ( row.norm )
        {
            ## row normalize
            ##
            x.m <- apply( x, 1, ctr, na.rm=na.rm )
            x.s <- apply( x, 1, loc, na.rm=na.rm )
            x <- ( x - x.m ) / x.s
        }
        if ( col.norm )
        {
            ## col normalize
            ##
            x.m <- apply( x, 2, ctr, na.rm=na.rm )
            x.s <- apply( x, 2, loc, na.rm=na.rm )
            x <- t( (t(x) - x.m)/x.s )
        }
    }
    return( x )
}
