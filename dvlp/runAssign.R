run_assign <- function(ES,   # ExpressionSet object 
                       eSig, # Gene signature used by ASSIGN (vector of valid gene symbols)
                       oDir, # Output directory where ASSIGN output will be stored
                       iter=3000 # Number of iterations to perform for the MCMC chain run
                       )
{
  if (class(ES)!="ExpressionSet") stop( "ExpressionSet expected: ", class(ES) )
  if (class(eSig)!="list") stop( "eSig must be a list: ", class(eSig) )
  
  cat("Running ASSIGN on Expression Set: ",substitute(ES),"\n")
  cat("Running with signature: ",substitute(eSig),"\n")
  cat("Number of genes in signature: ",length(eSig),"\n")
  cat("Saving output to directory: ",eval(oDir),"\n")
  
  ## Create output directory
  dir.create(oDir,showWarnings=FALSE)
  
  ## Start assign
  processed.data <- assign.preprocess(trainingData=NULL,
                                      testData=exprs(ES),
                                      trainingLabel=NULL,
                                      geneList=eSig,
                                      n_sigGene=NA,
                                      theta0=0.05,
                                      theta1=0.9)
  
  ## Run mcmc
  mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, 
                            Bg = processed.data$B_vector, 
                            X=processed.data$S_matrix, 
                            Delta_prior_p = processed.data$Pi_matrix, 
                            iter = iter, 
                            adaptive_B=TRUE, 
                            adaptive_S=TRUE, 
                            mixture_beta=TRUE)
  
  mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=1000, 
                                  iter=iter, adaptive_B=TRUE, 
                                  adaptive_S=TRUE,mixture_beta=TRUE)
  
  assign.output(processed.data=processed.data, 
                mcmc.pos.mean.testData=mcmc.pos.mean, 
                trainingData=NULL, 
                testData=ES, 
                trainingLabel=NULL, 
                testLabel=NULL, 
                geneList=eSig, 
                adaptive_B=TRUE, 
                adaptive_S=TRUE, 
                mixture_beta=TRUE, 
                outputDir=oDir)
  ## with output function
  output.data <- list(processed.data = processed.data, 
                      mcmc.pos.mean.testData = mcmc.pos.mean)
  save(output.data, file = paste(oDir,"output.rda",sep='/'))
  
  cat("Done!\n")  
}
