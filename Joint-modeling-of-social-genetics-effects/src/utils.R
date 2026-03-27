
#####                  FUNCTIONS TO SIMULATE GENOS DATA           #####################
##' Simulates genotypes with the sequential coalescent with recombination model
##' @return list
##' @author Timothee Flutre
##' @examples
##' set.seed(12345)
##' level_mig_rate <- "med"
##' genomes <- simulGenos(level_mig_rate=level_mig_rate)
##' M <- genomes$genos
##' out_pca <- pca(M)
##' plotPca(rotation=out_pca$rot.dat,
##'         prop.vars=out_pca$prop.vars,
##'         main=paste0("PCA (migration= ", level_mig_rate, ")"))
##' afs <- estimSnpAf(M)
##' plotHistAllelFreq(afs=afs,
##'                   main=paste0("Allele frequencies (migration=",
##'                   level_mig_rate, ")"))
##' mafs <- estimSnpMaf(M)
##' plotHistMinAllelFreq(mafs=mafs,
##'                      main=paste0("Minor allele frequencies (migration=",
##'                      level_mig_rate, ")"))
##'
##' ## estimate additive genetic relationships
##' A_vanraden <- estimGenRel(M, afs=afs, relationships="additive",
##'                           method="vanraden1", verbose=0)
##' hist(diag(A_vanraden), # 1 expected under HWE
##'      breaks="FD",
##'      main=paste0("Inbreeding coefficients (migration=",
##'      level_mig_rate, ")"))
##' hist(A_vanraden[upper.tri(A_vanraden)], # 0 expected under HWE
##'      main=paste0("Additive genetic relationships (migration=",
##'      level_mig_rate, ")"))
##' A_noia <- estimGenRel(M, afs=afs, relationships="additive",
##'                       method="noia", verbose=0)
##'
##' ## estimate LD between SNPs with a MAF high enough:
##' min_maf <- 0.15
##' length(snps_tokeep <- colnames(M[, mafs >= min_maf]))
##' out_LD <- estimLd(X=M[, snps_tokeep], snp.coords=genomes$snp.coords)
##' out_LD$dist <- distSnpPairs(data.frame(loc1=out_LD$loc1, loc2=out_LD$loc2),
##'                             genomes$snp.coords,
##'                             nb.cores=1, verbose=1)
simulGenos <- function(nb_genos=300, nb_chroms=10, chrom_len=10^5,
                       Ne=10^4, mut=10^(-8), rec=10^(-8),
                       nb_pops=10, level_mig_rate="med", ind.ids = NULL,
                       verbose=1){
  stopifnot(requireNamespace("rutilstimflutre"),
            requireNamespace("scrm"))
  mig_rates <- c("high"=10^2, "med"=10, "low"=0.5)
  genomes <- rutilstimflutre::simulCoalescent(
    nb.inds=nb_genos,
    ind.ids = ind.ids,
    nb.reps=nb_chroms,
    pop.mut.rate=4 * Ne * mut * chrom_len,
    pop.recomb.rate=4 * Ne * rec * chrom_len,
    chrom.len=chrom_len,
    nb.pops=nb_pops,
    mig.rate=mig_rates[level_mig_rate],
    verbose=verbose)
  return(genomes)
}


##' #check package function
##'
packageCheck <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required"))
  }
}


#################     FUNCTION TO SIMULATE YIELD      ##########################
##' Make the design matrix of the DBVs
##' @param df is dataframe
##' @param species the column with the species name
##' @param focal is col with the focal genotypes name
##' @param Spec column specified the species (Wheat or Pea) to compute DBVs matrix for
##' @skipUnusedCols logical value (default is TRUE); when true, column 0 value is remove in the output
##' @return matrix
##' @author Jemay Salomon
##' @examples
##' ## Create a fake data frame
##' dat <- data.frame(
##' focal = c("wheat01", "pea02", "wheat01"),
##'  neighbors = c("pea02", "wheat01", "pea01"),
##' name = c("wheat01-pea02", "wheat01-pea02", "wheat01-pea01"),
##'  stand = c("mixed", "mixed", "mixed"),
##'   species = c("wheat", "pea", "wheat")
##' )
##' ## Print dataframe
##' print(dat)
##' (Z_DBV=mkZDBV(df=dat, colD="focal", colS="species", species="wheat",
##'                    skipUnusedCols=TRUE))

mkZDBV <- function(df, colD, colS, species, skipUnusedCols = TRUE) {
  require(gtools)
  stopifnot(is.data.frame(df),
            colD %in% colnames(df),
            colS%in%colnames(df))

  idx <- which(df[[colS]] == species)
  genos <- sort(unique(as.character(df[[colD]][idx])))
  nbGenos <- length(genos)
  Z_DBV <- matrix(0, nrow=nrow(df), ncol=nbGenos)
  colnames(Z_DBV) <- genos

  for(geno in genos){
    idxGeno <- which(df[[colS]] == species & df[[colD]] == geno)
    Z_DBV[idxGeno, geno] <- 1
  }

  if (skipUnusedCols) {
    Z_DBV<- Z_DBV[, colSums(Z_DBV) != 0, drop = FALSE]
  }
  Z_DBV<- Z_DBV[, mixedsort(colnames(Z_DBV))]
  stopifnot(is.matrix(Z_DBV))
  return(Z_DBV)
}


##' Make the design matrix of the SIGV effects
##' @param df is dataframe
##' @param focal is col with the the focal genotypes name
##' @param stand name of the column indicating the plots status (mixed or pure)
##' @param species the column with the species name
##' @param Spec column specified the species (Wheat or Pea) to compute SBVs matrix for
##' @param model3Emma  Logical value (default is true); when true, the intra specific social effects (SIS) is added
##' @param weightW the intensity of sis  for wheat when model3Emma  is true; default is 0.25
##' @param weightP the intensity of sis  for pea when model3Emma  is true; default is 0.25
##' @skipUnusedCols: Logical value (default is TRUE); when true, column 0 value is remove in the output
##' For more details, see \href{https://doi.org/10.1016/j.fcr.2019.107571}{Forst et al (2019)}
##' @return a matrix
##' @author Jemay SALOMON
##' @examples
##' ## Create a dataframe
##'  dat <- data.frame(
##'    focal = c("wheat01", "pea02", "wheat01", "wheat01"),
##'    neighbors = c("pea02", "wheat01", "pea01", "wheat01"),
##'    name = c("wheat01-pea02", "wheat01-pea02", "wheat01-pea01", "wheat01"),
##'    stand = c("mixed", "mixed", "mixed", "pure"),
##'    species = c("wheat", "pea", "wheat", "wheat")
##'  )
##'  print(dat)
##'  (Z_SIGV=mkZSIGV(df=dat, focal="focal", Spec="pea", species="species", stand="stand",skipUnusedCols=TRUE))
# mkZSIGV <- function(df, focal, species, stand, Spec, weightP=0.25,
#                    weightW=0.25, mixed = FALSE, skipUnusedCols=TRUE){
#
#   stopifnot(is.data.frame(df),
#             focal %in% colnames(df),
#             is.logical(skipUnusedCols),
#             is.logical(mixed))
#
#   genos <- df[[focal]]
#   unique_genos <- unique(genos)
#   Z_SIGV <- matrix(0, nrow = nrow(df), ncol = length(unique_genos))
#   colnames(Z_SIGV) <- unique_genos
#   current_species=df[[species]]
#   current_stand=df[[stand]]
#
#   for(i in 1:nrow(Z_SIGV)){
#     if (Spec=="pea"){
#       if (current_species[i]=="pea"){
#         if(current_stand[i]=="pure"){
#           Z_SIGV[i, which(unique_genos == genos[i])] <- 1
#         } else if (current_stand[i]=="mixed"){
#           if(mixed){
#             Z_SIGV[i, which(unique_genos == genos[i])] <- 1*weightP
#           }
#         }
#       }
#     } else if (Spec=="wheat"){
#       if (current_species[i]=="wheat"){
#         if(current_stand[i]=="pure"){
#           Z_SIGV[i, which(unique_genos == genos[i])] <- 1
#         } else if (current_stand[i]=="mixed"){
#           if(mixed){
#             Z_SIGV[i, which(unique_genos == genos[i])] <- 1*weightW
#           }
#         }
#       }
#     }
#
#   }
#   if (skipUnusedCols) {
#     Z_SIGV<- Z_SIGV[, colSums(Z_SIGV) != 0, drop = FALSE]
#   }
#   Z_SIGV<- Z_SIGV[, mixedsort(colnames(Z_SIGV))]
#   stopifnot(is.matrix(Z_SIGV))
#   return(Z_SIGV)
# }
mkZSIGV <- function(df, focal, species,
                    stand, mkZForSpec,
                    weight=0.25,
                    mixed = FALSE, skipUnusedCols=TRUE){

  stopifnot(is.data.frame(df),
            focal %in% colnames(df),
            is.logical(skipUnusedCols),
            is.logical(mixed))

  genos <- df[[focal]]
  unique_genos <- unique(genos)
  Z_SIGV <- matrix(0, nrow = nrow(df), ncol = length(unique_genos))
  colnames(Z_SIGV) <- unique_genos
  current_species=df[[species]]
  current_stand=df[[stand]]

  for(i in 1:nrow(Z_SIGV)){
    if (current_species[i]==mkZForSpec){
      if(current_stand[i]=="pure"){
        Z_SIGV[i, which(unique_genos == genos[i])] <- 1
      } else if (current_stand[i]=="mixed"){
        if(mixed){
          Z_SIGV[i, which(unique_genos == genos[i])] <- 1*weight
        }
      }
    }
  }
  if (skipUnusedCols) {
    Z_SIGV<- Z_SIGV[, colSums(Z_SIGV) != 0, drop = FALSE]
  }
  Z_SIGV<- Z_SIGV[, mixedsort(colnames(Z_SIGV))]
  stopifnot(is.matrix(Z_SIGV))
  return(Z_SIGV)
}


getNeighbors <- function(df, name, sep){
  neighbors <- strsplit(as.character(df[[name]]), sep)
  names(neighbors) <- as.character(df[[name]])
  neighbors <- lapply(neighbors, unique)
  return(neighbors)
}

getGenos <- function(neighbors){
  return(mixedsort(unique(do.call(c, neighbors))))
}

##' Make the design matrix of the DBV-SBV interactions
##' @param df is dataframe
##' @param name name of the column containing the genotypes per pure stand or mixture
##' @param sep separator
##' @skipUnusedCols: Logical value (default is TRUE); when true, column 0 value is remove in the output
##' @return a  matrix
##' @author Jemay SALOMON
##' @examples
##' \dontrun{
##' ## Create a dataset
##' dat <- data.frame(
##'    focal = c("wheat01", "pea02", "wheat01"),
##'    neighbors = c("pea02", "wheat01", "pea01"),
##'    name = c("wheat01-pea02", "wheat01-pea02", "wheat01-pea01"),
##'    stand = c("mixed", "mixed", "mixed"),
##'    species = c("wheat", "pea", "wheat")
##'  )
##'
##' ## Print dataframe
##' print(dat)
##'
##'  (ZDBV.SBV=mkZDBV.SBV(df=dat, col="name", sep="-"))
###' }
mkZDBV.SBV <- function(df, sep, name, species,Spec="both",skipUnusedCols = TRUE) {
  require(gtools)
  stopifnot(is.data.frame(df),
            name %in% colnames(df),
            is.logical(skipUnusedCols))


  species <- df[[species]]

  neighbors <- getNeighbors(df, name, sep)
  genos <-getGenos(neighbors)
  colNames <- expand.grid(genos,genos)
  colNames <- apply(colNames, 1, function(x) paste((x), collapse = sep))
  colNames <- colNames[!duplicated(colNames)]

  ZDBV.SBV <- matrix(0, nrow = nrow(df), ncol = length(colNames))
  colnames(ZDBV.SBV) <- colNames


  allPairs <- lapply(neighbors, function(x) {
    if (length(x) == 1) {
      paste(rep(x, 2), collapse = sep)
    } else {
      tmp <- t(utils::combn((x), 2))
      paste(mixedsort(tmp[, 1]), mixedsort(tmp[, 2]), sep = sep)
    }
  })


  for (i in 1:nrow(ZDBV.SBV)) {
    conditions <- switch(
      Spec,
      "both" = TRUE,
      "pea" =   (species[i] =="pea"),
      "wheat" = (species[i] =="wheat"),
      stop("Invalid value for Spec")
    )

    ZDBV.SBV[i, allPairs[[i]]] <- as.integer(conditions)
  }

  for (name in 1:ncol(ZDBV.SBV)) {
    column_name_components <- unlist(strsplit(colnames(ZDBV.SBV)[name], sep))
    if (column_name_components[1] == column_name_components[2]) {
      ZDBV.SBV[, name] <- 0
    }
  }

  if (skipUnusedCols) {
    isEmpty <- colSums(ZDBV.SBV) == 0
    isPure <- sapply(strsplit(colnames(ZDBV.SBV), sep), anyDuplicated) == 2
    colsToKeep <- !isPure & !isEmpty
    ZDBV.SBV <- ZDBV.SBV[, colsToKeep, drop = FALSE]
  }
  ZDBV.SBV <- ZDBV.SBV[, mixedsort(colnames(ZDBV.SBV))]

  stopifnot(is.matrix(ZDBV.SBV))
  return(ZDBV.SBV)
}


##' @param x error variance for wheat species
##' @return a logical
##' @author Jemay Salomon
##' @examples
##'
isDecimal <- function(x) {
  return(x != trunc(x))
}

## Function to simulate design
##' @author Jemay SALOMON
simulDesign <- function(strategies, panelGenos, testerGenos,
                        propsole_onlypur, propsole_inter_50pur, propinter_onlytester1,
                        propsole_inter_50tester1, block, seed = NULL)
{

  nbPanelGenos<- length(panelGenos)
  propDtester1 = ifelse(strategies=="inter_only", propinter_onlytester1, propsole_inter_50tester1)

  ## stopifnot
  stopifnot(!is.null(panelGenos))
  if(!(strategies=="sole_only")){
    stopifnot(!(is.null(testerGenos)),
              !(isDecimal(nbPanelGenos* propDtester1)),
              !(isDecimal(nbPanelGenos* (1 - propDtester1))))
  }

  if(! is.null(panelGenos)) stopifnot(is.character(panelGenos))
  if(! is.null(testerGenos)) stopifnot(is.character(testerGenos))
  if(!is.factor(panelGenos)) panelGenos <- factor(panelGenos)
  if(!is.factor(testerGenos)) testerGenos <- factor(testerGenos)

  #set.seed
  if (!is.null(seed)) set.seed(seed)

  #testers
  testers <- sample(testerGenos, length(unique(testerGenos)), replace = FALSE)


  #blocks
  blocks <- LETTERS[1:block]
  levblocks <- blocks

  #helper function to sample genotypes
  getMixedGenotypes <- function(prop, panelGenos) {
    nbMixStand <- nbPanelGenos* prop
    propGenoDmix <- mixedsort(sample(panelGenos, nbMixStand, replace = FALSE))
    return(propGenoDmix)
  }

  #define global-local variables
  propDpur = ifelse(strategies=="sole_only", propsole_onlypur, propsole_inter_50pur)
  focal <- mixedsort(as.character(sample(panelGenos, size = nbPanelGenos*propDpur, replace = FALSE)))
  nbMixStand <- (nbPanelGenos*propDtester1) + (nbPanelGenos*(1-propDtester1))
  propGenoDmix1 <- mixedsort(as.character(getMixedGenotypes(propDtester1, panelGenos)))
  r_panelGenos <- panelGenos[!(panelGenos %in% propGenoDmix1)]
  propGenoDmix2 <- mixedsort(as.character(getMixedGenotypes(1-propDtester1, r_panelGenos)))

  #Loops to complete dataframe
  if(strategies=="sole_only"){ # strategies 1
    dat.block1 <- data.frame(stand=c(rep("pure",nbPanelGenos*propsole_onlypur)),
                             focal=c(mixedsort(focal)),
                             block=NA)
    dat <- dat.block1
    for(k in 2:block){
      dat <- rbind(dat, dat.block1)
    }

    dat$block <- as.factor(rep(blocks, each=nrow(dat.block1)))
  }else if(strategies== "inter_only"){

    #strategies 2
    dat.block1 <- data.frame(stand=c(rep("mixed",nbMixStand)),
                             focal=c(rep(mixedsort(propGenoDmix1), length(testers[1])),
                                     rep(mixedsort(propGenoDmix2), length(testers[2]))),
                             neighbors=c(rep(testers[1], length(propGenoDmix1)),
                                         rep(testers[2],length(propGenoDmix2))),
                             block= NA)
    dat <- dat.block1

    for(k in 2:block){
      dat <- rbind(dat, dat.block1)
    }
    dat$block <- rep(blocks, each=nrow(dat.block1))

  }else if(strategies== "sole_inter_50"){# strategies 3

    stopifnot(!isDecimal(nbMixStand/4))

    dat.block1 <- data.frame(stand=c(rep("pure",length(focal)),
                                     rep("mixed", nbMixStand)),

                             focal=c(focal,
                                     rep(propGenoDmix1, length(testers[1])),
                                     rep(propGenoDmix2, length(testers[2]))),

                             neighbors=c(as.factor(rep(NA, length(focal))),
                                         rep(testers[1], length(propGenoDmix1)),
                                         rep(testers[2], length(propGenoDmix2))),

                             block=NA)
    dat <- dat.block1

    dat$block <- c(rep(blocks[1], length(focal)/block),
                   rep(blocks[2], length(focal)/block),
                   rep(blocks[1], nbMixStand/4),
                   rep(blocks[2], nbMixStand/4),
                   rep(blocks[1], nbMixStand/4),
                   rep(blocks[2], nbMixStand/4))



  } else if(strategies=="inter_only_ctrl"){

    #strategies 2
    dat.block1 <- data.frame(stand=c(rep("mixed",length(panelGenos)*length(testerGenos))),
                             focal=c(rep(panelGenos, each=length(testerGenos))),
                             neighbors=c(rep(testerGenos, length(panelGenos))),
                             block= NA)
    dat <- dat.block1

    for(k in 2:block){
      dat <- rbind(dat, dat.block1)
    }
    dat$block <- rep(blocks, each=nrow(dat.block1))

  } else if(strategies=="sole_inter_ctrl"){
    dat.block1 <- data.frame(stand=c(rep("pure",length(panelGenos)),
                                     rep("mixed",length(panelGenos)*length(testerGenos))),
                             focal=c(panelGenos, rep(panelGenos, each=length(testerGenos))),
                             neighbors=c(as.factor(rep(NA, length(panelGenos))),
                                         rep(testerGenos, length(panelGenos))),
                             block=NA)
    dat <- dat.block1

    for(k in 2:block){
      dat <- rbind(dat, dat.block1)
    }
    dat$block <- rep(blocks, each=nrow(dat.block1))

  } else{
    stop("strategies must be  sole_only, inter_only, inter_only_ctrl, sole_inter_50 & sole_inter_ctrl")
    break
  }

  if (!(strategies == "sole_only")) {
    dat$name <- paste0(dat$focal, "_", dat$neighbors)
    dat$name <- gsub("NA_|_NA", "", dat$name)
    dat$name <- factor(dat$name)
    dat$neighbors<- factor(dat$neighbors, levels = c("pea01", "pea02"))
  }
  dat$block <- factor(dat$block, levels = levblocks)
  dat$stand <- as.factor(dat$stand)
  dat$focal <- as.factor(dat$focal)
  return(dat)
  invisible(gc()) #free up memory
}

#get sigma2 err
getSigma2Err <- function(sigma2Geno, h2, nbBlocks){
  ((1 - h2) / h2) * (nbBlocks*sigma2Geno)
}


## function to simulate yield in pur to be used in sole_only and sole_inter_50
simPur <- function(datPur, listZs, listXs,
                   listVCov,
                   paramsScenario,
                   explFacts, strategy,
                   simulId){

  set.seed(paramsScenario$simulSeeds[simulId])

  bEffs <- c(paramsScenario$block1SpeciesEffs["mu","pure", "wheat"],
             paramsScenario$block1SpeciesEffs["block1","pure","wheat"])

  out <- list()

  rownames(datPur) <- NULL

  datPur$species <- ifelse(datPur$focal %in% paramsScenario$panelGenos, "wheat", "pea")
  datPur$species <- as.factor(datPur$species)

  listZs[["Z_DBV"]] = mkZDBV(df=datPur, colD="focal",
                             colS="species", species="wheat",
                             skipUnusedCols=TRUE)

  ##ZSIS
  listZs[["Z_SIGV"]] = mkZSIGV(df = datPur, focal = "focal",
                               species = "species",
                               stand="stand",
                               mkZForSpec = "wheat",
                               skipUnusedCols = TRUE, mixed = FALSE)

  ## design matrix for pur
  listXs[["X"]] <- model.matrix(~ 1 + block, datPur,
                                contrasts.arg=list("block"="contr.sum"))


  ## errors
  Idn_panelGenos_pur <- diag(nrow(datPur))
  R <-paramsScenario[["sigma2_err_w"]]*Idn_panelGenos_pur
  e <- mvrnorm(n=1, mu=rep(0,nrow(Idn_panelGenos_pur)), Sigma=R)


  #print dimensions
  ##########################################
  # cat("dimX \n")
  # print(dim(listXs[["X"]]))
  #
  # cat("dimZD \n")
  # print(dim(listZs[["ZD"]]))
  #
  # cat("dimZS \n")
  # print(dim(listZs[["ZS"]]))
  #
  # cat("dimXBVs \n")
  # print(dim(explFacts$U_DBV_SBV))
  #
  # cat("length bEFFs \n")
  # print(length(bEffs))
  #
  # cat("length SIS \n")
  # print(length(explFacts$SIGVpanelGenos))
  #############################################


  ##Check
  stopifnot(paramsScenario$panelGenos==colnames(listZs[["Z_DBV"]]))
  stopifnot(paramsScenario$panelGenos==colnames(listZs[["Z_SIGV"]]))
  stopifnot(ncol(listXs[["X"]])==length(bEffs))
  stopifnot(ncol(listZs[["Z_DBV"]])==nrow(explFacts$U_DBV_SBV[,1]))

  y_panelGenos_pur <- listXs[["X"]] %*% bEffs + listZs[["Z_DBV"]]%*%explFacts$U_DBV_SBV[,1] +
    listZs[["Z_SIGV"]] %*% explFacts$SIGVpanelGenos + e


  #append
  y_panelGenos_pur  <- y_panelGenos_pur [,1]
  datPur$yield <- y_panelGenos_pur
  rownames(datPur) <- NULL


  out <- list(data = datPur,
              y = y_panelGenos_pur,
              purZs = listZs,
              purXs = listXs)

  return(out)
}


##Draw MVN Errors
simE =  function(n_1, d, paramsScenario, simulId){

  set.seed(paramsScenario$simulSeeds[simulId])

  stopifnot(requireNamespace("MixMatrix"))

  M_E <- matrix(0, n_1, d,dimnames=list(NULL, paramsScenario$levSpecies))
  Id_n1 <- diag(n_1)
  Sigma_E <- matrix(c(paramsScenario[["sigma2_err_w_mix"]], NA, NA, paramsScenario[["sigma2_err_p_mix"]]), 2, 2)
  Sigma_E[1,2] <- Sigma_E[2,1] <- paramsScenario[["cor_E"]] * sqrt(Sigma_E[1,1] * Sigma_E[2,2])
  cov_Er <- paramsScenario[["cor_E"]] * sqrt(Sigma_E[1,1] * Sigma_E[2,2])
  dimnames(Sigma_E) <- list(paramsScenario$levSpecies, paramsScenario$levSpecies)
  paramsScenario[["Sigma_E"]] <- Sigma_E
  ## E <- rutilstimflutre::rmatnorm(n=1, M_E, U = Id_n1, V = Sigma_E)[,,1]
  E <- MixMatrix::rmatrixnorm(n=1, mean=M_E, U=Id_n1, V=Sigma_E)

  return(E)
}


simulExplFactors  <- function(paramsScenario, kinship, simulId){

  stopifnot(is.list(kinship),
            all(c("K_f","Kmix") %in% names(kinship)))

  set.seed(paramsScenario$simulSeeds[simulId]) #reproductibility

  stopifnot(requireNamespace("MASS"),
            requireNamespace("MixMatrix"))

  if(!is.null(kinship$K_f)){
    if(!is.matrix(kinship$K_f)&&!is.numeric(kinship$K_f)){
      stop("K must be a numeric matrix")
    }
  }

  #chek
  stopifnot(colnames(kinship$K_f)==paramsScenario$panelGenos)
  stopifnot(rownames(kinship$K_f)==paramsScenario$panelGenos)

  GE <- list()

  #U_DBV_SBV
  nbPanelGenos <- length(paramsScenario$panelGenos)

  M = matrix(0, ncol=2, nrow = nbPanelGenos)
  rownames(M)= paramsScenario$panelGenos
  colnames(M)=c("DBV", "SBV")
  cov.DBVw.SBVw <- paramsScenario[["cor_DBVw_SBVw"]] *
    sqrt(paramsScenario[["sigma2.DBVw"]] *
           paramsScenario[["sigma2.SBVw"]])

  G = matrix(c(paramsScenario[["sigma2.DBVw"]], cov.DBVw.SBVw,
               cov.DBVw.SBVw, paramsScenario[["sigma2.SBVw"]]),
             byrow=FALSE, nrow=2, ncol=2)


  colnames(G)=c("DBV", "SBV")
  rownames(G)=c("DBV", "SBV")

  # print(kinship[["K_f"]][1:5, 1:5])

  paramsScenario[["varcov_DBV_SBV_w"]] = G
  U_DBV_SBV  <- rutilstimflutre::rmatnorm(1, M, kinship$K_f , G)[,,1]
  # U_DBV_SBV  <- MixMatrix::rmatrixnorm(n=1, mean=M, U=kinship$K_t, V=G)
  U_DBV_SBV <-round(data.frame(U_DBV_SBV), 3)
  U_DBV_SBV = as.matrix(U_DBV_SBV)

  GE[["U_DBV_SBV"]] = U_DBV_SBV

  #SIS
  Q <- paramsScenario[["sigma2.SIGVw"]] * kinship$K_f
  SIGVpanelGenos =  MASS::mvrnorm(n=1, mu=rep(0,nbPanelGenos), Sigma=Q)
  names(SIGVpanelGenos) <- paramsScenario$panelGenos
  GE[["SIGVpanelGenos"]] = SIGVpanelGenos

  #U_DBV_x_SBV
  GE[["Kmix"]] = kinship[["Kmix"]]
  M_DBV_x_SBV = matrix(0, ncol=2, nrow=nrow(kinship[["Kmix"]]))
  colnames(M_DBV_x_SBV)=c("DBV.SBV.wp", "DBV.SBV.pw")
  rownames(M_DBV_x_SBV)=rownames(kinship[["Kmix"]])

  I = matrix(c(paramsScenario[["sigma2.DBV.SBV.wp"]],
               paramsScenario[["cov_DBV_x_SBV"]],
               paramsScenario[["cov_DBV_x_SBV"]],
               paramsScenario[["sigma2.DBV.SBV.pw"]]),
             byrow=FALSE, nrow=2, ncol=2)

  # cat("var-cov-G\n")
  # print(I)
  # print(kinship[["Kmix"]][1:5, 1:5])

  colnames(I)=c("DBV.SBV.wp", "DBV.SBV.pw")
  rownames(I)=c("DBV.SBV.wp", "DBV.SBV.pw")
  U_DBV_x_SBV  <- rutilstimflutre::rmatnorm(1, M_DBV_x_SBV, kinship[["Kmix"]] , I)[,,1]
  U_DBV_x_SBV <-round(data.frame(U_DBV_x_SBV), 3)
  U_DBV_x_SBV = as.matrix(U_DBV_x_SBV)
  GE[["U_DBV_x_SBV"]] = U_DBV_x_SBV

  # cat("var-cov-from data\n")
  # print(cov(U_DBV_x_SBV))
  # print(head(U_DBV_x_SBV))
  # space(2)

  return((GE))

}


###simulate data in mixed stand
simMix <- function(datmix, listZs, listXs, listVCov,
                   paramsScenario,
                   explFacts, strategy, simulId){

  #Fixed effects of block and mu in mix
  set.seed(paramsScenario$simulSeeds[simulId])
  B_mix <- paramsScenario$block1SpeciesEffs[,"mixed",]

  #DBV-SBV pea
  U_p <- matrix(c(paramsScenario$DBVSBVEffs[["neighbors1"]], 0, 0,
                  paramsScenario$DBVSBVEffs[["focal1"]]), ncol = 2)

  datmix <- datmix %>% arrange(neighbors, block)
  rownames(datmix) <- NULL

  if (is.null(datmix)){
    stop("simMix function don't take null data")
  }


  mix <- droplevels(datmix[datmix$stand=="mixed",])
  nMix <- length(unique(mix$name))


  ## Reformat
  d <- ncol(B_mix) # number of species
  p <- nrow(B_mix) # number of fixed effects in Bmix
  n_1 <- nrow(datmix)
  n_2 <- n_1 * d

  #check
  stopifnot(all(dim(B_mix) == c(p, d)))

  #Xmix
  listXs[["X_mix"]] <- model.matrix(~ 1 + block, datmix,
                                    contrasts.arg=list(block="contr.sum"))

  #check
  stopifnot(all(dim(listXs[["X_mix"]]) == c(n_1, p)))


  ## append species for Z
  datmix$species <- "wheat"

  ## list of matrix
  listZs[["Z_BV"]] <- mkZDBV(df = datmix, colD= "focal", colS = "species", species = "wheat",
                             skipUnusedCols = TRUE)

  listZs[["Z_SIGV_mix"]] <-  mkZSIGV(df = datmix,
                                     focal = "focal",
                                     species = "species",
                                     stand="stand",
                                     mkZForSpec = "wheat", weight = 0,
                                     skipUnusedCols = F, mixed = TRUE)

  listZs[["Z_DBV_x_SBV"]] <- mkZDBV.SBV(df=datmix, sep="_", name="name", species="species",Spec="both",skipUnusedCols = TRUE)

  listZs[["ZDl_mix"]] <- cbind(model.matrix(~1+neighbors, data=datmix, contrasts.arg = list("neighbors"="contr.sum"))[, 2],
                               model.matrix(~1+neighbors, data=datmix, contrasts.arg = list("neighbors"="contr.sum"))[, 2])

  ## select

  ##Error
  E = simE(n_1 = nrow(datmix), d = 2, paramsScenario, simulId)

  ##Check
  stopifnot(paramsScenario$panelGenos==names(explFacts$SIGVpanelGenos))
  stopifnot(paramsScenario$panelGenos==colnames(listZs[["Z_BV"]]))
  stopifnot(paramsScenario$panelGenos==colnames(listZs[["Z_SIGV_mix"]]))
  stopifnot(names(explFacts$SIGVpanelGenos)==colnames(listZs[["Z_SIGV_mix"]]))


  # add interactions
  explFacts[["U_DBV_x_SBV"]] <- explFacts[["U_DBV_x_SBV"]][colnames(listZs[["Z_DBV_x_SBV"]]), ]
  stopifnot(identical(as.character(colnames(listZs[["Z_DBV_x_SBV"]])),
                      as.character(rownames(explFacts[["U_DBV_x_SBV"]]))))

  explFacts[["Kmix"]] <- explFacts[["Kmix"]][colnames(listZs[["Z_DBV_x_SBV"]]),
                                             colnames(listZs[["Z_DBV_x_SBV"]])]

  listVCov[["Kmix"]] <- explFacts[["Kmix"]]
  stopifnot(identical(as.character(colnames(listZs[["Z_DBV_x_SBV"]])),
                      as.character(rownames(explFacts[["Kmix"]]))))

  # update Kmix and check correspondance with Z_DBV_x_SBV colnames
  Y_mix_1 <- listXs[["X_mix"]] %*% B_mix + listZs[["ZDl_mix"]]%*%U_p +
    listZs[["Z_BV"]]%*%explFacts[["U_DBV_SBV"]]  +
    listZs[["Z_DBV_x_SBV"]]%*%explFacts[["U_DBV_x_SBV"]] + E

  SISmix <- listZs[["Z_SIGV_mix"]]%*%explFacts$SIGVpanelGenos
  if(strategy =="inter_only"||strategy == "inter_only_ctrl"){

    Y_mix_2 <- Y_mix_1

  } else{

    Y_mix_2 <- matrix(NA, ncol = 2, nrow = nrow(Y_mix_1))
    for(i in 1: nrow(Y_mix_2)){
      for(j in 1: nrow(Y_mix_2)){
        Y_mix_2[i, 1] <- Y_mix_1[i, 1]+SISmix[i]
        Y_mix_2[j, 2] <- Y_mix_1[j, 2]
      }
    }
  }
  colnames(Y_mix_2) <- c("yield.wheat", "yield.pea")

  datmix <- cbind(datmix, Y_mix_2)

  ## re-write X_mix for infer
  listXs[["X_mix"]] <- model.matrix(~1+block+neighbors, data = datmix,
                                    contrasts=list("block"="contr.sum",
                                                   "neighbors"="contr.sum"))
  #remove Z_DBV_x_SBV for inference
  if(all(abs(colSums(as.matrix(explFacts[["U_DBV_x_SBV"]]))) < 1e-10)) {
    listZs[["Z_DBV_x_SBV"]] <- NULL
  }

  out <- list(data = datmix,
              Y = Y_mix_2,
              mixZs= listZs,
              mixXs = listXs,
              mixVCov = listVCov)

  return(out)
}


##'function to generate dataset
simulPhenos <- function(paramsScenario, explFacts,
                        outDesign, strategy, simulId){

  if (is.null(outDesign)){
    stop("SimulPhenos function don't take null outDesign")

  }
  # range check
  stopifnot(is.data.frame(outDesign))


  #create a list to store return data and matrix
  data <- list()
  listZs <- list()
  listXs <- list()
  listVCov <- list()


  #setting up dataset in pur stand
  if(strategy %in% c("sole_only","sole_inter_50","sole_inter_ctrl")){

    #set datPur based on sole_inter_50 or sole_only
    datPur <- droplevels(outDesign[outDesign$stand=="pure", ])

    tmp <- simPur(datPur, listZs, listXs,
                  listVCov,paramsScenario,
                  explFacts, strategy,simulId)

    datPur <- tmp$data

    if ("neighbors"%in%names(datPur)){
      datPur <- datPur[, -which(names(datPur) == "neighbors")]
    }

    purZs <- tmp$purZs
    purXs <- tmp$purXs
    purVCov <- tmp$purVCov
    y = tmp$y
    rm(tmp)
  }

  #pure stand only
  if(strategy=="sole_only"){
    data <- list(data =  list(datPur = datPur, y=y),
                 listZs = purZs, listXs = purXs, listVCov = purVCov)
  }

  #set datmix based on inter_only or sole_inter_50
  if (strategy %in% c("inter_only","sole_inter_50", "inter_only_ctrl", "sole_inter_ctrl")){
    datmix <- droplevels(outDesign[outDesign$stand=="mixed",])
  }


  #mix + pure&mix
  if(strategy %in% c("inter_only","sole_inter_50", "inter_only_ctrl","sole_inter_ctrl")){
    tmp <- simMix(datmix, listZs, listXs, listVCov,
                  paramsScenario,
                  explFacts, strategy, simulId)


    datmix <- tmp$data
    mixZs = tmp$mixZs
    mixXs = tmp$mixXs
    mixVCov = tmp$mixVCov
    Y = tmp$Y

    # #rm fke species
    datmix <- datmix[, -which(names(datmix) == "species")]

    #remove tmp
    rm(tmp)

    #set return data
    if(strategy=="inter_only" || strategy=="inter_only_ctrl"){
      data <- list(data =  list(datmix = datmix, Y = Y),
                   listZs=mixZs, listXs = mixXs, listVCov = mixVCov)
    }

    #set return data
    if(strategy=="sole_inter_50"||strategy=="sole_inter_ctrl"){
      data <- list(data =  list(datmix = datmix, datPur = datPur,
                                Y = Y, y = y),
                   listZs = c(purZs, mixZs),
                   listVCov = c(purVCov, mixVCov),
                   listXs = c(purXs, mixXs))}

  }

  return(data)
  invisible(gc()) #free up memory
}


#####                  FUNCTIONS FOR INFERENCE            #####################

##' Convert a Matrix or List of Matrices to Sparse Matrix Format
##'
##' This function converts either a single dense matrix or a list of dense matrices
##' to a specified sparse matrix format from the \pkg{Matrix} package (e.g., "dgCMatrix", "dgTMatrix").
##' If the matrix is already sparse, it is returned unchanged with a warning.
##'
##' @param Z A single numeric matrix, a sparse matrix, or a list of matrices.
##' @param type Character string specifying the sparse format to convert to.
##'   Default is `"dgCMatrix"`. See \code{\link[Matrix]{Matrix}} for possible formats.
##'
##' @return If \code{Z} is a matrix, returns a single sparse matrix of the specified type.
##'   If \code{Z} is a list, returns a list of sparse matrices of the specified type,
##'   preserving the element names if present.
##'
##' @author Jemay Salomon
asSparseMatrix <- function(Z, type = "dgCMatrix") {

  require(Matrix)
  isSparse <- function(x) inherits(x, "sparseMatrix")
  if (!(is.matrix(Z) || is.list(Z) || isSparse(Z))) {
    stop("Z must be a matrix, a sparse matrix, or a list of matrices")
  }
  if (isSparse(Z)) {
    warning("Input is already a sparse matrix — returning as is.")
    return(Z)
  }
  if (is.matrix(Z)) {
    return(as(Z, type))
  }
  if (is.list(Z)) {
    return(lapply(Z, function(spm_mat) {
      if (!(is.matrix(spm_mat) || isSparse(spm_mat)))
        stop("All elements of the list must be matrices or sparse matrices")
      if (isSparse(spm_mat)) {
        warning("is are already sparse — returning it unchanged.")
        return(spm_mat)
      }
      as(spm_mat, type)
    }))
  }
}

#' Validate design matrices and variance-covariance list
#'
#' @description
#' This function validates that the dimensions and names of the
#' provided Z matrices match each other and the variance-covariance
#' matrix `K` in `listVCov`. It stops with an error if any mismatch
#' is found.
#'
#' @param listZs A named list of design matrices containing:
#'   \code{Z_DBV}, \code{Z_SIGV}, \code{Z_BV}, \code{Z_SIGV_mix}.
#' @param listVCov A named list containing at least the matrix \code{K}.
#' @param Y A list containing \code{pure} and \code{mix} vectors (for row counts).
#'
#' @author Jemay Salomon
#'
#' @return Invisibly returns TRUE if all checks pass, otherwise stops with an error.
#'
validateInput <- function(listZs, listVCov, listXs, Y) {

  # Extract K matrix
  K <- listVCov[["K"]]

  # Helper: safe name compare
  sn <- function(x, y) {
    identical(as.character(x), as.character(y))
  }

  # 1. Basic col/rowname matches
  if("y"%in%names(Y)){
    stopifnot(sn(colnames(listZs[["Z_DBV"]]), rownames(K)))
    stopifnot(sn(colnames(listZs[["Z_SIGV"]]), colnames(K)))
    stopifnot(ncol(listZs[["Z_SIGV"]]) == ncol(listZs[["Z_DBV"]]))
    stopifnot(nrow(listZs[["Z_SIGV"]]) == nrow(listZs[["Z_DBV"]]))
    stopifnot(nrow(listZs[["Z_SIGV"]]) == length(Y$y))}

  if("Z_SIGV_mix"%in%names(listZs))
    stopifnot(sn(colnames(listZs[["Z_SIGV_mix"]]), rownames(K)))

  # 3. K name consistency across all Z matrices
  for (mat_name in c("Z_BV", "Z_SIGV_mix")){
    if(mat_name%in%names(listZs)){
      stopifnot(sn(colnames(listZs[[mat_name]]), colnames(K)))
      stopifnot(sn(colnames(listZs[[mat_name]]), rownames(K)))
    }
  }

  stopifnot(nrow(listZs[["Z_BV"]]) == nrow(listZs[["Z_SIGV_mix"]]))
  stopifnot(nrow(listZs[["Z_BV"]]) == nrow(Y$Y))
  stopifnot(nrow(listXs[["X_mix"]]) == nrow(Y$Y))

  invisible(TRUE)
}

# Define formula based on strategy
formula <- function(strategy) {
  if (strategy == "sole_only") {
    return(yield ~ 1 + block)
  } else if (strategy == "inter_only" || strategy=="inter_only_ctrl") {
    return(Y ~ 1 + block)
  } else if (strategy == "sole_inter_50"||strategy=="sole_inter_ctrl") {
    return(list(yield ~ 1 + block, Y ~ 1 + block))
  } else {
    stop("Unknown strategy")
  }
}



## @param capture output from \code{capture.output(nlminb(..., control=list("trace"=1)))}
## from plantmix (upcoming R package)
traceFromNlminb <- function(capture) {
  stopifnot(is(capture, "character"))
  out <- NULL
  idx <- grep("^[ 0-9]+:", capture)
  if (length(idx) > 0) {
    out <- capture[idx]
    out <- strcapture(
      "^(\\s+[0-9]+):\\s+([.0-9]+):",
      out,
      data.frame(
        iter = numeric(),
        objfn = numeric()
      )
    )
  }
  return(out)
}

setMap <- function(listData, listPars, inMap=NULL){
  defltMap <- list()
  if(listData$model=="MMTMB"){
    if(length(listData$y)==0 || is.null(listData$y)){
      defltMap[["beta"]] = rep(NA, length(listPars[["beta"]]))
      defltMap[["SIGV"]] = rep(NA, length(listPars[["SIGV"]]))
      defltMap[["log_sd_e"]] = NA
      defltMap[["log_sd_sigv"]] = NA
    }
    if(length(listData$y_t)==0 || is.null(listData$y_t)){
      defltMap[["beta_t"]] = rep(NA, length(listPars[["beta_t"]]))
      defltMap[["log_sd_e_t"]] = NA
    }
    if(ncol(listData$listZs[["Z_DBV_x_SBV"]])<=1 && nrow(listData$listZs[["Z_DBV_x_SBV"]])<=1){
      defltMap[["DBV_x_SBV"]]= matrix(NA, nrow = ncol(listData$listZs[["Z_DBV_x_SBV"]]), ncol(listData[["Y"]]))
      defltMap[["log_sd_DBV_x_SBV"]] = rep(NA, ncol(listData[["Y"]]))
      defltMap[["unconstr_cor_DBV_x_SBV"]] = NA
    }
    if(ncol(listData$listBs[["Bs_SC_f"]])<=1 && nrow(listData$listBs[["Bs_SC_f"]])<=1){
      defltMap[["cs_SC_f"]] = rep(NA, length(listPars[["cs_SC_f"]]))
      defltMap[["log_lambda_u_sc_f"]] = rep(NA, length(listPars[["log_lambda_u_sc_f"]]))
      defltMap[["log_lambda_v_sc_f"]] = rep(NA, length(listPars[["log_lambda_v_sc_f"]]))
    }
    if(ncol(listData$listBs[["Bs_SC_t"]])<=1 && nrow(listData$listBs[["Bs_SC_t"]])<=1){
      defltMap[["cs_SC_t"]] = rep(NA, length(listPars[["cs_SC_t"]]))
      defltMap[["log_lambda_u_sc_t"]] = rep(NA, length(listPars[["log_lambda_u_sc_t"]]))
      defltMap[["log_lambda_v_sc_t"]] = rep(NA, length(listPars[["log_lambda_v_sc_t"]]))
    }
    if(ncol(listData$listBs[["Bs_IC_f"]])<=1 && nrow(listData$listBs[["Bs_IC_f"]])<=1){
      defltMap[["cs_IC_f"]] = rep(NA, length(listPars[["cs_IC_f"]]))
      defltMap[["log_lambda_u_ic_f"]] = rep(NA, length(listPars[["log_lambda_u_ic_f"]]))
      defltMap[["log_lambda_v_ic_f"]] = rep(NA, length(listPars[["log_lambda_v_ic_f"]]))
    }
    if(ncol(listData$listBs[["Bs_IC_t"]])<=1 && nrow(listData$listBs[["Bs_IC_t"]])<=1){
      defltMap[["cs_IC_t"]] = rep(NA, length(listPars[["cs_IC_t"]]))
      defltMap[["log_lambda_u_ic_t"]] = rep(NA, length(listPars[["log_lambda_u_ic_t"]]))
      defltMap[["log_lambda_v_ic_t"]] = rep(NA, length(listPars[["log_lambda_v_ic_t"]]))
    }
  }
  if(listData$model=="lmm"){
    if(ncol(listData[["Bs"]])<=1 && nrow(listData[["Bs"]])<=1){
      defltMap[["theta"]] = rep(NA, length(listPars[["theta"]]))
      defltMap[["log_lambda_u"]] = rep(NA, length(listPars[["log_lambda_u"]]))
      defltMap[["log_lambda_v"]] = rep(NA, length(listPars[["log_lambda_v"]]))
    }
  }
  defltMap <- lapply(defltMap, as.factor)
  if(is.null(inMap)){
    inMap <- defltMap
  } else{
    stopifnot(! is.null(names(inMap)))
    inMap[intersect(names(inMap), names(defltMap))] <- NULL
    inMap <- c(defltMap, inMap)
  }
  return(inMap)
}

##' Fit DBV-SBV models for interspecific mixtures. Only needed parameters to reproduce the paper analyses is described here. Other parameters can be used as default
##'
##' Fits DBV-SBV models for interspecific mixtures.
##' @param Y named list of response variables (a matrix Y and, optionally, a vector y and a vector y_t)
##' @param listX named list of design matrices for the fixed-effect explanatory factors (a matrix X_mix and, optionally, a matrix X and a matrix X_t)
##' @param listZ named list of design matrices of random-effect explanatory factors
##' @param listVCov named list of square, symmetric matrices used, after re-scaling, as variance-covariance matrices for the random-effect factors; named "K" for the BVs (DBVs and SBVs) and "Kmix" for the DBVxSBVs; these matrices must have dimnames (both rows and columns) coherent with the column names of the design matrices in \code{listZ}
##' @param REML logical
##' @param control named list of options (for advanced users)
##' @param verbose verbosity level (True or FALSE)
##' @param verbose silent(True or FALSE)
##' @return list
##' @author Jemay Salomon
##' @export
MMTMB <- function(Y, listXs, listZs, listVCov, listBs = NULL,
                  listPus = NULL, listPvs = NULL,
                  dllID = NULL,verbose=F,
                  silent = T, REML=TRUE,
                  Bayesian=FALSE,
                  sampling=NULL,
                  inferWithK=NULL,
                  sdrep=F,coef.fixed = F,
                  coef.rand = F,
                  control = NULL){

  # initialize
  if(is.null(listBs))
    listBs = list()
  if(is.null(listPus))
    listPus = list()
  if(is.null(listPvs))
    listPvs = list()

  requireNamespace(package="TMB")
  validateInput(listZs, listVCov, listXs, Y)

  # check bayesian if TRUE
  if(Bayesian){
    requireNamespace(package="tmbstan")
    requireNamespace(package="rstan")
    stopifnot(!is.null(sampling),
              is.list(sampling))}

  vcat <- function(...) if (verbose) message(...)
  out <- NULL

  vcat("checking input...")
  dllID = "MMTMB"
  try(dyn.unload(file.path(srcDir, TMB::dynlib(dllID))), silent = TRUE)
  if (!file.exists(file.path(srcDir, paste0(dllID, .Platform$dynlib.ext)))) {
    if(Sys.info()[["sysname"]] == "Windows"){
      FLAGSX = '-g -O1'
      TMB::compile(
        file = file.path(srcDir, paste0(dllID, ".cpp")),
        flags = FLAGSX,
        DLLFLAGS = "")
    } else {
      FLAGSX = '-g -O0 -Wall'
      TMB::compile(
        file=file.path(srcDir, paste0(dllID, ".cpp")),
        flags = FLAGSX)
    }
  }
  dyn.load(file.path(srcDir, TMB::dynlib(dllID)))
  out$dllID <- dllID

  #tester conditional
  if(!"y_t" %in% names(Y))
    Y[["y_t"]] <- numeric(0)
  if(!"X_t" %in% names(listXs))
    listXs[["X_t"]] <- matrix(0,1,1)

  #add focal conditional on y
  if(!"y" %in% names(Y))
    Y[["y"]] <- numeric(0)
  if(!"X" %in% names(listXs))
    listXs[["X"]]<- matrix(0,1,1)
  if(!"Z_DBV" %in% names(listZs))
    listZs[["Z_DBV"]] <- matrix(0,1,1)
  if(!"Z_SIGV" %in% names(listZs))
    listZs[["Z_SIGV"]] <- matrix(0,1,1)
  if(!"Z_SIGV_mix" %in% names(listZs))
    listZs[["Z_SIGV_mix"]] <- matrix(0,1,1)
  if(!"X" %in% names(listXs))
    listXs[["X"]]<- matrix(0,1,1)

  #add interactions conditional on DBV_x_SBV
  if(!"Z_DBV_x_SBV" %in% names(listZs))
    listZs[["Z_DBV_x_SBV"]] <- matrix(0,1,1)
  if(!"Kmix" %in% names(listVCov))
    listVCov[["Kmix"]] <- matrix(0,1,1)

  #SC
  if(!"Bs_SC_f" %in% names(listBs))
    listBs[["Bs_SC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Bs_SC_t" %in% names(listBs))
    listBs[["Bs_SC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pu_SC_f" %in% names(listPus))
    listPus[["Pu_SC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pu_SC_t" %in% names(listPus))
    listPus[["Pu_SC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pv_SC_f" %in% names(listPvs))
    listPvs[["Pv_SC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pv_SC_t" %in% names(listPvs))
    listPvs[["Pv_SC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")

  #IC
  if(!"Bs_IC_f" %in% names(listBs))
    listBs[["Bs_IC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Bs_IC_t" %in% names(listBs))
    listBs[["Bs_IC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pu_IC_f" %in% names(listPus))
    listPus[["Pu_IC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pu_IC_t" %in% names(listPus))
    listPus[["Pu_IC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pv_IC_f" %in% names(listPvs))
    listPvs[["Pv_IC_f"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(!"Pv_IC_t" %in% names(listPvs))
    listPvs[["Pv_IC_t"]] <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")

  # beta start
  if (is.null(control[["beta_start_t"]]))
    beta_t <- rep(0, ncol(listXs[["X_t"]]))
  else
    beta_t <- control[["beta_start_t"]]

  if (is.null(control[["beta_start_f"]]))
    beta_f <- rep(0, ncol(listXs[["X"]]))
  else
    beta_f <- control[["beta_start_f"]]

  if(is.null(control[["sim_with_spats"]]))
    control[["sim_with_spats"]] <- 0

  d = ncol(Y$Y)

  listZs=asSparseMatrix(listZs, "dgCMatrix")

  vcat("constructing TMB data and parameters...")

  ## list of data variables and parameter variables
  listData <- list(model = "MMTMB",
                   Y = Y$Y,
                   X_mix = listXs[["X_mix"]],
                   listZs = listZs,
                   listBs = listBs,
                   listPus = listPus,
                   listPvs = listPvs,
                   sim_with_spats = control[["sim_with_spats"]],
                   Id_mix=diag(nrow(listXs[["X_mix"]])),

                   y = Y$y,
                   y_t = Y$y_t,
                   X = listXs[["X"]],
                   X_t = listXs[["X_t"]],
                   Id = diag(nrow(listXs[["X"]])),
                   Id_t = diag(nrow(listXs[["X_t"]])),
                   A = listVCov[["K"]],
                   Amix= listVCov[["Kmix"]])

  listPars <- list(B = matrix(c(rep(0, 2*ncol(listXs[["X_mix"]]))),byrow=T, ncol=2),
                   BV = matrix(0, nrow = ncol(listZs[["Z_BV"]]), ncol = d),
                   log_sd_BV = rep(log(1), d),
                   log_sd_E = rep(log(1), d),
                   unconstr_cor_BV = rep(0, 1),
                   unconstr_cor_E = rep(0, 1),

                   beta = beta_f,
                   beta_t = beta_t,
                   SIGV = rep(0, ncol(listZs[["Z_SIGV"]])),
                   log_sd_sigv=log(1),
                   log_sd_e=log(1),
                   log_sd_e_t=log(1),
                   DBV_x_SBV = matrix(0, nrow = ncol(listZs[["Z_DBV_x_SBV"]]), ncol = d),
                   log_sd_DBV_x_SBV = rep(log(1), d),
                   unconstr_cor_DBV_x_SBV= rep(0, 1),

                   # coefs
                   cs_SC_f = rep(0, ncol(listBs[["Bs_SC_f"]])),
                   cs_SC_t = rep(0, ncol(listBs[["Bs_SC_t"]])),

                   cs_IC_f = rep(0, ncol(listBs[["Bs_IC_f"]])),
                   cs_IC_t = rep(0, ncol(listBs[["Bs_IC_t"]])),

                   # params
                   log_lambda_u_sc_f = log(1),
                   log_lambda_v_sc_f = log(1),
                   log_lambda_u_ic_f = log(1),
                   log_lambda_v_ic_f = log(1),

                   log_lambda_u_sc_t = log(1),
                   log_lambda_v_sc_t = log(1),
                   log_lambda_u_ic_t = log(1),
                   log_lambda_v_ic_t = log(1))

  rands <- c("BV")
  if(length(listData$y)>0)
    rands <-  union(rands, c("SIGV"))
  if(ncol(listData$listZs[["Z_DBV_x_SBV"]])>1)
    rands <- union(rands, c("DBV_x_SBV"))
  if(ncol(listData$listBs[["Bs_SC_f"]])>1)
    rands <- union(rands, c("cs_SC_f"))
  if(ncol(listData$listBs[["Bs_SC_t"]])>1)
    rands <- union(rands, c("cs_SC_t"))
  if(ncol(listData$listBs[["Bs_IC_f"]])>1)
    rands <- union(rands, c("cs_IC_f"))
  if(ncol(listData$listBs[["Bs_IC_t"]])>1)
    rands <- union(rands, c("cs_IC_t"))

  if(REML){
    if(length(listData$y)>0)
      rands <-  union(rands, c("beta"))
    if(length(listData$y_t)>0)
      rands <-  union(rands, c("beta_t"))
    rands <- union(rands, c("B"))
  }

  map <- setMap(listData, listPars, inMap = control[["map"]])
  vcat("Construct objective functions...")
  t0 <- proc.time()
  obj <- TMB::MakeADFun(data=listData,
                        parameters=listPars,
                        random=rands,
                        ## hessian=TRUE,
                        fit = TRUE,
                        map=map,
                        DLL=basename(dllID),
                        silent=silent)
  t1 <- proc.time()
  out$tmb_fit_time <- (t1-t0)[1]

  vcat("Optimization...")
  t0 <- proc.time()
  capture <- capture.output(
    fit <- nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = NULL,
      control = list("trace" = 1)
    )
  )
  t1 <- proc.time()
  out$nlminb_optim_time <- (t1-t0)[1]
  vcat("Preparing output...")
  fit$parfull <- obj$env$last.par.best
  out$obj <- obj
  out$fit <- fit
  out$trace <- traceFromNlminb(capture)
  out$map = map

  if(sdrep)
    out$sdrep <- sdreport(obj)

  if(coef.rand){
    randsZi <- c("BV")
    if(length(listData$y) > 0)
      randsZi <- c(randsZi, "SIGV")
    if(ncol(listZs[["Z_DBV_x_SBV"]]) > 1)
      randsZi <- c(randsZi, "DBV_x_SBV")
    out$coef.random <- tmbExtract(out, dllID = dllID, params = randsZi,
                                  reNames = NULL, paramsType = "random", path = srcDir)
  }

  if(coef.fixed){
    betaZi <- c("B")
    if(length(listData$y) > 0) betaZi <- c(betaZi, "beta")
    if(length(listData$y_t) > 0) betaZi <- c(betaZi, "beta_t")

    if(REML)
      out$coef.fixed <- tmbExtract(out, dllID = dllID, params = betaZi,
                                   reNames = NULL, paramsType = "random", path = srcDir)
    else
      out$coef.fixed <- tmbExtract(out, dllID = dllID, params = betaZi,
                                   reNames = NULL, paramsType = "paramsTmb", path = srcDir)
  }

  if(Bayesian){
    vcat("Bayesian sampling...")
    rstan::rstan_options(auto_write = FALSE)
    sink("/dev/null")
    mcmc <- tmbstan(obj,
                    chains = sampling[["chains"]],
                    iter = sampling[["iter"]],
                    control = sampling[["control"]],
                    init = sampling[["init"]],
                    open_progress = FALSE,
                    verbose = FALSE,
                    silent = TRUE,
                    show_messages = FALSE)
    sink()
    out$mcmc <- mcmc
  }

  ## add attribute to out
  class(out) <- "MMTMB"
  return(out)
  # invisible(gc()) #free up memory
}


##' simple linear mixed model. Only needed parameters to reproduce the paper analyses is described here. Other parameters can be used as default
##'
##' @param Y named list of response variable (only y is needed)
##' @param listX named list of design matrices for the fixed-effect explanatory factors (only X is needed)
##' @param listZ named list of design matrices of random-effect explanatory factors (only one Z is needed)
##' @param listVCov named list of square, symmetric matrices used, after re-scaling, as variance-covariance matrices for the random-effect factors; only "K" for the BVs in solecrop. these matrices must have dimnames (both rows and columns) coherent with the column names of the design matrices in \code{listZ}
##' @param REML logical
##' @param control named list of options (for advanced users)
##' @param verbose verbosity level (True or FALSE)
##' @param verbose silent(True or FALSE)
##' @return list
##' @author Jemay Salomon
##' @export
lmmTMB <- function(Y, listXs,listZs, listVCov,
                   dllID=NULL,
                   verbose = F, silent = T,
                   REML=TRUE,Bayesian=FALSE,
                   sampling=NULL,
                   inferWithK = TRUE,
                   log_lambda=NULL,
                   coef.fixed = T,
                   Bs=NULL,
                   Pu=NULL, Pv = NULL,
                   control = NULL){

  requireNamespace(package="TMB")

  # check bayesian control
  if(Bayesian){
    requireNamespace(package="tmbstan")
    requireNamespace(package="rstan")
    stopifnot(!is.null(sampling),
              is.list(sampling))}

  vcat <- function(...) if (verbose) message(...)

  out <- NULL

  #Compilation condions
  dllID = "MMTMB"
  try(dyn.unload(file.path(srcDir, TMB::dynlib(dllID))), silent = TRUE)
  if (!file.exists(file.path(srcDir, paste0(dllID, .Platform$dynlib.ext)))) {
    if(Sys.info()[["sysname"]] == "Windows"){
      FLAGSX = '-g -O1'
      TMB::compile(
        file = file.path(srcDir, paste0(dllID, ".cpp")),
        flags = FLAGSX,
        DLLFLAGS = "")
    } else {
      FLAGSX = '-g -O0 -Wall'
      TMB::compile(
        file=file.path(srcDir, paste0(dllID, ".cpp")),
        flags = FLAGSX)
    }
  }
  dyn.load(file.path(srcDir, TMB::dynlib(dllID)))
  out$dllID <- dllID

  # internal params
  q = ncol(listZs[["Z_DBV"]])

  listZs=asSparseMatrix(listZs, "dgCMatrix")

  #SC
  if(is.null(Bs))
    Bs <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(is.null(Pu))
    Pu <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")
  if(is.null(Pv))
    Pv <- as(as(as(matrix(0,1,1), "dMatrix"), "generalMatrix"), "TsparseMatrix")

  if(is.null(control))
    control[["sim_with_spats"]] <- 0

  ## list of data
  listData <- list(model= "lmm",
                   Bs = Bs,
                   Pu = Pu,
                   Pv = Pv,
                   sim_with_spats = control[["sim_with_spats"]],
                   listZs=listZs,
                   y=Y[["y"]],
                   X = listXs[["X"]],
                   A = listVCov[["K"]])

  ## list of parameters
  listPars <- list(beta=rep(0, ncol(listXs[["X"]])),
                   u=rep(0, q),
                   log_lambda_u = log(1),
                   log_lambda_v = log(1),
                   theta = rep(0, ncol(listData[["Bs"]])),
                   log_sd_u=log(1),
                   log_sd_e=log(1))

  rands <- c("u")
  if(REML)
    rands <- union(rands, c("beta"))

  if(ncol(listData$Bs)>1)
    rands <- union(rands, c("theta"))


  if(!inferWithK)
    listData[["model"]] = "uv_lmm_Id" # dont use for now


  vcat("Construct objective functions")
  cat("model: ",listData[["model"]], "\n")
  t0 <- proc.time()
  obj <- TMB::MakeADFun(data=listData,
                        parameters=listPars,
                        random=rands,
                        ## hessian=TRUE,
                        DLL=basename(dllID),
                        silent=silent)
  t1 <- proc.time()

  out$tmb_fit_time <- (t1-t0)[1]

  vcat("Optimization...")
  t0 <- proc.time()
  fit <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, hessian=NULL)
  t1 <- proc.time()
  out$nlimb_optim_time <- (t1-t0)[1]

  vcat("Preparing output...")
  out$obj <- obj
  fit$parfull <- obj$env$last.par.best
  out$fit <- fit
  out$coef.random  <-tmbExtract(out,
                                dllID = dllID,
                                params = c("u"),
                                reNames = NULL,
                                paramsType = "random",
                                path = srcDir)
  if(coef.fixed){
    if(REML)
      out$coef.fixed <-tmbExtract(out,
                                  dllID = dllID,
                                  params = c("beta"),
                                  reNames = NULL,
                                  paramsType = "random",
                                  path = srcDir)
    else
      out$coef.fixed <-tmbExtract(out,
                                  dllID = dllID,
                                  params = "beta",
                                  reNames = NULL,
                                  paramsType = "paramsTmb",
                                  path = srcDir)
  }

  if (Bayesian) {
    vcat("Bayesian sampling...")
    rstan::rstan_options(auto_write = FALSE)
    sink("/dev/null")
    mcmc <- tmbstan(obj,
                    chains = sampling$chains,
                    iter = sampling$iter,
                    control = sampling$control,
                    init = sampling$init,
                    open_progress = FALSE,
                    verbose = FALSE,
                    silent = TRUE,
                    show_messages = FALSE)
    sink()
    out$mcmc <- mcmc
  }

  class(out) <- "MMTMB"
  return(out)
}


##########     HELPER FUNCTIONS TO RENAME SIMULATION PARAMETERS     ############
##########                   IN INFERENCE                    ############

## Renames parameters
namesParams <- function(strategy, DBVxSBV) {
  if (strategy == "sole_only") {
    return(c("mu_w_pur", "block1_w_pur", "sigma2.cBV", "sigma2_err_w"))

  } else if (strategy == "inter_only" || strategy == "inter_only_ctrl") {
    NamesParams <- c("mu_w_mix", "block1_w_mix","neighbors1",
                     "mu_p_mix", "block1_p_mix", "focal1",
                     "sigma2.DBV", "sigma2.SBV","sigma2_err_w_mix",
                     "sigma2_err_p_mix")
    if(isTRUE(DBVxSBV))
      NamesParams <- c(NamesParams, "sigma2.DBV.SBV.wp","sigma2.DBV.SBV.pw")
    NamesParams <- c( NamesParams, "Cor_BV", "Cor_E")
    if(isTRUE(DBVxSBV))
      NamesParams <- c(NamesParams, "Cor_DBV_x_SBV")

    return(NamesParams)
  } else if (strategy == "sole_inter_50"||strategy=="sole_inter_ctrl") {
    NamesParams <- c("mu_w_mix", "block1_w_mix", "neighbors1",
                     "mu_p_mix","block1_p_mix", "focal1","mu_w_pur", "block1_w_pur",
                     "sigma2.DBV", "sigma2.SBV","sigma2.SIGV", "sigma2_err_w_mix",
                     "sigma2_err_p_mix","sigma2_err_w_pur")
    if(isTRUE(DBVxSBV))
      NamesParams <- c(NamesParams, "sigma2.DBV.SBV.wp","sigma2.DBV.SBV.pw")

    NamesParams <- c( NamesParams, "Cor_BV", "Cor_E")

    if(isTRUE(DBVxSBV))
      NamesParams <- c(NamesParams, "Cor_DBV_x_SBV")

    return(NamesParams)

  } else {
    stop("Unknown strategy")
  }
}

# Rename beta from tmbExtract based on strategy
betaParams <- function(strategy) {
  if (strategy == "sole_only") {
    return(c("mu_w_pur", "block1_w_pur"))

  } else if (strategy == "inter_only"|| strategy == "inter_only_ctrl") {
    return(c("mu_w_mix", "block1_w_mix","neighbors1", "mu_p_mix",
             "block1_p_mix", "focal1"))

  } else if (strategy == "sole_inter_50"||strategy=="sole_inter_ctrl") {
    return(c("mu_w_mix", "block1_w_mix","neighbors1",
             "mu_p_mix","block1_p_mix","focal1", "mu_w_pur", "block1_w_pur"))

  } else {
    stop("Unknown strategy")
  }
}


#variance parameters to extract
# get interactions
varToExtract <- function(strategy,DBVxSBV){
  if (strategy == "sole_only") {
    return(c("log_sd_u", "log_sd_e"))
  } else if (strategy == "inter_only"||strategy == "inter_only_ctrl") {
    vars <- c("log_sd_BV", "log_sd_E")
    if(isTRUE(DBVxSBV))
      vars <- c(vars, "log_sd_DBV_x_SBV")
    return(vars)

  } else if (strategy == "sole_inter_50"||strategy=="sole_inter_ctrl") {

    vars <- c("log_sd_BV", "log_sd_sigv", "log_sd_E", "log_sd_e")
    if(isTRUE(DBVxSBV))vars <- c(vars, "log_sd_DBV_x_SBV")

    return(vars)

  } else {
    stop("Unknown strategy")
  }
}


corToExtract <- function(DBVxSBV){
  corParams <- c("Cor_BV", "Cor_E")
  if(isTRUE(DBVxSBV))
    corParams <- c(corParams, "Cor_DBV_x_SBV")
  return(corParams)
}


##get_true_parameters from ScenarioParams
# TODO: get true interactions params
get_true_params <- function(strategy, truth, DBVxSBV) {

  if (strategy == "sole_only") {
    betaTr <- truth$paramsScenario$block1SpeciesEffs[,"pure", "wheat"]
    varTr <- c(truth$paramsScenario$sigma2.cBVw, truth$paramsScenario$sigma2_err_w)
    true <- c(betaTr, varTr)

  } else if (strategy == "inter_only"||strategy == "inter_only_ctrl") {

    betaTr <- c(truth$paramsScenario$block1SpeciesEffs[,"mixed", "wheat"],
                truth$paramsScenario$DBVSBVEffs$neighbors1,
                truth$paramsScenario$block1SpeciesEffs[,"mixed", "pea"],
                truth$paramsScenario$DBVSBVEffs$focal1)

    varTr <- c(truth$paramsScenario$sigma2.DBVw, truth$paramsScenario$sigma2.SBVw,
               truth$paramsScenario$sigma2_err_w_mix, truth$paramsScenario$sigma2_err_p_mix)

    if(isTRUE(DBVxSBV))
      varTr <- c(varTr, truth$paramsScenario$sigma2.DBV.SBV.wp,truth$paramsScenario$sigma2.DBV.SBV.pw)

    corTr <- c(truth$paramsScenario$cor_DBVw_SBVw, truth$paramsScenario$cor_E)
    if(isTRUE(DBVxSBV))
      corTr  <- c(corTr, truth$paramsScenario$cor_DBV_x_SBV)

    true <- c(betaTr, varTr, corTr)

  } else {
    betaTr <- c(truth$paramsScenario$block1SpeciesEffs[,"mixed", "wheat"],
                truth$paramsScenario$DBVSBVEffs$neighbors1,
                truth$paramsScenario$block1SpeciesEffs[,"mixed", "pea"],
                truth$paramsScenario$DBVSBVEffs$focal1,
                truth$paramsScenario$block1SpeciesEffs[,"pure", "wheat"])

    varTr <- c(truth$paramsScenario$sigma2.DBVw,
               truth$paramsScenario$sigma2.SBVw,
               truth$paramsScenario$sigma2.SIGVw,
               truth$paramsScenario$sigma2_err_w_mix,
               truth$paramsScenario$sigma2_err_p_mix,
               truth$paramsScenario$sigma2_err_w)

    if(isTRUE(DBVxSBV))
      varTr <- c(varTr, truth$paramsScenario$sigma2.DBV.SBV.wp,truth$paramsScenario$sigma2.DBV.SBV.pw)

    corTr <- c(truth$paramsScenario$cor_DBVw_SBVw,truth$paramsScenario$cor_E)
    if(isTRUE(DBVxSBV)) corTr  <- c(corTr, truth$paramsScenario$cor_DBV_x_SBV)

    true <- c(betaTr, varTr, corTr)
  }
  return(true)
}

#get estimated_values from TMB
# get interactions
get_estimated_params <- function(strategy, out,
                                 names_of_parameters,
                                 dllIDPath,
                                 DBVxSBV) {

  betaEs <- out$coef.fixed

  varEs <- tmbExtract(out,
                      dllID = "MMTMB",
                      params = varToExtract(strategy,
                                            DBVxSBV),
                      paramsType = "variance",
                      path = dllIDPath)

  if (!strategy %in% c("sole_only")) {
    corEs <- tmbExtract(out,
                        dllID = "MMTMB",
                        params = corToExtract(DBVxSBV),
                        paramsType = "correlation",path = dllIDPath)

    estimated <- c(betaEs, varEs, corEs)
  } else {
    estimated <- c(betaEs, varEs)
  }

  return(estimated)
}

#get parameters value
get_params_value <- function(strategy, out,
                             names_of_parameters,
                             truth, dllIDPath, DBVxSBV){

  #get estimated and true of parameters values
  estimated <- get_estimated_params(strategy, out,
                                    names_of_parameters=betaParams(strategy),
                                    dllIDPath,
                                    DBVxSBV)
  true <- get_true_params(strategy, truth, DBVxSBV)

  params <- cbind("estimated" = round(estimated, 2),
                  "true" = round(true, 2))
  #rownames parameters
  rownames(params) <- namesParams(strategy, DBVxSBV)

  return(params)
}

get_genetic_values <- function(strategy, truth, out,
                               dllIDPath, DBVxSBV) {

  # Load the DLL using dyn.load
  if(!is.null(dllIDPath))
    dyn.load(file.path(dllIDPath, TMB::dynlib(out$dllID)))
  else
    dyn.load(TMB::dynlib(out$dllID))

  DBVSBVSIGV <- list(DBV = NULL, SBV = NULL,
                     SIGV = NULL, DBVxSBV_w=NULL,
                     DBVxSBV_p=NULL)  # Initialized with NULL values

  dbv <- truth$explFacts$U_DBV_SBV[, 1]
  sigv <- truth$explFacts$SIGVpanelGenos
  sbv <- truth$explFacts$U_DBV_SBV[, 2]

  if (strategy == "sole_only") {
    DBVSBVSIGV$DBV <- cbind("true" = dbv,
                            "estimated" = tmbExtract(out,
                                                     dllID = "MMTMB",
                                                     params = c("u"),
                                                     reNames = NULL, paramsType = "random",
                                                     path = dllIDPath))


    DBVSBVSIGV$SBV <- cbind("true" = sbv,
                            "estimated" = rep(NA, length(sbv)))


    DBVSBVSIGV$SIGV <- cbind("true" = sigv,
                             "estimated" = rep(NA, length(sbv)))


  } else if (strategy == "inter_only"||strategy == "inter_only_ctrl") {
    DBVSBVSIGV$DBV <- cbind("true" = dbv,
                            "estimated" = tmbExtract(out,
                                                     dllID = "MMTMB",
                                                     params = c("BV"),
                                                     reNames = NULL,
                                                     paramsType = "random",
                                                     path = dllIDPath)[1:length(dbv)])

    DBVSBVSIGV$SBV <- cbind("true" = sbv,
                            "estimated" = tmbExtract(out,
                                                     dllID = "MMTMB",
                                                     params = c("BV"),
                                                     reNames = NULL,
                                                     paramsType = "random",
                                                     path = dllIDPath)[-(1:length(dbv))])

    DBVSBVSIGV$SIGV <- cbind("true" = sigv,
                             "estimated" = rep(NA, length(sigv)))
    if(isTRUE(DBVxSBV)){

      dbv.sbv.w <- truth$explFacts$U_DBV_x_SBV[, 1]
      dbv.sbv.p <- truth$explFacts$U_DBV_x_SBV[, 2]

      DBVSBVSIGV$DBVxSBV_w <- cbind("true" = dbv.sbv.w,
                                    "estimated" = tmbExtract(out,
                                                             dllID = "MMTMB",
                                                             params = c("DBV_x_SBV"),
                                                             reNames = NULL,
                                                             paramsType = "random",
                                                             path = dllIDPath)[1:length(dbv.sbv.w)])
      DBVSBVSIGV$DBVxSBV_p <- cbind("true" = dbv.sbv.p,
                                    "estimated" = tmbExtract(out,
                                                             dllID = "MMTMB",
                                                             params = c("DBV_x_SBV"),
                                                             reNames = NULL,
                                                             paramsType = "random",
                                                             path = dllIDPath)[-(1:length(dbv.sbv.p))])
    }

  } else {
    DBVSBVSIGV$DBV <- cbind("true" = dbv,
                            "estimated" = tmbExtract(out,
                                                     dllID = "MMTMB",
                                                     params = c("BV"),
                                                     reNames = NULL,
                                                     paramsType = "random",
                                                     path = dllIDPath)[1:length(dbv)])

    DBVSBVSIGV$SBV <- cbind("true" = sbv,
                            "estimated" = tmbExtract(out,
                                                     dllID = "MMTMB",
                                                     params = c("BV"),
                                                     reNames = NULL,
                                                     paramsType = "random",
                                                     path = dllIDPath)[-(1:length(sbv))])
    DBVSBVSIGV$SIGV <- cbind("true" = sigv,
                             "estimated" = tmbExtract(out,
                                                      dllID = "MMTMB",
                                                      params = c("SIGV"),
                                                      reNames = NULL, paramsType = "random",
                                                      path = dllIDPath))

    if(isTRUE(DBVxSBV)){

      dbv.sbv.w <- truth$explFacts$U_DBV_x_SBV[, 1]
      dbv.sbv.p <- truth$explFacts$U_DBV_x_SBV[, 2]

      DBVSBVSIGV$DBVxSBV_w <- cbind("true" = dbv.sbv.w,
                                    "estimated" = tmbExtract(out,
                                                             dllID = "MMTMB",
                                                             params = c("DBV_x_SBV"),
                                                             reNames = NULL,
                                                             paramsType = "random",
                                                             path = dllIDPath)[1:length(dbv.sbv.w)])
      DBVSBVSIGV$DBVxSBV_p <- cbind("true" = dbv.sbv.p,
                                    "estimated" = tmbExtract(out,
                                                             dllID = "MMTMB",
                                                             params = c("DBV_x_SBV"),
                                                             reNames = NULL,
                                                             paramsType = "random",
                                                             path = dllIDPath)[-(1:length(dbv.sbv.p))])
    }

  }
  rownames(DBVSBVSIGV$DBV) <- truth$paramsScenario$panelGenos
  return(DBVSBVSIGV)
}



#load data for Inference
TMBInputData <- function(outAll, strategy){

  #initialize
  out <- list()
  if (strategy == "sole_only") {
    out$y <- outAll$data$y
  } else if (strategy %in% c("inter_only","inter_only_ctrl")){
    out$Y <- outAll$data$Y
  } else if (strategy %in% c("sole_inter_50","sole_inter_ctrl")){
    out$Y <- outAll$data$Y
    out$y <- outAll$data$y
  }
  return(out)
}




##########     FUNCTION TO EXTRACT PARAMETERS FROM TMB OUTPUT       ############
#' @title ExtractRandTmb
#'
#' @description
#'  Function to extracts random coefficients from a list containing TMB::MakeADFun and nlminb bObjects.
#'
#'
#' @param tmbObj A list that contains the TMB::MakeADFun and nlminb objects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @return A vector of chosen parameters.
#' @author Jemay Salomon
## @examples
#'@export
ExtractRandTmb <- function(tmbObj,
                           params = NULL,
                           dllID,
                           reNames = NULL,
                           path=NULL) {

  #require tmb packages and make summary
  requireNamespace(package="TMB")

  # Load the DLL using dyn.load
  if(!is.null(path))
    dyn.load(file.path(path, TMB::dynlib(dllID)))
  else
    dyn.load(TMB::dynlib(dllID))

  sdreporttmbObj <- TMB::sdreport(tmbObj$obj) # change f to obj

  if (is.null(params)) {
    randEffs <- summary(sdreporttmbObj, select = "random")[, "Estimate"]
  } else {
    if (!is.character(params)) stop("The 'random' argument must be a character vector.")

    randEffs <- lapply(params, function(rand) {
      idx <- which(rownames(summary(sdreporttmbObj, select = "random")) == rand)
      if (length(idx) == 0) stop(paste("Random effect '", rand, "' not found in summary."))
      summary(sdreporttmbObj, select = "random")[idx, "Estimate"] })
    randEffs <- unlist(randEffs, recursive = TRUE, use.names = TRUE)
    randEffs <- as.numeric(randEffs)
    names(randEffs) <- reNames
  }

  return(unlist(randEffs))
}



#' @title ExtractParamsTmb
#'
#'@description
#' Function to extracts specified parameters from a list containing TMB::MakeADFun and nlminb Objects.
#'
#' @param tmbObj A list that contains the TMB::MakeADFun and nlminb tmbObjects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @return A vector of chosen parameters.
#' @author Jemay Salomon
## @examples
#'@export
ExtractParamsTmb <- function(tmbObj,
                             params = NULL,
                             dllID,
                             reNames = NULL,
                             path = NULL) {
  #require tmb packages and make summary
  requireNamespace(package="TMB")

  # Load the DLL using dyn.load
  if(!is.null(path))
    dyn.load(file.path(path, TMB::dynlib(dllID)))
  else
    dyn.load(TMB::dynlib(dllID))


  if (!is.list(tmbObj)) {
    stop("out must be a list")
  }

  if(is.null(params)){
    parameters <- tmbObj$fit$par
  } else {
    tmbParams <- lapply(params, function(param) {
      if (!any(grepl(paste0("^", param, "$"), names(tmbObj$fit$par)))) {
        stop(paste(param, " not found in out$fit$par"))}
      idx <- grepl(paste0("^", param, "$"), names(tmbObj$fit$par))
      return(tmbObj$fit$par[idx])})

    parameters <- (unlist(tmbParams))

    if (!is.null(reNames)) {
      stopifnot(length(reNames)==length(parameters))
      names(parameters) <- reNames
    }
  }

  return(parameters)
}



#' @title ExtractVarTmb
#'
#'@description
#' Function to extract  variances parameters from a list containing TMB::MakeADFun and nlminb Objects
#'
#' @param tmbObj A list that contains the TMB::MakeADFun and nlminb tmbObjects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @return A vector of chosen parameters.
#' @author Jemay Salomon
## @examples
#'
#'@export
ExtractVarTmb <- function(tmbObj,
                          params=NULL,
                          dllID,
                          reNames = NULL,
                          path = NULL) {

  #require tmb packages and make summary
  requireNamespace(package="TMB")

  # Load the DLL using dyn.load
  if(!is.null(path))
    dyn.load(file.path(path, TMB::dynlib(dllID)))
  else
    dyn.load(TMB::dynlib(dllID))

  if(is.null(params)){
    stop("Params must be specified")}
  idx <- ExtractParamsTmb(tmbObj, dllID = dllID,
                          params, reNames, path = path)

  var <- exp((idx))^2

  return(var)
}


#' @title  ExtractCorTmb
#'
#'@description
#' Function to extract  correlation parameters from a list containing TMB::MakeADFun and nlminb Objects
#'
#' @param tmbObj A list that contains the TMB::MakeADFun and nlminb tmbObjects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @return A vector of chosen parameters.
#' @author Jemay Salomon
## @examples
#'
#'@export
ExtractCorTmb <- function(tmbObj,
                          params = NULL,
                          dllID,
                          reNames = NULL,
                          path = NULL) {

  #require tmb packages and make summary
  requireNamespace(package="TMB")

  # Load the DLL using dyn.load
  if(!is.null(path))
    dyn.load(file.path(path, TMB::dynlib(dllID)))
  else
    dyn.load(TMB::dynlib(dllID))

  #range check
  if (!is.list(tmbObj)) {
    stop("out must be a list")
  }

  #Get report object from report(f)
  objReport = tmbObj$obj$report() # => TODO: to change in f (for simul)=> obj(real data)

  #set the output list
  out <- list()

  #set conditions parameters
  if(is.null(params)){
    Names <- list()
    for (param in 1: length(objReport)){
      Names[[param]] <- names(objReport)[[param]]
      out[[param]] <- objReport[[param]][2]
      names(out) <- Names}

  } else {
    for(param in params) {
      if (is.null(objReport[[param]])) stop(paste(param, "not found in obj$f$report()"))
      out[[param]] <- objReport[[param]][2]}

    if (!is.null(reNames)) {
      stopifnot(length(reNames) == length(out))
      names(out) <- reNames
    }
  }

  return(unlist(out))
}


#' @title  ExtractStdTmb
#'
#'@description
#' Function to extract  standard error parameters from a list containing TMB::MakeADFun and nlminb Objects
#'
#' @param tmbObj A list that contains the TMB::MakeADFun and nlminb tmbObjects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @return A matrix of chosen parameters with their standard error
#' @author Jemay Salomon
## @examples
#'
#'@export
ExtractStdTmb <- function(tmbObj,
                          params = NULL,
                          dllID,
                          reNames = NULL,
                          path = NULL){

  #require tmb packages and make summary
  requireNamespace(package="TMB")

  # Load the DLL using dyn.load
  if(!is.null(path))
    dyn.load(file.path(path, TMB::dynlib(dllID)))
  else
    dyn.load(TMB::dynlib(dllID))


  #range check
  if(!is.list(tmbObj)) stop("tmbOj must be a list")


  sdreporttmbObj <- TMB::sdreport(tmbObj$f) # change f to obj

  if (is.null(params)) {
    stdEffs <- summary(sdreporttmbObj)[, "Std. Error"]
  } else {
    if (!is.character(params)) stop("The 'params' argument must be a character vector.")

    tmp <- lapply(params, function(std) {
      idx <- which(rownames(summary(sdreporttmbObj)) == std)
      if (length(idx) == 0) stop(paste("std effect '", std, "' not found in summary."))
      summary(sdreporttmbObj)[idx, "Std. Error"] })
  }

  out <- unlist(tmp)

  if(!is.null(reNames)){
    stopifnot(length(out)==length(reNames))
    names(out) <- reNames
  }

  return(cbind("Std. Error" = out))
}


#' @title tmbExtract
#'
#'@description
#' Macro function to extract TMB parameters of specified types.
#'
#' @param tmbObj A list containing TMB::MakeADFun and nlminb tmbObjects.
#' @param params Parameter names to extract. If NULL, all parameters will be extracted.
#' @param reNames A vector of names to rename parameters. If NULL, the original TMB names will be retained.
#' @param path The path trough your dllID-object location. If it is in the current working directory, set it to NULL
#' @param paramsType Specifies the type of TMB parameters to extract (e.g., "paramsTmb", "random", "variance", "correlation").
#' @return A vector of selected parameters.
#' @author Jemay Salomon
#' @export
tmbExtract <- function(tmbObj,
                       params = NULL,
                       reNames = NULL,
                       dllID,
                       path = NULL,
                       paramsType){


  #set arguments parameters
  argsTmb <- list(tmbObj = tmbObj,
                  params = params,
                  dllID = dllID,
                  reNames = reNames,
                  path = path)

  #set conditions
  if (paramsType == "paramsTmb") {
    out <- do.call(ExtractParamsTmb, argsTmb)
  } else if (paramsType == "random") {
    out <- do.call(ExtractRandTmb, argsTmb)
  } else if (paramsType == "variance") {
    out <- do.call(ExtractVarTmb, argsTmb)
  } else if (paramsType == "correlation") {
    out <- do.call(ExtractCorTmb, argsTmb)
  } else if (paramsType == "std. error") {
    out <- do.call(ExtractStdTmb, argsTmb)
  } else {
    stop("Invalid paramsType")
  }

  #return
  return(out)

}

#compute rmse
rmse <- function(actual, predicted) {
  return(sqrt(mean((actual - predicted)^2)))
}


paramsStrat <- function(strategiesToRun, strategy, outDir, scenarioName,
                        simulId, num.cores,
                        scenario_parameters,
                        method, dllIDPath,
                        DBVxSBV=F){

  # Use mclapply to run computations in parallel
  parameters <- mclapply(strategiesToRun, function(strategy) {
    inF <- file.path(outDir, "infer", scenarioName, simulId, strategy, "out_tmb.rds")
    if (file.exists(inF)){
      out <- readRDS(inF)
      truth <- list(paramsScenario = readRDS(file.path(outDir, "simul", scenarioName, "paramsScenario.rds")),
                    explFacts = readRDS(file.path(outDir, "simul", scenarioName, simulId, "explFacts.rds")))
      # Get parameters values
      params <- get_params_value(strategy, out, names_of_parameters,
                                 truth, dllIDPath=dllIDPath, DBVxSBV)

      # Create a data frame to store the parameters for this iteration
      data_infer <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = rownames(params),
        True = params[, "true"],
        Estimated = params[, "estimated"],
        Correlation = NA,
        stringsAsFactors = FALSE
      )

      # Append data
      scenario_parameters <<- data.table::rbindlist(list(scenario_parameters, data_infer))

      # Get genetic values
      DBVSBVSIGV <- get_genetic_values(strategy, truth, out, dllIDPath, DBVxSBV)

      # Calculate correlations
      if(is.null(method))
        method <- "pearson"

      if(strategy!="sole_only")
        cor_DBV <- cor(DBVSBVSIGV$DBV[, "true"], DBVSBVSIGV$DBV[, "estimated"], method = method)
      else
        cor_DBV <- NA

      cor_SBV <- NA
      cor_SIGV <- NA
      cor_cBV <- NA
      cor_DBVxSBV_w <- NA
      cor_DBVxSBV_p <- NA

      if(strategy=="sole_only")
        cor_cBV <- cor((DBVSBVSIGV$DBV[, "true"] + DBVSBVSIGV$SIGV[, "true"]) ,
                       DBVSBVSIGV$DBV[, "estimated"], method = method)

      # Compute correlations based on strategy
      if (strategy %in% c("inter_only", "inter_only_ctrl", "sole_inter_50", "sole_inter_ctrl")) {
        cor_SBV <- cor(DBVSBVSIGV$SBV[, "true"], DBVSBVSIGV$SBV[, "estimated"], method = method)
        if(isTRUE(DBVxSBV)){
          cor_DBVxSBV_w <- cor(DBVSBVSIGV$DBVxSBV_w[, "true"], DBVSBVSIGV$DBVxSBV_w[, "estimated"], method = method)
          cor_DBVxSBV_p <- cor(DBVSBVSIGV$DBVxSBV_p[, "true"], DBVSBVSIGV$DBVxSBV_p[, "estimated"], method = method)
        }

      }

      if (strategy %in% c("sole_inter_50", "sole_inter_ctrl")) {
        cor_SIGV <- cor(DBVSBVSIGV$SIGV[, "true"], DBVSBVSIGV$SIGV[, "estimated"], method = method)

        cor_cBV <- cor((DBVSBVSIGV$DBV[, "true"] + DBVSBVSIGV$SIGV[, "true"]) ,
                       (DBVSBVSIGV$DBV[, "estimated"] + DBVSBVSIGV$SIGV[, "estimated"]), method = method)
        if(isTRUE(DBVxSBV)){
          cor_DBVxSBV_w <- cor(DBVSBVSIGV$DBVxSBV_w[, "true"], DBVSBVSIGV$DBVxSBV_w[, "estimated"], method = method)
          cor_DBVxSBV_p <- cor(DBVSBVSIGV$DBVxSBV_p[, "true"], DBVSBVSIGV$DBVxSBV_p[, "estimated"], method = method)
        }
      }

      # Create data frames for correlations
      correlation_df_DBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_DBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_DBV,
        stringsAsFactors = FALSE
      )

      correlation_df_SBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_SBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_SBV,
        stringsAsFactors = FALSE
      )

      correlation_df_SIGV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_SIGV",
        True = NA,
        Estimated = NA,
        Correlation = cor_SIGV,
        stringsAsFactors = FALSE
      )

      # Create data frames for correlations
      correlation_df_cBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_cBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_cBV,
        stringsAsFactors = FALSE
      )

      out <- list(data_infer, correlation_df_DBV, correlation_df_SBV,
                  correlation_df_SIGV, correlation_df_cBV)

      if(isTRUE(DBVxSBV)){
        correlation_df_DBVxSBV_w <- data.frame(
          Scenario = scenarioName,
          SimulID = simulId,
          Strategy = strategy,
          Parameter = "cor_DBVxSBV_w",
          True = NA,
          Estimated = NA,
          Correlation = cor_DBVxSBV_w,
          stringsAsFactors = FALSE
        )

        correlation_df_DBVxSBV_p <- data.frame(
          Scenario = scenarioName,
          SimulID = simulId,
          Strategy = strategy,
          Parameter = "cor_DBVxSBV_p",
          True = NA,
          Estimated = NA,
          Correlation = cor_DBVxSBV_p,
          stringsAsFactors = FALSE
        )

        # Return parameter data frames
        out$correlation_df_DBVxSBV_w = correlation_df_DBVxSBV_w
        out$correlation_df_DBVxSBV_p = correlation_df_DBVxSBV_p
      }

      return(out)

    }}, mc.cores = num.cores)

}

corStrat <- function(strategiesToRun, strategy, outDir, scenarioName,
                     simulId, num.cores, asreml, tmb, dllIDPath, DBVxSBV){

  results <- mclapply(strategiesToRun, function(strategy) {

    if(asreml)
      inF <- file.path(outDir, "infer", scenarioName, simulId, strategy, "out_asreml.rds")

    if(tmb)
      inF <- file.path(outDir, "infer", scenarioName, simulId, strategy, "out_tmb.rds")

    if (file.exists(inF)){
      out <- readRDS(inF)
      truth <- list(paramsScenario = readRDS(file.path(outDir, "simul", scenarioName, "paramsScenario.rds")),
                    explFacts = readRDS(file.path(outDir, "simul", scenarioName, simulId, "explFacts.rds")))

      if(asreml)
        DBVSBVSIGV <- get_genetic_values_asreml(strategy, truth, out)
      else
        DBVSBVSIGV <- get_genetic_values(strategy, truth, out,  dllIDPath, DBVxSBV)


      genos <- truth$paramsScenario$panelGenos

      True_DBV <- DBVSBVSIGV$DBV[,"true"]
      True_SBV <- DBVSBVSIGV$SBV[,"true"]
      True_SIGV <- DBVSBVSIGV$SIGV[,"true"]


      est_DBV <- DBVSBVSIGV$DBV[, "estimated"]
      est_SBV <- rep(NA, length(genos))
      est_SIGV <- rep(NA, length(genos))

      if (!strategy %in% c("sole_only")) {
        est_SBV <- DBVSBVSIGV$SBV[, "estimated"]
        # if(isTRUE(DBVxSBV)){
        #   True_DBVxSBV_w <- DBVSBVSIGV$DBVxSBV_w[,"true"]
        #   True_DBVxSBV_p <- DBVSBVSIGV$DBVxSBV_p[,"true"]
        #   est_DBVxSBV_w <- DBVSBVSIGV$DBVxSBV_w[, "estimated"]
        #   est_DBVxSBV_p <- DBVSBVSIGV$DBVxSBV_p[, "estimated"]
        # }
      }

      if (strategy %in% c("sole_inter_50", "sole_inter_ctrl")) {
        est_SIGV <- DBVSBVSIGV$SIGV[, "estimated"]
      }

      dat <- data.frame(
        Scenario = scenarioName,
        genotypes = genos,
        SimulID = simulId,
        Strategy = strategy,
        True_DBV = True_DBV,
        Estimated_DBV = est_DBV,
        True_SBV = True_SBV,
        Estimated_SBV = est_SBV,
        True_SIGV = True_SIGV,
        Estimated_SIGV = est_SIGV
      )

      # if(isTRUE(DBVxSBV)){ ## ils ne sont pas de meme tailles => erreur
      #   dat <- data.frame(
      #     Scenario = scenarioName,
      #     genotypes = genos,
      #     SimulID = simulId,
      #     Strategy = strategy,
      #     True_DBV = True_DBV,
      #     Estimated_DBV = est_DBV,
      #     True_SBV = True_SBV,
      #     Estimated_SBV = est_SBV,
      #     True_SIGV = True_SIGV,
      #     Estimated_SIGV = est_SIGV,
      #     True_DBVxSBV_w = True_DBVxSBV_w,
      #     True_DBVxSBV_p = True_DBVxSBV_p,
      #     Estimated_DBVxSBV_w = est_DBVxSBV_w,
      #     Estimated_DBVxSBV_p =est_DBVxSBV_p
      #   )
      # }

      return(dat)
    }}, mc.cores = num.cores)
}

##' FTPR Function: Calculate True Positive Rate or False Positive Rate
##'
##' This function computes either the True Positive Rate (TPR) or the False Positive Rate (FPR)
##' based on the predicted and true indices for a set of genotypes. The function allows users to specify
##' the proportion of genotypes to retain and compare the overlap between true and predicted rankings.
##'
##' @param data A data frame containing genotype data and corresponding indices.
##' @param colGenos A string specifying the column name in `data` that contains genotype identifiers.
##' @param propToKeep A numeric value (between 0 and 1) indicating the proportion of genotypes to retain based on the index.
##' @param colTrueIdx A string specifying the column name in `data` that contains the true index values for the genotypes.
##' @param colPredIdx A string specifying the column name in `data` that contains the predicted index values for the genotypes.
##' @param TPR A logical value (TRUE/FALSE). If TRUE, the function returns the True Positive Rate (TPR). If FALSE, it returns the False Positive Rate (FPR).
##'
##' @return A numeric value representing the TPR or FPR based on the specified input.
##'
##' @examples
##' # Example usage:
##' set.seed(123)
##' data <- data.frame(
##'   Genos = 1:100,
##'   TrueIdx = runif(100),
##'   PredIdx = runif(100)
##' )
##' tpr <- FTPR(data, "Genos", 0.2, "TrueIdx", "PredIdx", TPR = TRUE)
##' fpr <- FTPR(data, "Genos", 0.2, "TrueIdx", "PredIdx", TPR = FALSE)
##'
##' @seealso \code{\link{binaryClassif()}}
FTPR <- function(data, colGenos, propToKeep, colTrueIdx, colPredIdx, TPR = TRUE){

  nToKeep <- ceiling(propToKeep * nrow(data))

  genosTrueIdx <- data[order(data[[colTrueIdx]], decreasing = TRUE), ][1:nToKeep, colGenos]
  genosPredIdx <- data[order(data[[colPredIdx]], decreasing = TRUE), ][1:nToKeep, colGenos]

  genosOffTrue <- data[order(data[[colTrueIdx]], decreasing = TRUE), ][(nToKeep+1):(nrow(data)), colGenos]
  genosOffPred <- data[order(data[[colPredIdx]], decreasing = TRUE), ][(nToKeep+1):(nrow(data)), colGenos]

  TP <- length(intersect(genosTrueIdx, genosPredIdx))
  TN <- length(intersect(genosOffTrue , genosOffPred))
  FP <- length(setdiff(genosPredIdx, genosTrueIdx))
  FN <- length(setdiff(genosTrueIdx, genosPredIdx))


  if(TPR){
    pr <- TP / nToKeep
  } else {
    pr <- FP / (FP + TN)
  }
  return(pr)
}

# Selection
selectIdx <- function(tmp, strategy){

  tmp$fixIpure <- tmp$True_DBV + tmp$True_SIGV # BVsc
  tmp$fixImix <- tmp$True_DBV +  tmp$True_SBV # BVic

  if (strategy %in% c("sole_only")) {
    # tmp$trueIpure <- tmp$True_DBV + tmp$True_SIGV
    # tmp$trueImix <- tmp$True_DBV + tmp$True_SBV
    tmp$estimIpure <- tmp$Estimated_DBV # which is cBV
    tmp$estimImix <- tmp$Estimated_DBV # which is cBV
  }
  if (strategy %in% c("inter_only", "inter_only_ctrl")) {
    # tmp$trueIpure <- tmp$True_DBV
    # tmp$trueImix <- tmp$True_DBV + tmp$True_SBV
    tmp$estimIpure <- tmp$Estimated_DBV
    tmp$estimImix <- tmp$Estimated_DBV + tmp$Estimated_SBV
  }
  if (strategy %in% c("sole_inter_50", "sole_inter_ctrl")) {
    # tmp$trueIpure <- tmp$True_DBV + tmp$True_SIGV
    # tmp$trueImix <- tmp$True_DBV + tmp$True_SBV
    tmp$estimIpure <- tmp$Estimated_DBV + tmp$Estimated_SIGV
    tmp$estimImix <- tmp$Estimated_DBV + tmp$Estimated_SBV
  }
  return(tmp)
}


##' Print space
##' @param x number of interline to print (numeric and not decimal)
##' @return interline space
##' @author Jemay SALOMON
##' @examples
##' space(2)
space <- function(x){
  stopifnot(is.numeric(x))
  stopifnot(x == trunc(x))
  for (i in 1:x) {
    cat("\n")
  }
}


##' IdxMixSole Function
##'
##' This function calculates a composite index for each genotype by combining DBV, SBV, and SIS values,
##' and returns a boolean vector where `FALSE` indicates the selected genotypes based on the index,
##' and `TRUE` indicates the remaining genotypes.
##'
##' @param data A data frame containing the relevant data.
##' @param colGenos A string specifying the column name in `data` that contains genotype identifiers.
##' @param colDBV A string specifying the column name in `data` that contains the DBV values.
##' @param colSBV A string specifying the column name in `data` that contains the SBV values.
##' @param colSIGV A string specifying the column name in `data` that contains the SIS values.
##' @param propToKeep A numeric value (between 0 and 1) indicating the proportion of genotypes to retain.
##' @param w_sole A numeric weight (between 0 and 1) that indicates how much weight to give to SIS in the calculation.
##'
##' @return A list of a boolean vector where `FALSE` indicates the selected genotypes based on the index & the Idx value
##' and `TRUE` indicates the rest. The function also adds a new column `IdxMixSole` to the data frame.
##'
##' @examples
##' # Example usage:
##' set.seed(123)
##' data <- data.frame(
##'   Genos = 1:10,
##'   DBV = runif(10),
##'   SBV = c(0.5, 0.5, 0.3, 0.3, 0.7, 0.1, 0.9, 0.4, 0.4, 0.8),
##'   SIS = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
##' )
##' results <- IdxMixSole(data, "Genos", "DBV", "SBV", "SIS", 0.5, 0.7)[["Idx"]]
##' print(results)
##' results <- IdxMixSole(data, "Genos", "DBV", "SBV", "SIS", 0.5, 0.7)[["boolVec"]]
##' print(results)

IdxMixSole <- function(data, colGenos, colDBV, colSBV, colSIGV,
                       propToKeep=NULL, w_sole=NULL) {

  Genos <- data[[colGenos]]
  colDBV <- data[[colDBV]]
  colSBV <- data[[colSBV]]
  colSIGV <- data[[colSIGV]]

  Idx <- colDBV +
    ifelse(!is.na(colSBV), (1 - w_sole) * colSBV, 0) +
    ifelse(!is.na(colSIGV), w_sole * colSIGV, 0)

  data$IdxMixSole <- Idx

  if(!is.null(propToKeep) && !is.null(w_sole)){
    nToKeep <- round(propToKeep * nrow(data))
    genosSelected <- data[order(data$IdxMixSole, decreasing = TRUE), ][1:nToKeep, colGenos]
    boolVec <- rep(TRUE, length(Genos))
    boolVec[Genos %in% genosSelected] <- FALSE
  } else {
    boolVec=NULL
  }

  return(out <- list(boolVec=boolVec, Idx=Idx))
}


##' CullIdx Function
##'
##' This function identifies a subset of genotypes to retain based on specified proportions,
##' and returns a boolean vector where `FALSE` indicates the selected genotypes
##' and `TRUE` indicates the remaining genotypes.
##'
##' @param data A data frame containing the relevant data.
##' @param colGenos A string specifying the column name in `data` that contains genotype identifiers.
##' @param colIsole A string specifying the column name in `data` that contains the "Isole" values.
##' @param colImix A string specifying the column name in `data` that contains the "Imix" values.
##' @param propToKeep A numeric value (between 0 and 1) indicating the proportion of genotypes to retain.
##' @param prop.sole A numeric value (between 0 and 1) indicating the proportion of `propToKeep`
##' that should be selected based on the "Isole" column.
##'
##' @return A boolean vector where `FALSE` indicates the selected genotypes and `TRUE` indicates the non-selected.
###' @author Jemay SALOMON
##' @examples
##' # Example usage:
##' set.seed(1234)
##' data <- data.frame(Genos = 1:10, Isole = runif(10), Imix = runif(10))
##' CullIdx(data, "Genos", "Isole", "Imix", 0.4, 0.8)
##'
CullIdx <- function(data, colGenos, colIsole, colImix, propToKeep, prop.sole) {

  Genos <- data[[colGenos]]
  Isole <- data[[colIsole]]
  Imix  <- data[[colImix]]

  nToKeep <- round(propToKeep * nrow(data))
  nSoleToKeep <- round(prop.sole * nToKeep)
  nMixToKeep <- nToKeep - nSoleToKeep

  genosSole <- data[order(Isole, decreasing = TRUE), ][1:nSoleToKeep, colGenos]
  genosMix <- data[order(Imix, decreasing = TRUE), ][1:nToKeep, colGenos]
  genosMixSelected <- setdiff(genosMix, genosSole)[0:nMixToKeep]

  genosSelected <- c(genosSole, genosMixSelected)
  boolVec <- rep(TRUE, length(Genos))
  boolVec[Genos %in% genosSelected] <- FALSE
  return(boolVec)
}

##' Hypothesis testing
##'
##' Return the number of true positives, false positives, true negatives,
##' false negatives, true positive proportion (sensitivity), false positive
##' proportion, accuracy, true negative proportion (specificity), false
##' discovery proportion, false negative proportion, positive predictive
##' value (precision) and Matthews correlation coefficient.
##' More details on Wikipedia (\href{http://en.wikipedia.org/wiki/Sensitivity_and_specificity}{[1]}, \href{http://en.wikipedia.org/wiki/Matthews_correlation_coefficient}{[2]}).
##' @param known.nulls vector of booleans (TRUE if the null is true)
##' @param called.nulls vector of booleans (TRUE if the null is accepted); should be in the same order as the othr vector!
##' @return vector with names
##' @author Timothee Flutre
binaryClassif <- function(known.nulls, called.nulls){
  ## http://en.wikipedia.org/wiki/Sensitivity_and_specificity
  ##
  ##                                  CALLED
  ##                     Accepted null     Rejected null
  ##
  ##       true null         TN (U)            FP (V)          n0
  ## TRUTH
  ##       false null        FN (T)            TP (S)          n1
  ##
  ##                         a                 r               n
  stopifnot(is.vector(known.nulls), is.vector(called.nulls),
            length(known.nulls) == length(called.nulls),
            sum(! is.logical(known.nulls)) == 0,
            sum(! is.logical(called.nulls)) == 0)

  n <- length(known.nulls) # total number of tests
  n0 <- sum(known.nulls)   # nb of true nulls
  n1 <- n - n0             # nb of "false nulls" (i.e. "true alternatives")
  a <- sum(called.nulls)   # nb of accepted nulls ("called not significant")
  r <- n - a               # nb of rejected nulls ("called significant", "discoveries")

  ## true positive = reject a false null
  tp <- sum(which(! called.nulls) %in% which(! known.nulls))

  ## false positive = reject a true null (type I error, "false alarm")
  fp <- sum(which(! called.nulls) %in% which(known.nulls))

  ## true negatives = accept a true null
  tn <- sum(which(called.nulls) %in% which(known.nulls))

  ## false negatives = accept a false null (type II error, "miss")
  fn <- sum(which(called.nulls) %in% which(! known.nulls))

  tpp <- tp / n1        # true positive prop (sensitivity)
  fnp <- fn / n1        # false negative prop
  tnp <- tn / n0        # true negative prop (specificity) = 1 - fpp
  fpp <- fp / n0        # false positive prop
  fdp <- fp / r         # false discovery prop
  ppv <- tp / r         # positive predictive value (precision)
  acc <- (tp + tn) / n  # accuracy
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

  return(c(n=n, n0=n0, n1=n1,
           r=r, tp=tp, fp=fp,
           a=a, tn=tn, fn=fn,
           tpp=tpp, fnp=fnp, tnp=tnp, fpp=fpp, fdp=fdp, ppv=ppv, acc=acc,
           mcc=mcc))
}


##' Data.frame
##'
##' Convert the factor columns into character.
##' @param x data.frame
##' @return data.frame
##' @author Timothee Flutre
##' @export
convertFactorColumnsToCharacter <- function(x){
  stopifnot(is.data.frame(x))
  idx <- sapply(x, is.factor)
  x[idx] <- lapply(x[idx], as.character)
  return(x)
}

##' Data.frame
##'
##' Convert the factor columns into character.
##' @param x data.frame
##' @return data.frame
##' @author Timothee Flutre
##' @export
convertCharacterColumnsToFactor <- function(x){
  stopifnot(is.data.frame(x))
  idx <- sapply(x, is.character)
  x[idx] <- lapply(x[idx], as.factor)
  return(x)
}

###############################################################################
###############################################################################
##                                                                      ######
##                        ASREML                                        #####
##############################################################################

#bivariate
mvn_lmmASREML <- function(Y, fEffs, vEffs, listXs=NULL, listZs=NULL, listVCov=NULL,
                          dllID = NULL, contrasts = NULL,
                          verbose = TRUE, silent = FALSE, data) {

  # Compute ainverse matrix
  if (verbose)
    cat("Computing ainverse matrix...\n")

  # Initialize output
  out <- list()

  # Extract phenotype names
  phenoNames <- names(Y)
  if (length(phenoNames) < 2) stop("Y must have at least two phenotypes.")
  yNames1 <- phenoNames[1]
  yNames2 <- phenoNames[2]

  # Dynamically construct the fixed effects formula
  fixedEffectsTerms <- c("trait", paste0("at(trait):", unlist(fEffs)))
  fEffsFormula <- paste(fixedEffectsTerms, collapse = " + ")

  # Check for consistency in listZs and listVCov dimensions
  stopifnot(
    colnames(listZs[["ZD"]]) == colnames(listVCov[["K"]]),
    colnames(listZs[["ZD"]]) == rownames(listVCov[["K"]])
  )

  K <- listVCov[["K"]]

  if(verbose)
    cat("construct randoms terms ...\n")
  # Construct the random effect formula dynamically
  randomFormula <- as.formula(
    paste0("~ corgh(trait, init = c(0.1, 1, 1)):~vm(", vEffs, ", K)")
  )

  if(verbose)
    cat("fitting the model...\n")
  # Fit the model
  t0 <- proc.time()
  fit <- asreml(
    fixed = as.formula(paste("cbind(", yNames1, ",", yNames2, ") ~ ", fEffsFormula)),
    random = randomFormula,
    residual = ~ id(units):corgh(trait),
    data = data,
    na.action = na.method(x = "include", y = "include"),
    maxit = 20, trace=F
  )

  # Store fit in output
  fit <- update(fit)
  t1 <- proc.time()

  out$fit_time <- (t1-t0)[1]
  out$fit <- fit

  return(out)
}


#univariate
uv_lmmASREML <- function(Y, fEffs, vEffs, listXs = NULL,
                         listZs = NULL, listVCov = NULL,
                         dllID = NULL, verbose = TRUE, data) {

  # Compute ainverse matrix
  if (verbose)
    cat("Computing ainverse matrix...\n")


  # Initialize output
  out <- list()

  # Extract phenotype name
  phenoName <- names(Y)[1]  # Only the first phenotype is used for univariate
  if (is.null(phenoName)) stop("Y must have at least one phenotype.")

  # Dynamically construct the fixed effects formula
  fixedEffectsTerms <- c("1", unlist(fEffs))
  fEffsFormula <- paste(fixedEffectsTerms, collapse = " + ")

  # Check consistency in listZs and listVCov (if provided)
  if (!is.null(listZs) && !is.null(listVCov)) {
    if (verbose)
      cat("Checking consistency of listZs and listVCov...\n")
    stopifnot(
      colnames(listZs[["ZD"]]) == colnames(listVCov[["K"]]),
      colnames(listZs[["ZD"]]) == rownames(listVCov[["K"]])
    )
  }

  K <- listVCov[["K"]]

  vEffs <- "focal"

  # Construct the random effect formula dynamically
  if (verbose)
    cat("Constructing random effects formula...\n")
  randomFormula <- as.formula(
    paste0("~ vm(", vEffs, ", K)")
  )


  # Fit the model
  if (verbose)
    cat("Fitting the model...\n")
  t0 <- proc.time()
  fit <- asreml(
    fixed = as.formula(paste(phenoName, "~", fEffsFormula)),
    random = randomFormula,
    residual = ~ idv(units),
    data = data,
    na.action = na.method(x = "include", y = "include"),
    maxit = 20, trace=F
  )

  # Store fit in output
  out$fit <- update(fit)
  t1 <- proc.time()
  out$fit_time <- (t1-t0)[1]
  return(out)
}


#multivariate
uv_mvn_lmmASREML <- function(Y, fEffs, vEffs, listXs = NULL,
                             listZs = NULL, listVCov = NULL,
                             dllID = NULL,verbose = TRUE, data) {

  out <- list()
  return(out)

}


loadASRemlData  <- function(outAll, strategy){

  vEffs <- "focal"

  if (strategy == "sole_only") {

    Y <- list(yield="yield")
    data <- outAll$data$datPur
    fEffs <- list(block="block")

  } else if (strategy %in% c("inter_only","inter_only_ctrl")){
    Y <- list(yield.wheat="yield.wheat", yield.pea="yield.pea")
    data <- outAll$data$datmix
    fEffs <- list(block="block", neighbors="neighbors")

  } else if (strategy %in% c("sole_inter_50","sole_inter_ctrl")){

    #not implement
  }

  return(list(vEffs = vEffs, fEffs=fEffs, Y=Y, data=data))

}


paramsStrat_asreml <- function(strategiesToRun, strategy, outDir, scenarioName,
                               simulId, num.cores, scenario_parameters, method){

  # Use mclapply to run computations in parallel
  parameters <- mclapply(strategiesToRun, function(strategy) {
    inF <- file.path(outDir, "infer", scenarioName, simulId, strategy, "out_asreml.rds")
    if (file.exists(inF)){
      out <- readRDS(inF)
      truth <- list(paramsScenario = readRDS(file.path(outDir, "simul", scenarioName, "paramsScenario.rds")),
                    explFacts = readRDS(file.path(outDir, "simul", scenarioName, simulId, "explFacts.rds")))

      # Get parameters values
      params <- get_params_value_asreml(strategy, out, truth)

      # Create a data frame to store the parameters for this iteration
      data_infer <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = rownames(params),
        True = params[, "true"],
        Estimated = params[, "estimated"],
        Correlation = NA,
        stringsAsFactors = FALSE
      )

      # Append data
      scenario_parameters <- data.table::rbindlist(list(scenario_parameters, data_infer))

      # Get genetic values
      DBVSBVSIGV <- get_genetic_values_asreml(strategy, truth, out)

      # Calculate correlations
      if(is.null(method))
        method <- "pearson"

      if(strategy!="sole_only")
        cor_DBV <- cor(DBVSBVSIGV$DBV[, "true"], DBVSBVSIGV$DBV[, "estimated"], method = method)
      else
        cor_DBV <- NA

      cor_SBV <- NA
      cor_SIGV <- NA
      cor_cBV <- NA

      if(strategy=="sole_only")
        cor_cBV <- cor((DBVSBVSIGV$DBV[, "true"] + DBVSBVSIGV$SIGV[, "true"]) ,
                       DBVSBVSIGV$DBV[, "estimated"], method = method)

      # Compute correlations based on strategy
      if (strategy %in% c("inter_only", "inter_only_ctrl", "sole_inter_50", "sole_inter_ctrl")) {
        cor_SBV <- cor(DBVSBVSIGV$SBV[, "true"], DBVSBVSIGV$SBV[, "estimated"], method = method)
      }

      if (strategy %in% c("sole_inter_50", "sole_inter_ctrl")) {
        cor_SIGV <- cor(DBVSBVSIGV$SIGV[, "true"], DBVSBVSIGV$SIGV[, "estimated"], method = method)

        cor_cBV <- cor((DBVSBVSIGV$DBV[, "true"] + DBVSBVSIGV$SIGV[, "true"]) ,
                       (DBVSBVSIGV$DBV[, "estimated"] + DBVSBVSIGV$SIGV[, "estimated"]), method = method)
      }

      # Create data frames for correlations
      correlation_df_DBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_DBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_DBV,
        stringsAsFactors = FALSE
      )

      correlation_df_SBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_SBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_SBV,
        stringsAsFactors = FALSE
      )

      correlation_df_SIGV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_SIGV",
        True = NA,
        Estimated = NA,
        Correlation = cor_SIGV,
        stringsAsFactors = FALSE
      )

      # Create data frames for correlations
      correlation_df_cBV <- data.frame(
        Scenario = scenarioName,
        SimulID = simulId,
        Strategy = strategy,
        Parameter = "cor_cBV",
        True = NA,
        Estimated = NA,
        Correlation = cor_cBV,
        stringsAsFactors = FALSE
      )

      # Return parameter data frames
      return(list(data_infer, correlation_df_DBV, correlation_df_SBV,
                  correlation_df_SIGV, correlation_df_cBV))
    }}, mc.cores = num.cores)

}


##get_true_parameters from ScenarioParams
get_true_params_asreml <- function(strategy, truth) {

  if (strategy == "sole_only") {
    betaTr <- truth$paramsScenario$block1SpeciesEffs[,"pure", "wheat"]
    varTr <- c(truth$paramsScenario$sigma2.cBVw, truth$paramsScenario$sigma2_err_w)
    true <- c(betaTr, varTr)


  } else if (strategy == "inter_only"||strategy == "inter_only_ctrl") {
    betaTr <- c(truth$paramsScenario$block1SpeciesEffs[,"mixed",],
                truth$paramsScenario$DBVSBVEffs$neighbors1,
                truth$paramsScenario$DBVSBVEffs$focal1)

    varTr <- c(truth$paramsScenario$sigma2.DBVw, truth$paramsScenario$sigma2.SBVw,
               truth$paramsScenario$sigma2_err_w_mix, truth$paramsScenario$sigma2_err_p)
    corTr <- c(truth$paramsScenario$cor_DBVw_SBVw, truth$paramsScenario$cor_E)

    true <- c(betaTr, varTr, corTr)


  } else {
    betaTr <- c(truth$paramsScenario$block1SpeciesEffs[,"pure","wheat"],
                truth$paramsScenario$block1SpeciesEffs[,"mixed", ],
                truth$paramsScenario$DBVSBVEffs$neighbors1,
                truth$paramsScenario$DBVSBVEffs$focal1)

    varTr <- c(truth$paramsScenario$sigma2.DBVw, truth$paramsScenario$sigma2.SBVw,
               truth$paramsScenario$sigma2.SIGVw, truth$paramsScenario$sigma2_err_w_mix,
               truth$paramsScenario$sigma2_err_p, truth$paramsScenario$sigma2_err_w)

    corTr <- c(truth$paramsScenario$cor_DBVw_SBVw, truth$paramsScenario$cor_E)

    true <- c(betaTr, varTr, corTr)
  }


  return(true)
}


#get estimated_values from asreml
get_estimated_params_asreml <- function(strategy, out) {

  coef.fixed <- summary.asreml(out$fit, coef=T)$coef.fixed
  varcomp <-    summary.asreml(out$fit)$varcomp

  if(strategy =="sole_only"){
    mu_w_pur <- coef.fixed[1,1]
    block_A <- coef.fixed["block_A",1]
    block_B <-coef.fixed["block_B",1]
    sigma2.cBV <-varcomp["vm(focal, K)",1]
    sigma2_err_w <- varcomp["units!units",1]
    block1_w_pur <- (block_A-block_B)/2
    betaEs <- c("mu_w_pur"=mu_w_pur, "block1_w_pur"=block1_w_pur)
    varEs <- c("sigma2.cDBV"=sigma2.cBV, "sigma2_err_w"=sigma2_err_w)
    estimated <- c(betaEs, varEs)
  }

  if (strategy %in% c("inter_only", "inter_only_ctrl")) {
    mu_w_mix <-coef.fixed["trait_yield.wheat",1]
    mu_p_mix <-coef.fixed["trait_yield.pea",1]
    block_A_w <- coef.fixed["at(trait, 'yield.wheat'):block_A",1]
    block_B_w <-coef.fixed["at(trait, 'yield.wheat'):block_B",1]
    block_A_p <- coef.fixed["at(trait, 'yield.pea'):block_A",1]
    block_B_p <-coef.fixed["at(trait, 'yield.pea'):block_B",1]

    neighbors_pea01 <-coef.fixed["at(trait, 'yield.wheat'):neighbors_pea01",1]
    neighbors_pea02 <-coef.fixed["at(trait, 'yield.wheat'):neighbors_pea02",1]
    focal_pea01 <-coef.fixed["at(trait, 'yield.pea'):neighbors_pea01", 1]
    focal_pea02 <-coef.fixed["at(trait, 'yield.pea'):neighbors_pea02", 1]

    sigma2.DBV <- varcomp["trait:vm(focal, K)!trait_yield.wheat", "component"]
    sigma2.SBV <-varcomp["trait:vm(focal, K)!trait_yield.pea", "component"]
    sigma2_err_w_mix <- varcomp["units:trait!trait_yield.wheat", "component"]
    sigma2_err_p_mix <- varcomp["units:trait!trait_yield.pea", "component"]
    Cor_BV <- varcomp["trait:vm(focal, K)!trait!yield.pea:!trait!yield.wheat.cor", "component"]
    Cor_E <-varcomp["units:trait!trait!yield.pea:!trait!yield.wheat.cor", "component"]
    focal1 = (focal_pea01 - focal_pea02)/2
    neighbors1 = (neighbors_pea01 -neighbors_pea02)/2
    block1_w_mix <- (block_A_w-block_B_w)/2
    block1_p_mix <- (block_A_p-block_B_p)/2

    betaEs <-c("mu_w_mix"=mu_w_mix, "block1_w_mix"=block1_w_mix,"neighbors1"=neighbors1,
               "mu_p_mix"=mu_p_mix, "block1_p_mix"=block1_p_mix, "focal1"=focal1)

    varEs <- c("sigma2.DBV"=sigma2.DBV, "sigma2.SBV"=sigma2.SBV,
               "sigma2_err_w_mix"=sigma2_err_w_mix, "sigma2_err_p_mix"=sigma2_err_p_mix)

    corEs <- c("Cor_BV"=Cor_BV,"Cor_E"=Cor_E)

    estimated <- c(betaEs, varEs, corEs)
  }

  if (strategy %in% c("sole_inter_50", "sole_inter_ctrl")){
    estimated <- NA
  }


  return(estimated)
}


#get parameters value
get_params_value_asreml <- function(strategy, out, truth){

  #get estimated and true of parameters values
  estimated <- get_estimated_params_asreml(strategy, out)
  true <- get_true_params(strategy, truth)
  params <- cbind("estimated" = round(estimated, 2),
                  "true" = round(true, 2))
  #rownames parameters
  rownames(params) <- namesParams(strategy, DBVxSBV)

  return(params)
}

######################################################################################
########################              simulate               ######################
######################################################################################

##' Simulate from a TMB fitted model
##'
##' @method simulate MMTMB
##' @param object fitted TMB (must have $obj and $fit$parfull)
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param ... extra arguments (not used)
##' @return list of length nsim, each entry is a list with $y and $Y, and random effects
##' @importFrom stats simulate
##' @author Jemay SALOMON
##' @export
##'
simulate.MMTMB <- function(object, nsim = 1, seed = NULL, ...) {

  # Handle RNG state as in stats::simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  # Each simulation returns a list — keep only names that exist
  ret <- replicate(
    nsim,
    {
      sims <- object$obj$simulate(par = object$fit$parfull)
      out  <- list()

      if ("yobs" %in% names(sims))           out$y          <- sims$yobs
      if ("yobs_t" %in% names(sims))         out$y_t        <- sims$yobs_t
      if ("BV_sim" %in% names(sims))         out$BV         <- sims$BV_sim
      if ("Yobs" %in% names(sims))           out$Y          <- sims$Yobs
      if ("DBV_x_SBV_sim" %in% names(sims))  out$DBV_x_SBV  <- sims$DBV_x_SBV_sim
      if ("SIGV_sim" %in% names(sims))       out$SIGV       <- sims$SIGV_sim
      if ("u_sim" %in% names(sims))          out$u          <- sims$u_sim

      out
    },
    simplify = FALSE
  )
  names(ret) <- sprintf("simul%03d", seq_len(nsim))
  attr(ret, "seed") <- RNGstate
  ret
}

