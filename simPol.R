#!/usr/bin/env Rscript

# Remove scientific notation in the entire R session
options(scipen = 100)

# load libraries
library(optparse)

#### parse command line arguments ####
parser <-
  OptionParser(prog = "./simPol.R",
               description = "Simulate RNA polymerase dynamics with transcription parameters.")

parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print messages [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print no messages")
# parser <- add_option(parser, c("-p", "--threads"), action="store", type="integer",
#                      default=1, help="Number of threads are used [default %default], currently only supports single thread",
#                      metavar="integer")

parser <- add_option(parser, c("-k", "--tssLen"), action="store", type="integer",
                     default=50, help="define the mean of pause sites across cells [default %defaultbp]",
                     metavar="integer")
parser <- add_option(parser, c("--kSd"), action="store", type="double",
                     default=0, help="define the standard deviation of pause sites across cells [default %default]",
                     metavar="double")
parser <- add_option(parser, c("--kMin"), action="store", type="integer",
                     default=17, help="upper bound of pause site allowed [default %defaultbp]",
                     metavar="integer")
parser <- add_option(parser, c("--kMax"), action="store", type="integer",
                     default=200, help="lower bound of pause site allowed [default %defaultbp]",
                     metavar="integer")

# parser <- add_option(parser, c("-m", "--ttsLen"), action="store", type="integer",
#                      default=250, help="define the length of termination region [default %defaultbp]",
#                      metavar="integer")
parser <- add_option(parser, c("--geneLen"), action="store", type="integer",
                     default=2e3, help="define the length of the whole gene [default %defaultbp]",
                     metavar="integer")

parser <- add_option(parser, c("-a", "--alpha"), action="store", type="double",
                     default=1, help="initiation rate [default %default event per min]",
                     metavar="double")
parser <- add_option(parser, c("-b", "--beta"), action="store", type="double",
                     default=1, help="pause release rate [default %default event per min]",
                     metavar="double")
# parser <- add_option(parser, c("-g", "--gamma"), action="store", type="double",
#                      default=1, help="termination rate [default %default event per min]",
#                      metavar="double")

parser <- add_option(parser, c("-z", "--zeta"), action="store", type="double",
                     default=2000, help="the mean of elongation rates across sites [default %defaultbp per min]",
                     metavar="double")
parser <- add_option(parser, c("--zetaSd"), action="store", type="double",
                     default=1000,
                     help="the standard deviation of elongation rates across sites [default %default]",
                     metavar="double")
parser <- add_option(parser, c("--zetaMax"), action="store", type="double",
                     default=2500,
                     help="the maximum of elongation rates allowed [default %defaultbp per min]",
                     metavar="double")
parser <- add_option(parser, c("--zetaMin"), action="store", type="double",
                     default=1500,
                     help="the minimum of elongation rates allowed [default %defaultbp per min]",
                     metavar="double")
parser <- add_option(parser, c("--zetaVec"), action="store", type="character",
                     default=NULL,
                     help="a file contains vector to scale elongation rates. All cells share the same set of parameters [default %default]",
                     metavar="character")

parser <- add_option(parser, c("-n", "--cellNum"), action="store", type="integer",
                     default=10, help="Number of cells being simulated [default %default]",
                     metavar="integer")
parser <- add_option(parser, c("-s", "--polSize"), action="store", type="integer",
                     default=33, help="Polymerase II size [default %defaultbp]",
                     metavar="integer")
parser <- add_option(parser, c("--addSpace"), action="store", type="integer",
                     default=17, help="Additional space in addition to RNAP size [default %defaultbp]",
                     metavar="integer")
parser <- add_option(parser, c("-t", "--time"), action="store", type="double",
                     default=0.1, help="Total time of simulating data in a cell [default %default min]",
                     metavar="double")

parser <- add_option(parser, c("--divergent"), action="store", type="logical",
                     default=TRUE, help="Generate a bigwig file on the opposite strand, reflecting divergent transcription [default %default]",
                     metavar="logical")
parser <- add_option(parser, c("--auxiliary"), action="store", type="logical",
                     default=TRUE, help="Generate auxiliary files for sampling or troubleshoot [default %default]",
                     metavar="logical")

parser <- add_option(parser, c("--continue"), action="store", type="logical",
                     default=FALSE, help="Continue simulations based on existing files from previous round [default %default]",
                     metavar="logical")
parser <- add_option(parser, c("--prefix"), action="store", type="character",
                     default=NULL, help="Prefix for file name [default %default]",
                     metavar="character")
parser <- add_option(parser, c("-d", "--outputDir"), action="store", type="character",
                     default=".", help="Directory for saving results [default %default]",
                     metavar="character")

cmd_args <- parse_args(parser)
# print_help(parser)

#### load additional libraries ####
suppressPackageStartupMessages({
  library(Matrix)
  library(rtracklayer)})

#### assign arguments ####
prefix <- cmd_args$prefix
result_dir <- cmd_args$outputDir

dir.create(result_dir, showWarnings=F, recursive=T)

# setting cores
# threads <- cmd_args$threads

# key parameters
k <- cmd_args$tssLen
ksd <- cmd_args$kSd
kmin <- cmd_args$kMin
kmax <- cmd_args$kMax

# m <- cmd_args$ttsLen
# l <- cmd_args$geneLen - k -m
l <- cmd_args$geneLen - k

alpha <- cmd_args$alpha # initiation rate
beta <- cmd_args$beta # pause release rate
# gamma <- cmd_args$gamma # termination rate

zeta <- cmd_args$zeta # elongation rate
zeta_sd <- cmd_args$zetaSd
zeta_max <- cmd_args$zetaMax
zeta_min <- cmd_args$zetaMin

zeta_vec <- cmd_args$zetaVec

cell_num <- cmd_args$cellNum
s <- cmd_args$polSize # Pol II size
h <- cmd_args$addSpace # additional space due to unknown steric hindrance
steric_hindrance <- s + h

delta_t <- 1 * 1e-4 # should be a small number
t <- cmd_args$time
steps <- t / delta_t

# total_sites <- k + l + m + 1
total_sites <- cmd_args$geneLen + 1 # total gene length + a free state

continue <- cmd_args$continue

steps_to_record <- 100 # number of final steps to take records

# make string to record params for figure and animation names
name_params <- paste0("k", k, "ksd", ksd, "kmin", kmin, "kmax", kmax,  "l", l,
                      # "m", m, 
                      "a", alpha, "b", beta,
                      # "g", gamma,
                      "z", zeta, "zsd", zeta_sd, "zmin", zeta_min, "zmax", zeta_max, 
                      "t", t, "n", cell_num, "s", s, "h", h)
# name_params <- stringr::str_remove_all(name_params, "\\.")

if (cmd_args$verbose) {
  message("The simulation is now running with parameters: ", name_params, ".")
  # message("It can take a couple of hours if thousands of cells are simulated.")
  # message('Check more info about the default parameters with "./simPol.R --help".')
}

# rename files if prefix is provided
if (!is.null(prefix)) name_params <- prefix

#### simulate one gene in a cell ####
# read in RNAP position and probability vector if continue simulation based on previous round
if (continue) {
  m_pos <- readRDS(file.path(result_dir, paste0(name_params, ".RDS")))
  pos_init <- m_pos[[length(m_pos)]]$pos
  prob_init <- readRDS(file.path(result_dir, paste0(name_params, "_prob_init.RDS")))
}

if (!exists("pos_init")) {
  # initialize a vector to hold Pol II presence and absence
  pos_init <- Matrix(data = 0, nrow = total_sites, ncol = cell_num)
  pos_init[1, ] <- 1
}

if (!exists("prob_init")) {
  # construct a probability matrix to control RNAP movement
  # Generate pause sites located from kmin to kmax with sd = ksd
  # https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
  y <- vector()
  while (length(y) < cell_num) {
    x <- rnorm(cell_num * 2, mean = k + 1, sd = ksd)
    x <- x[(x >= kmin + 1) & (x <= kmax + 1)]
    y <- c(y, x)
  }
  y <- round(y[1:cell_num]) 
  
  # a matrix of probabilities to control transition from state to state
  # cols are cells, rows are positions
  # ignore gamma for now
  
  if (is.null(zeta_vec)) {
    
    zv <- vector()
    while (length(zv) < total_sites * cell_num) {
      x <- rnorm(total_sites * cell_num * 2, mean = zeta, sd = zeta_sd)
      x <- x[(x >= zeta_min) & (x <= zeta_max)]
      zv <- c(zv, x)
    }
    zv <- zv[1:(total_sites * cell_num)] 
    
    prob_init <- Matrix(data = zv * delta_t, nrow = total_sites, ncol = cell_num)
    
  } else {
    message("Scaling vector exists, use it to scale elongation rates.")
    zv <- as.double(readLines(zeta_vec))
    if (length(zv) >= total_sites) {
      zv <- zv[1:length(total_sites)]
    } else if (length(zv) == total_sites - 1) {
      zv <- c(mean(zv), zv)
    } else {
      stop("Vector for scaling zeta is too short, check total length of the vector!")
    }
    prob_init <- Matrix(data = zeta * zv * delta_t, nrow = total_sites, ncol = cell_num)
  }
  
  prob_init[1, ] <- alpha * delta_t
  prob_init[cbind(y, 1:cell_num)] <- beta * delta_t
  
}

# total length for all sites across all cells
total_length <- total_sites * cell_num
# initialize RNAP position
pos <- pos_init
# initialize a list to record RNAP positions
m_pos <- list()
# index for the list to record positions
j <- 1

# start simulation
used_time <-
  system.time(
    for (i in 1:steps) {
      # calculate the transition probability
      prob <- prob_init * pos
      
      #### determine whether polymerase can move or not ####
      # criteria 1, probability larger than random draw
      # find which positions has probability larger than 0
      pos_pending <- which(prob > 0)
      # randomly draw numbers for comparison
      draws <- runif(n = length(pos_pending), min = 0, max = 1)
      c1 <- prob[pos_pending] > draws
      
      # criteria 2, enough space ahead to let polymerase advance
      if (length(c1) > 1) {
        # check if space ahead is larger than polymerase size
        c2 <- diff(pos_pending, lag = 1, differences = 1) > steric_hindrance
        # and always allow the polymerase at the end to move
        c2 <- c(c2, TRUE) } else {
        # allow it to move if there is only one RNAP  
          c2 <- c(TRUE)
        }
      
      # always allow the last polymerase in every cell to move
      c2[cumsum(colSums(pos > 0))] <- TRUE
      
      # decide which polymerase can eventually move
      pos_origin <- pos_pending[c1 & c2]
      # set original positions to 0
      pos[pos_origin] <- 0
      # set new positions where polII are moving to 1
      pos_moving <- pos_origin + 1
      pos[pos_moving[pos_moving <= total_length]] <- 1
      # ensure there are always polymerases waiting to be initialized (i.e., first row equals 1)
      pos[1, ] <- 1
      # record additional info for studying steric hindrance
      if (i > (steps - steps_to_record)) {
        m_pos[[j]] <-
          list("pos" = pos, "c1" = c1, "c2" = c2, "pos_pending" = pos_pending)
        j <- j + 1
      }
    }
  )

if (cmd_args$verbose) message("Total used time for ", name_params, " is ",
                              round(used_time[3] / 60, digits = 2), " mins.")

# combine data from all cells and save a BigWig file
res_all <- rowSums(pos)
# get rid of position 1, which is always equal to cell number
res_all <- res_all[-1]

# select a starting point for our simulated gene
start_point <- 0.99 * 1e6

# generate bigwig for positive strand
bw_grng_pos <-
  GRanges(seqnames = "chr1",
          IRanges(start = (start_point + 1):(start_point + total_sites - 1), width = 1),
          score = res_all,
          strand = "+",
          seqlengths = c("chr1" = cmd_args$geneLen * 10)  + start_point)

if (cmd_args$divergent) {
  # generate bigwig for negative strand
  # sample 75% number of cells to use as signals on negative strand
  cells_neg <- sample(1:cell_num, size = round(cell_num * 0.75), replace = TRUE)
  pos_neg <- pos[, cells_neg, drop = FALSE]
  # combine data from all cells and save a BigWig file
  res_all_neg <- rowSums(pos_neg)
  # get rid of position 1, which is always equal to cell number
  res_all_neg <- res_all_neg[-1]
  # reverse the order of values for putting them on negative strand
  res_all_neg <- rev(res_all_neg)
  
  bw_grng_neg <-
    GRanges(seqnames = "chr1",
            IRanges(start = seq((start_point - (total_sites + 18)), by = 1, length.out = total_sites -1),
                    width = 1),
            score = - res_all_neg,
            strand = "-",
            seqlengths = c("chr1" = cmd_args$geneLen * 10)  + start_point)
  # export bigwig file
  export(bw_grng_neg, file.path(result_dir, paste0(name_params, ".minus.bw")))
}

# export bigwig file
export(bw_grng_pos, file.path(result_dir, paste0(name_params, ".plus.bw")))

if (cmd_args$auxiliary) {
  # save res for later examination
  saveRDS(object = m_pos, file = file.path(result_dir, paste0(name_params, "_pos.RDS")))
  saveRDS(object = prob_init, file = file.path(result_dir, paste0(name_params, "_prob_init.RDS")))
}
