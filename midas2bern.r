#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# This file is part of This program.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(argparser)

#' MIDAS merged SNPs to per site distributions
#'
#' @param midas_dir 
#' @param meta 
#' @param outdir 
#' @param write_tables 
#' @param depth_thres 
#' @param npop_thres 
#' @param max_sites 
#'
#' @return
#' @export
#' @author Sur from Fraser Lab.
midas2sitesdist <- function(midas_dir, meta, outdir = "./",
                            write_tables = FALSE,
                            depth_thres = 5,
                            npop_thres = 3,
                            max_sites = 0){

  # midas_dir <- args$midas_dir
  # meta <- meta
  # outdir <- args$outdir
  # write_tables <- TRUE
  # depth_thres <- args$depth_thres
  # npop_thres <- args$n_thres
  # max_sites <- args$max_sites
  
  
  if(missing(midas_dir) || missing(meta)){
    stop("ERROR: You must provide a midas_dir & meta(data) tibble.",
         call. = TRUE)
  }
  if(!is_tibble(meta)){
    stop("ERROR: meta must be a tibble", call. = TRUE)
  }
  if(!all(c("pt", "start", "end") %in% colnames(meta))){
    stop("ERROR: columns pt, start &* end must be present in meta",
         call. = TRUE)
  }else{
    meta <- meta %>%
      select(pt, start, end)
  }
  
  ### Read and process midas data
  midas_map <- meta %>%
    pivot_longer(-pt,values_to = "sample", names_to = "timepoint") %>%
    rename(Group = pt)
  cat("Reading data...\n")
  Dat <- HMVAR::read_midas_data(midas_dir = midas_dir,
                                map = midas_map)
  
  cat("Keeping only samples from patients with start & end")
  sample_ids <- setdiff(colnames(Dat$freq), "site_id")
  midas_map <- midas_map %>%
    dplyr::filter(sample %in% sample_ids) %>%
    dplyr::group_by(Group) %>%
    dplyr::filter(dplyr::n() == 2)
  
  
  if(nrow(midas_map) < 1){
    warning("Samples in midas_dir do not include any complete pop...\n")
    print(midas_map, n = nrow(midas_map))
    return(list(Sites = NULL, Pops = NULL))
  }
  
  Dat$freq <- Dat$freq %>%
    dplyr::select(site_id, midas_map$sample)
  Dat$depth <- Dat$depth %>%
    dplyr::select(site_id, midas_map$sample)
  
  if(max_sites){
    cat("!!Keeping only max_sites...\n")
    Dat$info <- Dat$info[ 1:max_sites, ]
    Dat$freq <- Dat$freq %>%
      dplyr::filter(site_id %in% Dat$info$site_id)
    Dat$depth <- Dat$depth %>%
      dplyr::filter(site_id %in% Dat$info$site_id)
  }
  
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq,
                                     depth = Dat$depth,
                                     info = Dat$info %>%
                                       dplyr::select(site_id, ref_id, ref_pos),
                                     map = midas_map %>%
                                       dplyr::rename(pt = Group), 
                                     depth_thres = depth_thres)
  
  if(nrow(Dat) < 1){
    warning("No sites passing depth threshold...\n")
    print(midas_map, n = nrow(midas_map))
    print(Dat)
    return(list(Sites = NULL, Pops = NULL))
  }
  
  cat("Calculating change at every position in every population...\n")
  Dat <- Dat %>%
    split(.$pt) %>%
    map_dfr(function(d){
      d %>%
        pivot_wider(id_cols = c("site_id", "ref_id", "ref_pos", "pt"),
                    names_from = "timepoint",
                    values_from = c("sample", "freq", "depth")) %>%
        filter(!(is.na(freq_start) | is.na(freq_end)))
    }) %>%
    group_by(site_id) %>%
    filter(n() >= npop_thres) %>%
    ungroup() %>%
    mutate(maf_change = freq_end - freq_start)
  
  if(nrow(Dat) < 1){
    warning("No sites in enough populations...\n")
    print(midas_map, n = nrow(midas_map))
    print(Dat)
    return(list(Sites = NULL, Pops = NULL))
  }
  
  cat("Calculating per-site distribution\n")
  Sites <- Dat %>%
    group_by(site_id) %>%
    summarise(n_decrease = sum(maf_change < 0),
              n_equal = sum(maf_change == 0),
              n_increase = sum(maf_change > 0),
              n_patients = length(pt),
              .groups = "drop")
  
  cat("Calculating per-population distribution")
  Pts <- Dat %>%
    group_by(pt) %>%
    summarise(n_decrease = sum(maf_change < 0),
              n_equal = sum(maf_change == 0),
              n_increase = sum(maf_change > 0),
              n_sites = length(site_id),
              .groups = 'drop') %>%
    # filter(n_sites >= args$prop_thres * max(n_sites)) %>%
    mutate(p_increase = n_increase / n_sites,
           p_decrease = n_decrease / n_sites,
           p_change = (n_increase + n_decrease) / n_sites)
  
  
  if(write_tables){
    # Prepare output dir
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    cat("Writing tables..\n")
    write_tsv(Sites, file.path(outdir, "sites.tsv"))
    write_tsv(Pts, file.path(outdir, "pops.tsv"))
  }
  
  return(list(Sites = Sites, Pops = Pts))
}

process_arguments <- function(){
  p <- arg_parser(paste("Prepare data from MIDAS merged for bern."))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("MIDAS merged SNPs directory"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("Mapping file"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--depth_thres",
                     help = paste(""),
                     type = "numeric",
                     default = 5)
  p <- add_argument(p, "--n_thres",
                    help = paste(""),
                    type = "numeric",
                    default = 3)
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--max_sites",
                    help = paste("Zero is all sites"),
                    type = "numeric",
                    default = 0)
                  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
# args <- list(midas_dir = "MGYG-HGUT-02478/",
#              map = "metadata_day0_day14_fos.tsv",
#              depth_thres = 5,
#              n_thres = 3,
#              outdir = "output",
#              max_sites = 1000)

# This is taking from orginal pipeline and cleaning & organizing. It is only going
# to take two timepoints. Should simplify future analysis.
# This can probably go in HMVAR, but bern model needs to be its own thing
library(tidyverse)

cat("Reading mapping file...\n")
meta <- read_tsv(args$map,
                 col_type = cols(pt = col_character(),
                                 start = col_character(),
                                 end = col_character()))


Res <- midas2sitesdist(midas_dir = args$midas_dir,
                       meta = meta,
                       outdir = args$outdir,
                       write_tables = TRUE,
                       depth_thres = args$depth_thres,
                       npop_thres = args$n_thres,
                       max_sites = args$max_sites)
# Res
cat(nrow(Res$Sites), "sites in", nrow(Res$Pops), "populations remained\n")


