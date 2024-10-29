#!/usr/bin/env nextflow
// Copyright (C) 2021 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Params
params.midas_dir = ''
// params.map = '' // TO DO: We need one map per species in the end
params.map_dir = ''
params.outdir = 'output'

// Process MIDAS params
params.depth_thres = 5
params.npop_thres = 3
params.max_sites = 0

// bern params
params.q_thres = 0.1
params.iter = 3000
params.warmup = 2000
params.chains = 4
params.vp = 5
params.vq = 5

// Process params
midas_dir = file(params.midas_dir)
map_dir = file(params.map_dir)
// mapfile = file(params.map)

MIDAS = Channel.fromPath("$midas_dir/*", type: 'dir', maxDepth: 1)
  .map{ midas_dir -> tuple(midas_dir.name,
    file(midas_dir)) }

MAPS = Channel.fromPath("$map_dir/*", type: 'file', maxDepth: 1)
  .map{ mapfile -> tuple(mapfile.name.replaceAll(/\.tsv/, ''),
    file(mapfile))}

process midas2bern {
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/sites", mode: 'rellink',
    pattern: "output/sites.tsv", saveAs: {"${spec}.tsv"}
  publishDir "$params.outdir/pops", mode: 'rellink',
    pattern: "output/pops.tsv", saveAs: {"${spec}.tsv"}

  input:
  // tuple spec, file(midas_dir) from MIDAS
  // file mapfile from mapfile
  tuple spec, file(midas_dir), file(mapfile) from MIDAS.join(MAPS)
  val depth_thres from params.depth_thres
  val npop_thres from params.npop_thres
  val max_sites from params.max_sites

  output:
  tuple spec, file("output/sites.tsv") optional true into SITES
  tuple spec, file("output/pops.tsv") optional true into POPS

  """
  Rscript ${workflow.projectDir}/midas2bern.r \
    $midas_dir \
    $mapfile \
    --depth_thres $depth_thres \
    --n_thres $npop_thres \
    --outdir output \
    --max_sites $max_sites
  """
}


process bern_mix{
  cpus params.chains
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/stan_models", mode: 'rellink',
    pattern: "output/m1.stan.rdat", saveAs: {"${spec}.stan.rdat"}
  publishDir "$params.outdir/p_directional", mode: 'rellink',
    pattern: "output/p_directional.tsv.gz", saveAs: {"${spec}.tsv.gz"}
  publishDir "$params.outdir/model_summaries", mode: 'rellink',
    pattern: "output/model_summaries.tsv.gz", saveAs: {"${spec}.tsv.gz"}
  publishDir "$params.outdir/CHECK_RHAT", mode: 'rellink',
    pattern: "CHECK_RHAT", saveAs: {"$spec"}

  input:
  tuple spec, file(sites_file) from SITES
  val q_thres from params.q_thres
  val min_patients from params.npop_thres
  val iter from params.iter
  val warmup from params.warmup
  val chains from params.chains
  val vp from params.vp
  val vq from params.vq

  output:
  file "output/m1.stan.rdat" optional true
  file "output/p_directional.tsv.gz" optional true
  file "output/model_summaries.tsv.gz" optional true
  file "CHECK_RHAT" optional true

  """
  Rscript $workflow.projectDir/bern_mix.r \
    $sites_file \
    --q_thres $q_thres \
    --min_patients $min_patients \
    --outdir output \
    --iter $iter \
    --warmup $warmup \
    --chains $chains \
    --vp $vp \
    --vq $vq
  """

}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '200h'
  memory = '10G'
  withLabel: 'r'{
    module = 'R/3.6.1'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
