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
params.input = ''
params.outdir = 'output'
params.q_thres = 0.1
params.min_patients = 5
params.iter = 3000
params.warmup = 2000
params.chains = 4
params.vp = 5
params.vq = 5

// Process params
indir = file(params.input)

SITES = Channel.fromPath("$indir/*")
  .map{ sites_file -> tuple(sites_file.name.replaceAll(/\.tsv\.gz$/, ''),
    file(sites_file)) }

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

  input:
  tuple spec, file(sites_file) from SITES
  val q_thres from params.q_thres
  val min_patients from params.min_patients
  val iter from params.iter
  val warmup from params.warmup
  val chains from params.chains
  val vp from params.vp
  val vq from params.vq

  output:
  file "output/m1.stan.rdat" optional true
  file  "output/p_directional.tsv.gz" optional true
  file  "output/model_summaries.tsv.gz" optional true

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
  memory = '5G'
  withLabel: 'r'{
    module = 'R/4.1.0'
    // module = "R/4.0.2:v8/8.4.371.22" // Make sure you have ~/.R/Makevars with CXX14
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
