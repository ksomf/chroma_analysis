# Chroma analysis pipeline

This is a collection of code written for running my PhD analysis pipelines.

* [Shared libraries](#shared-libraries)
* [Chroma parameter generation](#chroma-parameter-generation)

## Shared libraries

- [hl.h][hl.h]: Single header library that houses all my common utilities and is used by the C code within this project. It includes:
	- Macros for assertions and utilising SSE/AVX intrinsics
	- Simple memory and program pool defitions
	- Matrix multiplication library for dealing with colour matrices using SIMD instructions
	- Small string library
	- A JSON construction and XML parsing library to fascilitate XML (of the Chroma of spec variety) to JSON files
	- The obligatory codebase hash map implementation
	- A GUI app base with dynamic loading of newly compiled code

- [hl.py][hl.py]: Common functions required across my analysis scripts, used mainly for handling bootstraps, concurrent maps and fits and nullspaces.

## Chroma parameter generation

Like every single chroma user I know, i ended up writting a script which writes the multi thousand line xml parameter input for me from a small parameter file. As part of the generation I also used these templates to create job scripts for running the simulation on various clusters, with their respective queue system parameters.

- [chroma_parameter_generators/generate_chroma_jobs.py][chroma_parameter_generators/generate_chroma_jobs.py]: Will scan for `params.yml` in the current directory which contains the propagators to generate, or trajectories to wilson flow and will create run scripts for the current HPC cluster detected as well as mustache templates using [chroma_parameter_generators/generate_chroma_input_xml.py][generate_chroma_input_xml.py] to input the sample and propagator sources for a fully functional chroma parameter file. These jobs can be used in combination with a `filelist` containing paths to all the trajectories to run the simulation on.
- [chroma_parameter_generators/generate_template_for_next_trajectory.py][generate_template_for_next_trajectory.py]: Is called by an active job to determine the next sample to run the simulation on, and create the final chroma parameter file for it. 
- There is also [chroma_parameter_generators/wilson_plaq_job.sh][wilson_plaq_job.sh], [chroma_parameter_generators/wilson_plaq_runner.sh][wilson_plaq_runner.sh], and [chroma_parameter_generators/merge_plaq_results.py][merge_plaq_results.py] for calculating wilson plaquttes for the trajectories for downstram analysis.
