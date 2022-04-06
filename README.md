# Chroma analysis pipeline

This is a collection of code written for running my PhD analysis pipelines.

* [Shared libraries](#shared-libraries)
* [Chroma parameter generation](#chroma-parameter-generation)
* [Chroma result extraction](#chroma-result-extraction)
* [Choma result analysis](#chroma-result-analysis)

## Shared libraries

- [hl.h](hl.h): Single header library that houses all my common utilities and is used by the C code within this project. It includes:
	- Macros for assertions and utilising SSE/AVX intrinsics
	- Simple memory and program pool defitions
	- Matrix multiplication library for dealing with colour matrices using SIMD instructions
	- Small string library
	- A JSON construction and XML parsing library to fascilitate XML (of the Chroma of spec variety) to JSON files
	- The obligatory codebase hash map implementation
	- A GUI app base with dynamic loading of newly compiled code
	- A lot adapted following the excellent handmade hero

- [hl.py](hl.py): Common functions required across my analysis scripts, used mainly for handling bootstraps, concurrent maps and fits and nullspaces.

## Chroma parameter generation

Like every single chroma user I know, i ended up writting a script which writes the multi thousand line xml parameter input for me from a small parameter file. As part of the generation I also used these templates to create job scripts for running the simulation on various clusters, with their respective queue system parameters.

- [generate_chroma_jobs.py](chroma_parameter_generators/generate_chroma_jobs.py): Will scan for `params.yml` in the current directory which contains the propagators to generate, or trajectories to wilson flow and will create run scripts for the current HPC cluster detected as well as mustache templates using [chroma_parameter_generators/generate_chroma_input_xml.py][generate_chroma_input_xml.py] to input the sample and propagator sources for a fully functional chroma parameter file. These jobs can be used in combination with a `filelist` containing paths to all the trajectories to run the simulation on.
- [generate_template_for_next_trajectory.py](chroma_parameter_generators/generate_template_for_next_trajectory.py): Is called by an active job to determine the next sample to run the simulation on, and create the final chroma parameter file for it. 
- There is also [wilson_plaq_job.sh](chroma_parameter_generators/wilson_plaq_job.sh), [wilson_plaq_runner.sh](chroma_parameter_generators/wilson_plaq_runner.sh), and [merge_plaq_results.py](chroma_parameter_generators/merge_plaq_results.py) for calculating wilson plaquttes for the trajectories for downstram analysis.

## Chroma result extraction 

I used a small utility script to exit the chroma output formats as quickly as possible. In combination this is an extraction of the run parameters from the chroma psuedo-xml output and real matricies. The xml is parsed and written in JSON, with a AoS to SoA transformation grouping by parameter not by sample. In addition, all momentum and propagators get converted to human readable format. As a result the downstream analysis becomes much more staightforward to write. The compile script creats multiple executables with hard coded momentum extraction modes.

## Chroma result analysis

The first step after extracting the data is to perform coherent bootstraps followed by effectivemass ratios; the somewhat more expensive but easily cacheable components of the calculation.

- [generate_feynhell_analysis.py](chroma_result_analysis/generate_feynhell_analysis.py): Creates the pre-processing scripts for Feynhellman propagators extracted chroma runs, mainly in charge of generating coherent bootstraps across different momentum transfers
- [generate_gluon_run.py](chroma_result_analysis/generate_gluon_run.py): Creates the pre-processing script for gluon verteces.

These are made up of several components:

- [hl_bootstrap](chroma_result_analysis/hl_bootstrap): Create and save coherent bootstraps
- [hl_effmass](chroma_result_analysis/hl_effmass): Generate effective masses, including all the critical ratios, and all the threepoint functions
- [hl_fit](chroma_result_analysis/hl_fit): Fit gluon threepoint functions and extract form factors
- [hl_threepoint_gen](chroma_result_analysis/hl_threepoint_gen): Generate gluon threepoint functions from two point and action
- [hl_plot.py](chroma_result_analysis/hl_plot.py): Create plot yaml files that can be turned into plots using [hl_plot](chroma_result_analysis/hl_plotz). This allows the generation of whole new plot asthetics for my entire thesis in record time.

Finally this allows the extraction of the compton amplitudes using:
- [extract_compton_amplitude](chroma_result_analysis/analysis_codes/extract_compton_amplitude): Does the final analysis combining the effectivemass ratios across all the simulations, producing the core results of my thesis.
