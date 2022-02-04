# Chroma analysis pipeline

This is a collection of code written for running my PhD analysis pipelines.

* [Shared libraries](#shared-libraries)

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
