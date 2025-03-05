using Pkg; Pkg.activate(".")
using HDF5
using LLRParsing

base_dir = "/home/fabian/Downloads/llr_parser_test_data/su3_4x20_8/"
repeats, replica_dirs = get_repeat_and_replica_dirs(base_dir)