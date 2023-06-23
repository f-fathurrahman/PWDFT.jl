#!/usr/bin/env julia
# From: https://github.com/JuliaCI/Coverage.jl

using Coverage
using Printf

# process '*.cov' files
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
#coverage = append!(coverage, process_folder("deps"))  # useful if you want to analyze more than just src/

# process '*.info' files, if you collected them
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = ( joinpath(pwd(), "src", ""),
                    joinpath(pwd(), "deps", "") )
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))

# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)

println()
@printf("covered_lines = %d\n", covered_lines)
@printf("total_lines   = %d\n", total_lines)
@printf("Percentage    = %.1f%%\n", covered_lines/total_lines * 100)
