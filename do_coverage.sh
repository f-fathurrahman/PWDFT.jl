#!/bin/bash
julia --code-coverage=user --inline=no test/runtests.jl
