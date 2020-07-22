#!/bin/bash

lualatex -shell-escape PWDFT_docs.tex
bibtex PWDFT_docs
makeindex PWDFT_docs
lualatex -shell-escape PWDFT_docs.tex
