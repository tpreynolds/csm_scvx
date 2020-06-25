#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=# 

using Debugger
break_on(:error)
include("utils.jl")
# include("types.jl")
include("scvx.jl")
# include("convexify.jl")
include("hello.jl")