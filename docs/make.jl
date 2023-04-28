using MoM_Basics
using Documenter

DocMeta.setdocmeta!(MoM_Basics, :DocTestSetup, :(using MoM_Basics); recursive=true)

makedocs(;
    modules=[MoM_Basics],
    authors="deltaeecs <1225385871@qq.com> and contributors",
    repo="https://github.com/deltaeecs/MoM_Basics.jl/blob/{commit}{path}#{line}",
    sitename="MoM_Basics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://deltaeecs.github.io/MoM_Basics.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/deltaeecs/MoM_Basics.jl",
)
