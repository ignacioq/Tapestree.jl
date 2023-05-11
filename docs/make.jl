using Documenter, Tapestree

makedocs(
         sitename = "Tapestree.jl",
         modules  = [Tapestree],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
    repo="github.com/ignacioq/Tapestree.jl",
)