using Documenter, Tapestree

makedocs(
         sitename = "Tapestree.jl",
         modules  = [Tapestree],
         pages=[
                "Home" => "index.md"
                "Manual" => [
                    "Installation" => "man/installation.md",
                    "INSANE" => "man/insane.md",
                    "TRIBE" => "man/tribe.md",
                    "ESSE" => "man/esse.md",
                    ]
               ],
        authors="Ignacio Quintero")

deploydocs(;
    repo="github.com/ignacioq/Tapestree.jl",
    push_preview = true,
)