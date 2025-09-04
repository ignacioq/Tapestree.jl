using Documenter, Tapestree

makedocs(
         sitename = "Tapestree.jl",
         modules  = [Tapestree],
         format   = Documenter.HTML(
            assets = ["assets/favicon.ico"],
            sidebar_sitename = false
         ),
         pages=[
                "Home" => "index.md"
                "Installation" => "installation.md"
                "Quick start" => "quick_start.md"
                "Manual" => [
                    "INSANE" => [
                        "Contents" => "man/insane/contents.md",
                        "Input and structures" => "man/insane/io.md",
                        "Constant birth-death" => "man/insane/cbd.md",
                        "Birth-death diffusion" => "man/insane/bdd.md",
                        "Constant fossilized birth-death" => "man/insane/cfbd.md",
                        "Fossilized birth-death diffusion" => "man/insane/fbdd.md",
                        "Occurrence birth-death" => "man/insane/obd.md",
                        "Diffused Brownian motion" => "man/insane/dbm.md",
                        "Processing" => "man/insane/processing.md",
                        "Plotting" => "man/insane/iplots.md"
                        ],
                    "TRIBE" => "man/tribe.md",
                    "ESSE" => "man/esse.md",
                    ]
               ],
        authors   = "Ignacio Quintero",
        checkdocs = :none)

deploydocs(
    repo="github.com/ignacioq/Tapestree.jl.git",
)