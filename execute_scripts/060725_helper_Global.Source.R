# For the interractive session in Minerva
.libPaths(c("/hpc/packages/minerva-centos7/rpackages/4.2.0/site-library", 
            "/hpc/packages/minerva-centos7/rpackages/bioconductor/3.15", 
            .libPaths()))

# Load MultiWAS

libs <- .libPaths()
libs[3] <- "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"
.libPaths(libs)


# Load MultiWAS
library(MultiWAS)

global.source <- "/hpc/users/anyfam01/Global.Scripts/Global.Source.R"


path.to.scripts <- "/hpc/users/anyfam01/Global.Scripts"
global.scripts <- list.files(path.to.scripts, full.names = T)

# remove present script to avoid endless loop
global.scripts <- global.scripts[-grep("Global.Source.R", global.scripts)]

p.paths <- "/hpc/users/anyfam01/Global.Scripts/Global.Paths.R"
paths <- p.paths
source(paths)

p.global.functions <- "/hpc/users/anyfam01/Global.Scripts/Global.Functions.R"
global.functions <- p.global.functions
source(global.functions)


p.gsea.fisher <- "/hpc/users/anyfam01/Global.Scripts/GSEA_Fisher_OR_Functions.R"
p.camera.for.all <- "/hpc/users/anyfam01/Global.Scripts/source_camera.for.all.R"

load('/hpc/users/anyfam01/Global.Resources/viridis_pallete.RData')







