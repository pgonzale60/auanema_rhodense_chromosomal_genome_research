library(dplyr)
library(ggplot2)
library(HiCExperiment)
library(HiContacts)
library(rtracklayer)



hic <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = 'mcool', resolution = 32000)
hic_X <- import("analyses/hi-C/cooltools/nxAuaRhod1_1.mcool", format = 'mcool', resolution = 32000, focus = 'SUPER_X')


fig_hic <- patchwork::wrap_plots(
  plotMatrix(hic),
  plotMatrix(hic_X, show_grid = T #, borders = GRS
             )
)


ggsave("report/figures/nxAuaRhod1.1_hic_panels.pdf",
       fig_hic, width = 10, height = 6)
