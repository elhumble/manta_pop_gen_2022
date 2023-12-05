# Script to collate and plot figures for manuscript
 
library(patchwork)
library(grDevices)
library(extrafont)
library(ggplot2)
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Figure 1          #         
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

dapc_alf <- readRDS("figs/DAPC_alf.RDS")
dapc_bir <- readRDS("figs/DAPC_bir.RDS")

m_alf <- readRDS("figs/map_alfredi.RDS")
m_bir <- readRDS("figs/map_birostris.RDS")

ggsave("figs/Figure_1.png", m_alf + dapc_alf + m_bir + dapc_bir + 
         plot_layout(ncol = 2, 
                     heights = c(1, 1, 1, 1), 
                     widths = c(2, 1, 2, 1)), 
       width = 10, height = 12, dpi = 600)

ggsave("figs/Figure_1.pdf", m_alf + dapc_alf + m_bir + dapc_bir + 
         plot_layout(ncol = 2, 
                     heights = c(1, 1, 1, 1), 
                     widths = c(2, 1, 2, 1)), 
       width = 10, height = 12, dev = cairo_pdf)


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Figure 2          #         
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

#install.packages("ggplotify")

fst_alf <- readRDS("figs/fst_alf.RDS")
fst_bir <- readRDS("figs/fst_bir.RDS")

ibd_alf <- readRDS("figs/ibd_alf.RDS")
ibd_bir <- readRDS("figs/ibd_bir.RDS")

fst_alf + ibd_alf + 
  fst_bir + ibd_bir + 
  plot_layout(ncol = 2,
              widths = c(1, 1.5, 1, 1.5))

ggsave("figs/Figure_2.png", fst_alf + ibd_alf + 
         fst_bir + ibd_bir + 
         plot_layout(ncol = 2,
                     widths = c(1, 1.5, 1, 1.5)), 
       width = 7, height = 6)

ggsave("figs/Figure_2.pdf", fst_alf + ibd_alf + 
         fst_bir + ibd_bir + 
         plot_layout(ncol = 2,
                     widths = c(1, 1.5, 1, 1.5)), 
       width = 7, height = 6, dev = cairo_pdf)
 