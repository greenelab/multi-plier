FROM jtaroni/multi-plier:recount

RUN R -e "devtools::install_github('ebecht/MCPcounter', ref = 'a79614eee002c88c64725d69140c7653e7c379b4', subdir = 'Source', dependencies = TRUE)"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/VennDiagram_1.6.20.tar.gz')"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/cowplot/cowplot_0.9.2.tar.gz')"
RUN R -e "devtools::install_github('topepo/caret/pkg/caret', ref = '6546939345fe10649cefcbfee55d58fb682bc902')"
RUN R -e "devtools::install_version('e1071', version = '1.6-8')"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/ggsignif_0.4.0.tar.gz')"
RUN R -e "BiocInstaller::biocLite('SRAdb')"
