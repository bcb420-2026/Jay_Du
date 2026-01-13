FROM risserlin/bcb420-base-image:winter2026

# Install BiocManager if not available
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

# Install required packages (DESeq2, pheatmap, enrichplot)
RUN R -e 'BiocManager::install(c("DESeq2", "enrichplot"))'
RUN R -e 'install.packages("pheatmap")'

# Set working directory
WORKDIR /home/rstudio/Jay_Du
