FROM neurogenomicslab/mungesumstats:1.7.10

# Install necessary databases (GRCh38)
RUN Rscript -e "options(timeout = 10000); BiocManager::install('SNPlocs.Hsapiens.dbSNP155.GRCh38')" && \
    Rscript -e "options(timeout = 10000); BiocManager::install('BSgenome.Hsapiens.NCBI.GRCh38')"

# Install necessary databases (GRCh37)
RUN Rscript -e "options(timeout = 10000); BiocManager::install('SNPlocs.Hsapiens.dbSNP155.GRCh37')" && \
    Rscript -e "options(timeout = 10000); BiocManager::install('BSgenome.Hsapiens.1000genomes.hs37d5')"

# Install vroom
# Install necessary databases (GRCh37)
RUN Rscript -e "install.packages('vroom')"

WORKDIR /home/
CMD ["bash"]