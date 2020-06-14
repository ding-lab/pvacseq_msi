FROM ensemblorg/ensembl-vep:release_97.4

LABEL maintainer="Wen-Wei Liang <liang.w@wustl.edu>"

# Install plugins
RUN perl ~/src/ensembl-vep/INSTALL.pl -n -a p -g Downstream

RUN cp /opt/vep/.vep/Plugins/Downstream.pm ~/src/ensembl-vep/modules/

RUN curl https://raw.githubusercontent.com/griffithlab/pVACtools/master/tools/pvacseq/VEP_plugins/Wildtype.pm -o ~/src/ensembl-vep/modules/Wildtype.pm


