


#=======================================================================
# Download HIPPIE
#=======================================================================
mkdir -p HIPPIE
wget -P HIPPIE http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt


#=======================================================================
# TADs from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014

RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK"


for CELL in ${RAO_CELLS} ; do
  
  # download
  wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
  
  # unzip 
  gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
      
  tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
    |cut -f 1-3 \
    | sed -e 's/^/chr/' \
    > Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 

done

