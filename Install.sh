wget http://119.3.41.228:8080/XCNV/data/CDTS_percentile.txt.gz -O ./data/CDTS_percentile.txt.gz
gunzip ./data/CDTS_percentile.txt.gz
wget http://119.3.41.228:8080/XCNV/data/hg19_ljb26_all_converted.vcf.gz -O ./data/CDTS_percentile.txt.gz
gunzip ./data/CDTS_percentile.txt.gz
rscript=`which Rscript`
script_path=`pwd`
tar zxvf ./tools/bedtools2.tar.gz -C tools
make
cp -p ./tools/bedtools2/bin/bedtools ./tools
${rscript} ./data/install_packages.R
echo -e "#!"${rscript}"\nscript.path=\""${script_path}"\"" | cat - ./data/XCNV.R > ./bin/XCNV
chmod 777 ./bin/XCNV
