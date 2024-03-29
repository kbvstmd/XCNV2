if [ -f "./data/CDTS_percentile.txt" ];then
  filezise=`ls -l ./data/CDTS_percentile.txt | awk '{ print $5}'`
  if [ "$filezise" = 7644586817 ];then
        echo "file exists!"
        else
        rm -f ./data/CDTS_percentile.txt
    wget -c http://119.3.41.228:8080/XCNV/data/CDTS_percentile.txt.gz -O ./data/CDTS_percentile.txt.gz
    gunzip ./data/CDTS_percentile.txt.gz
  fi
  else
  wget http://119.3.41.228:8080/XCNV/data/CDTS_percentile.txt.gz -O ./data/CDTS_percentile.txt.gz
  gunzip ./data/CDTS_percentile.txt.gz
fi
if [ -f "./data/hg19_ljb26_all_converted.vcf" ];then
  filezise=`ls -l ./data/hg19_ljb26_all_converted.vcf | awk '{ print $5}'`
  if [ "$filezise" = 2242773397 ];then
        echo "file exists!"
        else
        rm -f ./data/hg19_ljb26_all_converted.vcf
        wget -c http://119.3.41.228:8080/XCNV/data/hg19_ljb26_all_converted.vcf.gz -O ./data/hg19_ljb26_all_converted.vcf.gz
        gunzip ./data/hg19_ljb26_all_converted.vcf.gz
  fi
  else
  wget http://119.3.41.228:8080/XCNV/data/hg19_ljb26_all_converted.vcf.gz -O ./data/hg19_ljb26_all_converted.vcf.gz
  gunzip ./data/hg19_ljb26_all_converted.vcf.gz
fi
rscript=`which Rscript`
script_path=`pwd`
tar zxvf ./tools/bedtools2.tar.gz -C tools
make
cp -p ./tools/bedtools2/bin/bedtools ./tools
${rscript} ./data/install_packages.R
echo -e "#!"${rscript}"\nscript.path=\""${script_path}"\"" | cat - ./data/XCNV.R > ./bin/XCNV
chmod 777 ./bin/XCNV
