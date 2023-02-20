## module load r-bundle-bioconductor/3.8-foss-2018b-r-3.5.1

## singularity run docker://r-base

## singularity pull --arch amd64 library://alefrol/default/bioconductor_3.13:latest
## singularity run docker://r-base

## https://stackoverflow.com/questions/36651091/how-to-install-packages-in-linux-centos-without-root-user-with-automatic-depen
## yumdownloader --destdir ~/rpm --resolve libmpc-devel mpfr-devel gmp-devel
## for file in ~/rpm/*.rpm; do echo $file; cd ~/centos && rpm2cpio $file | cpio -id; cd -; done

## intsalling R.4.1.2 on the cluster
## cd ~/progs/; wget https://cran.r-project.org/src/base/R-4/R-4.1.2.tar.gz
## tar xf R-4.1.2.tar.gz
## cd R-4.1.2
## export  LDFLAGS=" -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5  -lgfortran "

## ./configure --prefix=$HOME/centos

#

## LDFLAGS=-L/users/anna.hakobyan/centos/usr/lib/gconv -L/users/anna.hakobyan/centos/usr/lib/audit -L/users/anna.hakobyan/centos/usr/lib/dracut -L/users/anna.hakobyan/centos/usr/lib/fipscheck -L/users/anna.hakobyan/centos/usr/lib/krb5 -L/users/anna.hakobyan/centos/usr/lib/nss -L/users/anna.hakobyan/centos/usr/lib/openssl -L/users/anna.hakobyan/centos/usr/lib/pkgconfig -L/users/anna.hakobyan/centos/usr/lib/sasl2 -L/users/anna.hakobyan/centos/usr/lib/tmpfiles2 -L/users/anna.hakobyan/centos/usr/lib/sse2 -L/users/anna.hakobyan/centos/usr/lib/gconv -L/users/anna.hakobyan/centos/usr/lib64/pkgconfig

## export CPATH=$HOME/centos/include:$CPATH


## rpm2cpio ./ncurses-libs-6.1-9.20180224.el8.x86_64.rpm | cpio -idmv

## ./configure --prefix=$HOME/centos  --enable-R-shlib --with-readline=no  CPPFLAGS="-I$HOME/centos/include/"  LDFLAGS="-L/users/anna.hakobyan/centos/lib  -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5  -lgfortran -L/users/anna.hakobyan/centos/usr/lib64/"


# building singularity https://gist.github.com/rohitfarmer/af486e346900534274cdc1153a764c90
library(TCGAbiolinks)
library(SummarizedExperiment)

## library(tidyverse)
library(tictoc)
GDC_projects = getGDCprojects()

tcga_project_ids = grep("TCGA-", GDC_projects$id, value = TRUE)

projectid = tcga_project_ids[1]
filename="runtime.dat"
for (projectid in tcga_project_ids[30:length(tcga_project_ids)]) {
    tic()
    system.time({proj_exps_query <- GDCquery(
        project = projectid,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification", 
        workflow.type = "STAR - Counts"
        #   barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
    )})
    
    
    system.time({
        GDCdownload(proj_exps_query)
        project_data <- GDCprepare(proj_exps_query)
    })
    
    project_assay = assay(project_data)
    write.table(project_assay, file = file.path("RNA_seq", paste0(projectid, ".star.counts.dat")),
                row.names = TRUE, col.names = NA, quote = FALSE,
                sep = "\t")
    
    tt = toc()
    write(paste(projectid, tt$toc, "\n"), file = filename, append = TRUE)
}



# cut -d $'\t' -f1 TCGA-GBM.star.counts.dat > tcga.all.exp.counts.dat;
# 
# for file in *.star.counts.dat; 
# do echo $file; 
# paste tcga.all.exp.counts.dat <(cut -d$'\t' -f2- "$file") >tcga.all.exp.counts.dat_
# mv tcga.all.exp.counts.dat_ tcga.all.exp.counts.dat
# done


