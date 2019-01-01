@set SRC=data\
@set BIN=Release\PersistentRunner_D64.exe


%BIN% --nthreads 16 --maxdim 1 --format point-cloud --outdir output --input %SRC%pcl_1.txt --multi 0

%BIN% --nthreads -1 --maxdim 1 --format point-cloud --outdir output --input %SRC%pcl_list.txt --multi 1