# delay-variant-embed
Analysis code for "Topological time-series analysis with delay-variant embedding" manuscript.

Preprint at https://arxiv.org/abs/1803.00208

Authors:

* Quoc Hoan Tran: zoro@biom.t.u-tokyo.ac.jp
* Yoshihiko Hasegawa: hasegawa@biom.t.u-tokyo.ac.jp

Hasegawa Lab, Department of Information and Communication Engineering,
Graduate School of Information Science and Technology,
The University of Tokyo

## Compute three-dimensional persistence diagrams and persistence kernels for time-series data

### Compile by Visual Studio 2015
Download ph-compute folder and launch ```DelayVariantTopo\TopolologyAnalysis.sln```

### Release/PersistenceRunner_D64.exe
Calculate persistence diagrams from point clouds.

    Options:
        --nthreads : number of threads using for multi-processing (-1 for all possible threads)
        --modulus  : coefficient in the prime field Z/<p>Z to compute homology
        --maxdim   : maximum holes dimension to compute homology
        --thres    : maximum diameter to construct Rips complex
        --format   : format of the input (point-cloud, lower-distance, upper-distance, distance, dipha)
        --outdir   : directory for the output
        --input    : input file or file contains list of input files
        --multi    : compute single input file (0) or multi files in one input file  (1) 
        --help     : print usage option

    Example: (see more at ph-compute\DelayVariantTopo\run-compute-persistent.bat)

        Release\PersistentRunner_D64.exe --nthreads 16 --maxdim 1 --format point-cloud --outdir output --input data\pcl_1.txt --multi 0

        Release\PersistentRunner_D64.exe --nthreads -1 --maxdim 1 --format point-cloud --outdir output --input data\pcl_list.txt --multi 1

### Release/TimeSeriesRunner_D64.exe
Calculate 3-dimensional persistence diagrams for time-series data via delay-variant embedding.

    Options:
        --nthreads : number of threads using for multi-processing (-1 for all possible threads)
        --emd      : embedded dimension for delay-variant embedding
        --modulus  : coefficient in the prime field Z/<p>Z to compute homology
        --maxdim   : maximum holes dimension to compute homology
        --thres    : maximum diameter to construct Rips complex
        --outdir   : directory for the output
        --input    : input time-series file or folder
        --multi    : compute single series (0) or multi series in one input file  (1) 
        --taumax   : compute up to maximum value of tau
        --npoints  : embed up to number of points, default=0, normal delay embedding
        --scale    : downsampling scale for time series, default=1
        --nskip    : skip a number of time series
        --help     : print usage option

Example of time-series data for input (csv format):

    0.0099134,0.033429,-0.0057636,-0.029279,-0.052795,0.0099134,-0.021441,-0.029279,-0.060633,-0.091987
    -0.013602,-0.060633,-0.013602,-0.07631,-0.060633,-0.12334,-0.084149,-0.068472,-0.037118,-0.07631
    -0.060633,-0.12334,-0.060633,-0.091987,-0.068472,-0.10766,-0.099826,-0.14686,-0.029279,0.0099134
    
For multi mode (multi=1), the programm will process one line as one time series (see more at timeseries/ECGFiveDays/run-compute-delay-PH.bat).

### Release/DiagramDistanceRunner_D64.exe
Calculate kernel for 3-dimensional persistence diagrams.

    Options:
        --timehold : parameter sigmal in the kernel (=0 for optimal value)
        --timetau  : ratio of xi vs. sigmal in the kernel
        --left     : (left) input as list of barcodes
        --right    : (right) input as list of barcodes
        --dim      : dimension of holes to compute kernel
        --skipinf  : skip holes which have infinity death-scale (default=True)
        --infval   : replace infinity death-scales with a default value
        --thres    : threshold to skip holes with death-birth < thres (default=0 to use all holes)
        --output   : output filename for gram matrix
        --posfix   : posfix for output file
        --method   : method to compute kernel
                     0: L2_inner_multi_sse, 1: L2_squared_distance, 
                     2: Slice Wasserstein, 3:L2_inner_multi_nosse, 4:riemmannian_metric
        --theta_ndirs   : number of direction in wasserstein slice distance for theta
        --phi_ndirs     : number of direction in wasserstein slice distance for phi
        --alltau   : use all tau (True), or single tau (False)
        --opttau   : specify single tau (enable if alltau=0)
                     if opttau > 0  -> calculate single-delay kernel with tau=opttau
                     if opttau <= 0 -> calculate single-delay kernel with tau specified in barcode list file
        --help     : print usage option

Example of list of barcodes for input:

    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_0_dim_1.txt;tau=8
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_1_dim_1.txt;tau=9
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_2_dim_1.txt;tau=6
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_3_dim_1.txt;tau=9
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_4_dim_1.txt;tau=7
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_5_dim_1.txt;tau=9
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_6_dim_1.txt;tau=8
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_7_dim_1.txt;tau=16
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_8_dim_1.txt;tau=7
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_9_dim_1.txt;tau=3
    ph(m=5,s=1)/delay_barcode_ECGFiveDays_time_series_10_dim_1.txt;tau=13
