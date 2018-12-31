@set DATA=ECGFiveDays
@set INPUT=%DATA%_time_series.txt
@set BIN=..\..\ph-compute\DelayVariantTopo\Release\

@set NPROCS=16
@set DIM=1
@set SCALE=1
@set MULTI=1
@set TAUMAX=20

@set EMD=4
@set DST=ph(m=%EMD%,s=%SCALE%)
%BIN%TimeSeriesRunner_D64.exe --nprocs %NPROCS% --emd %EMD% -i  %INPUT% -d %DIM% --multi %MULTI% -o %DST% --taumax %TAUMAX% --scale %SCALE%

