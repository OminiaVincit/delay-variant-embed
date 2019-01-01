@set DATA=ECGFiveDays
@set INPUT=%DATA%_time_series.txt
@set BIN=..\..\ph-compute\DelayVariantTopo\Release\

@set DIM=1
@set SCALE=1
@set MULTI=1
@set TAUMAX=20

@set EMD=5
@set DST=ph(m=%EMD%,s=%SCALE%)
%BIN%TimeSeriesRunner_D64.exe --nthreads -1 --emd %EMD% -i  %INPUT% -d %DIM% --multi %MULTI% -o %DST% --taumax %TAUMAX% --scale %SCALE%

