@set DATA=ECGFiveDays
@set BIN=..\..\ph-compute\DelayVariantTopo\Release\
@set DST=kernel

@set SCALE=1
@set EMD=5
@set DIM=1
@set OPTTAU=0

@set TIMEHOLE=0.0
@set TIMETAU=1.0
@set ALLTAU=1

@set POSFIX=multi_scale
@set THRES=0.0
@set METHOD=0
@set TAUMAX=0

@set LEFT=%DATA%_(m=%EMD%,s=%SCALE%)_barcodes_dim_%DIM%.txt
@set OUT=%DST%\%DATA%_(m=%EMD%,s=%SCALE%)_dim_%DIM%_kernel_method_%METHOD%_timehole_%TIMEHOLE%_timetau_%TIMETAU%_thres_%THRES%_alltau_%ALLTAU%_taumax_%TAUMAX%_opttau_%OPTTAU%_%POSFIX%.txt

%BIN%DiagramDistanceRunner_D64.exe --dim %DIM% --timehole %TIMEHOLE% --timetau %TIMETAU% --alltau %ALLTAU% --opttau %OPTTAU% -l %LEFT% -o %OUT% --threshold %THRES% --method %METHOD% --taumax %TAUMAX%

