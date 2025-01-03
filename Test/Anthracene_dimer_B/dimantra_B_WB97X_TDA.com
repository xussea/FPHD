%chk=dimantra_B_WB97X_TDA
%rwf=dimantra_B_WB97X_TDA
%Mem=6Gb
%Nproc=6
#P WB97XD 6-31G** nosymm scf=direct gfinput gfoldprint iop(6/7=3) TDA=(nstates=60,singlets)

A

0 1
C      -1.569991    6.189055    3.351213
C      -0.767361    6.951257    2.573187
C      -0.374166    6.507291    1.272799
C       0.445777    7.263477    0.438674
C      -0.827504    5.225324    0.824008
C      -1.674837    4.465528    1.682964
C      -2.031882    4.929347    2.902422
H      -1.825398    6.509096    4.221204
H      -0.428882    7.802493    2.878512
H       0.771291    8.121330    0.754115
H      -1.984436    3.579401    1.351889
H      -2.641377    4.385518    3.494679
C       0.827504    6.806276   -0.824008
C      -0.445777    4.768123   -0.438674
C       0.374166    5.524309   -1.272799
C       1.674837    7.566072   -1.682964
H      -0.771291    3.910270   -0.754115
C       0.767361    5.080343   -2.573187
C       2.031882    7.102253   -2.902422
H       1.984436    8.452199   -1.351889
C       1.569991    5.842545   -3.351213
H       0.428882    4.229107   -2.878512
H       2.641377    7.646082   -3.494679
H       1.825398    5.522504   -4.221204
C       5.846291    3.181155   -3.351213
C       5.043661    3.943357   -2.573187
C       4.650466    3.499391   -1.272799
C       3.830523    4.255577   -0.438674
C       5.103804    2.217424   -0.824008
C       5.951137    1.457628   -1.682964
C       6.308182    1.921447   -2.902422
H       6.101698    3.501196   -4.221204
H       4.705182    4.794593   -2.878512
H       3.505009    5.113430   -0.754115
H       6.260736    0.571501   -1.351889
H       6.917677    1.377618   -3.494679
C       3.448796    3.798376    0.824008
C       4.722077    1.760223    0.438674
C       3.902134    2.516409    1.272799
C       2.601463    4.558172    1.682964
H       5.047591    0.902370    0.754115
C       3.508939    2.072443    2.573187
C       2.244418    4.094353    2.902422
H       2.291864    5.444299    1.351889
C       2.706309    2.834645    3.351213
H       3.847418    1.221207    2.878512
H       1.634923    4.638182    3.494679
H       2.450902    2.514604    4.221204
