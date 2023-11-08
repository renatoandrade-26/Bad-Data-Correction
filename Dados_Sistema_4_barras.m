%% Dados de entrada do sistema de distribuição com 4 barras

%% Dados gerais do sistema

% base de tensão V = 1.0 kV; S = 100MVA

Sbase = 100; % MVA
Vbase = 1.0; % kV

K = (Vbase^2)/Sbase;

nbus = 4;


%% Dados de barras

% Tipos de barras:  1 - Barra Slack; 
%                   2 - Barra PV;
%                   3 - Barra PQ;

% OBS - Potências em MW

%          |Bus | Type |    Vsp |   theta |   PGi |    QGi |    PLi    |   QLi   |  Qmin | Qmax | Shunt 
busdata = [   1     1       1.000     0      0.0000    0.000    0          0        0       0     0;     
              2     3       1.000     0      0.0000    0.000    2.00       1.00     0       0     0;
              3     3       1.000     0      0.0000    0.000    4.00       1.00     0       0     0;   
              4     2       0.980     0      9.0000    0.000    4.00       0.00     0       0     0;
              ];
          
%% Dados de linhas

%         |  From |  To   |   R        |   X     |    B/2  |  X'mer  |   | Defasagem (°) |
%         |  Bus  | Bus   |  pu        |  pu     |    pu   | TAP (a) |   | Defasagem (°) |
linedata =  [1      2       0.2      	0.1          0        1             0
             2      3       0.2     	0.1          0        1             0             
             2      4       0.1 	    0.05         0        1             0
             ];