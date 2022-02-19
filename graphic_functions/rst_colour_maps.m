function [RGB]=rst_colour_maps(N_colours)

% This is a simple utility to return colour maps as RGB
% For N_colours < 12, colours are taken from ColorBrewer.org
% http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
% Using the scheme that gives the 1st four colour as colour bind safe
% For N_colours > 12, it samplea at equidistance from the xrain colour 
% scale taken from https://github.com/CPernet/brain_colours
%
% FORMAT RGB = rst_colour_maps(N_colours)
% INPUT        N_colours, how many triplets are needed
% OUTPUT       RGB a N*3 matrix of colour triplets
%
% Cyril Pernet - RST toolbox
% ---------------------------

if N_colours <= 11
    
    % colorbrewer = [166  206 227
     colorbrewer = [31 120 180
        178 223 138
        51  160 44
        251 154 153
        227 26  28
        253 191 111
        255 127 0
        202 178 214
        106 61  154
        255 255 153
        177 89  40];
    
    colorbrewer = colorbrewer/255;
    RGB = colorbrewer(1:N_colours,:);
    
else
    
    xrain = [0    0.0045         0
        0.0167    0.0055    0.0044
        0.0377    0.0066    0.0424
        0.0563    0.0079    0.0715
        0.0711    0.0093    0.0927
        0.0832    0.0110    0.1107
        0.0924    0.0133    0.1282
        0.0990    0.0159    0.1456
        0.1040    0.0188    0.1628
        0.1089    0.0210    0.1797
        0.1139    0.0225    0.1963
        0.1192    0.0233    0.2126
        0.1245    0.0233    0.2286
        0.1301    0.0227    0.2443
        0.1357    0.0213    0.2598
        0.1414    0.0194    0.2749
        0.1470    0.0175    0.2898
        0.1524    0.0157    0.3045
        0.1577    0.0139    0.3190
        0.1629    0.0123    0.3334
        0.1678    0.0104    0.3478
        0.1725    0.0090    0.3621
        0.1770    0.0076    0.3765
        0.1812    0.0063    0.3909
        0.1852    0.0052    0.4055
        0.1889    0.0041    0.4202
        0.1922    0.0032    0.4352
        0.1951    0.0024    0.4504
        0.1977    0.0017    0.4658
        0.1997    0.0012    0.4815
        0.2014    0.0009    0.4974
        0.2025    0.0007    0.5136
        0.2032    0.0007    0.5300
        0.2033    0.0009    0.5466
        0.2029    0.0015    0.5633
        0.2019    0.0023    0.5802
        0.2003    0.0036    0.5970
        0.1982    0.0054    0.6138
        0.1955    0.0077    0.6305
        0.1922    0.0107    0.6470
        0.1885    0.0148    0.6631
        0.1842    0.0197    0.6788
        0.1794    0.0257    0.6939
        0.1743    0.0330    0.7083
        0.1689    0.0420    0.7220
        0.1632    0.0515    0.7348
        0.1573    0.0614    0.7465
        0.1514    0.0718    0.7572
        0.1456    0.0824    0.7666
        0.1398    0.0933    0.7748
        0.1345    0.1045    0.7816
        0.1294    0.1160    0.7871
        0.1249    0.1276    0.7911
        0.1209    0.1394    0.7937
        0.1176    0.1514    0.7949
        0.1149    0.1634    0.7947
        0.1129    0.1754    0.7931
        0.1115    0.1874    0.7902
        0.1106    0.1993    0.7861
        0.1102    0.2113    0.7809
        0.1101    0.2230    0.7747
        0.1102    0.2347    0.7676
        0.1105    0.2462    0.7597
        0.1109    0.2575    0.7511
        0.1112    0.2686    0.7419
        0.1115    0.2796    0.7323
        0.1116    0.2903    0.7223
        0.1116    0.3008    0.7121
        0.1114    0.3111    0.7017
        0.1110    0.3212    0.6913
        0.1104    0.3311    0.6808
        0.1096    0.3408    0.6704
        0.1086    0.3503    0.6601
        0.1075    0.3596    0.6499
        0.1063    0.3688    0.6399
        0.1049    0.3778    0.6300
        0.1035    0.3866    0.6204
        0.1019    0.3953    0.6109
        0.1003    0.4039    0.6017
        0.0987    0.4123    0.5927
        0.0971    0.4206    0.5838
        0.0955    0.4289    0.5751
        0.0939    0.4370    0.5666
        0.0923    0.4450    0.5582
        0.0907    0.4529    0.5500
        0.0891    0.4608    0.5420
        0.0876    0.4686    0.5340
        0.0859    0.4763    0.5262
        0.0845    0.4839    0.5184
        0.0829    0.4915    0.5108
        0.0815    0.4991    0.5032
        0.0800    0.5065    0.4958
        0.0786    0.5140    0.4884
        0.0772    0.5214    0.4811
        0.0759    0.5288    0.4738
        0.0745    0.5361    0.4666
        0.0731    0.5434    0.4594
        0.0717    0.5507    0.4523
        0.0702    0.5579    0.4452
        0.0689    0.5651    0.4381
        0.0673    0.5723    0.4311
        0.0658    0.5795    0.4241
        0.0643    0.5866    0.4172
        0.0629    0.5938    0.4103
        0.0612    0.6009    0.4034
        0.0598    0.6080    0.3965
        0.0581    0.6151    0.3896
        0.0566    0.6222    0.3828
        0.0549    0.6292    0.3760
        0.0532    0.6363    0.3692
        0.0516    0.6433    0.3625
        0.0499    0.6503    0.3557
        0.0482    0.6574    0.3490
        0.0465    0.6644    0.3423
        0.0448    0.6714    0.3356
        0.0429    0.6784    0.3290
        0.0413    0.6854    0.3223
        0.0395    0.6924    0.3157
        0.0377    0.6993    0.3091
        0.0360    0.7063    0.3024
        0.0341    0.7133    0.2959
        0.0324    0.7203    0.2893
        0.0308    0.7273    0.2828
        0.0292    0.7342    0.2762
        0.0276    0.7412    0.2697
        0.0261    0.7482    0.2632
        0.0246    0.7551    0.2567
        0.0232    0.7621    0.2502
        0.0218    0.7691    0.2437
        0.0205    0.7760    0.2373
        0.0193    0.7830    0.2308
        0.0183    0.7900    0.2243
        0.0176    0.7969    0.2179
        0.0173    0.8039    0.2115
        0.0176    0.8109    0.2051
        0.0186    0.8178    0.1987
        0.0204    0.8247    0.1924
        0.0233    0.8317    0.1860
        0.0276    0.8386    0.1796
        0.0336    0.8455    0.1733
        0.0421    0.8524    0.1670
        0.0518    0.8592    0.1607
        0.0630    0.8661    0.1544
        0.0753    0.8729    0.1481
        0.0890    0.8796    0.1418
        0.1037    0.8863    0.1356
        0.1197    0.8929    0.1293
        0.1368    0.8994    0.1231
        0.1550    0.9058    0.1169
        0.1742    0.9121    0.1107
        0.1945    0.9183    0.1044
        0.2157    0.9244    0.0983
        0.2378    0.9303    0.0921
        0.2609    0.9360    0.0859
        0.2849    0.9415    0.0797
        0.3096    0.9467    0.0736
        0.3351    0.9517    0.0674
        0.3612    0.9565    0.0612
        0.3879    0.9609    0.0551
        0.4152    0.9651    0.0490
        0.4428    0.9689    0.0427
        0.4707    0.9723    0.0367
        0.4987    0.9754    0.0309
        0.5269    0.9780    0.0258
        0.5550    0.9802    0.0212
        0.5829    0.9820    0.0171
        0.6106    0.9833    0.0135
        0.6378    0.9841    0.0101
        0.6645    0.9845    0.0073
        0.6906    0.9843    0.0050
        0.7159    0.9836    0.0030
        0.7403    0.9825    0.0013
        0.7638    0.9808         0
        0.7862    0.9786         0
        0.8075    0.9759         0
        0.8275    0.9728         0
        0.8464    0.9691         0
        0.8640    0.9650         0
        0.8803    0.9605         0
        0.8953    0.9556         0
        0.9090    0.9502         0
        0.9214    0.9445         0
        0.9327    0.9385         0
        0.9427    0.9321         0
        0.9517    0.9255         0
        0.9596    0.9185         0
        0.9665    0.9113         0
        0.9725    0.9039         0
        0.9777    0.8963         0
        0.9821    0.8885         0
        0.9859    0.8806         0
        0.9891    0.8725         0
        0.9917    0.8642         0
        0.9939    0.8559         0
        0.9957    0.8475         0
        0.9972    0.8390         0
        0.9984    0.8304         0
        0.9993    0.8217         0
        1.0000    0.8130         0
        1.0000    0.8043         0
        1.0000    0.7954         0
        1.0000    0.7866         0
        1.0000    0.7776         0
        1.0000    0.7687         0
        1.0000    0.7597         0
        1.0000    0.7506         0
        1.0000    0.7415         0
        1.0000    0.7324         0
        1.0000    0.7232         0
        1.0000    0.7140         0
        1.0000    0.7047         0
        1.0000    0.6954         0
        1.0000    0.6860         0
        1.0000    0.6766         0
        1.0000    0.6671         0
        1.0000    0.6576         0
        1.0000    0.6479         0
        1.0000    0.6383         0
        1.0000    0.6285         0
        1.0000    0.6187         0
        1.0000    0.6088         0
        1.0000    0.5988         0
        1.0000    0.5887         0
        1.0000    0.5786         0
        1.0000    0.5683         0
        1.0000    0.5579         0
        1.0000    0.5474         0
        1.0000    0.5368         0
        1.0000    0.5261         0
        1.0000    0.5152         0
        1.0000    0.5042         0
        1.0000    0.4930         0
        1.0000    0.4816         0
        1.0000    0.4701         0
        1.0000    0.4583         0
        1.0000    0.4464         0
        1.0000    0.4342         0
        1.0000    0.4217         0
        1.0000    0.4089         0
        1.0000    0.3959         0
        1.0000    0.3824         0
        1.0000    0.3686         0
        1.0000    0.3544         0
        1.0000    0.3396         0
        1.0000    0.3243         0
        1.0000    0.3083         0
        1.0000    0.2915         0
        1.0000    0.2738         0
        1.0000    0.2549         0
        1.0000    0.2347         0
        1.0000    0.2126         0
        1.0000    0.1881         0
        1.0000    0.1600         0
        1.0000    0.1263         0
        1.0000    0.0811         0
        1.0000    0.0026         0];
    
     RGB = xrain(1:floor(256/N_colours):256,:);
     RGB = RGB(1:N_colours,:);
   
end


