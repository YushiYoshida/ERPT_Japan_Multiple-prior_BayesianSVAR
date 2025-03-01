# ERPT_Japan_Multiple-prior_BayesianSVAR
Matlab codes and the dataset for Yushi Yoshida and Weiyang Zhang, 2025(103312), Journal of International Money and Finance

The original main code is "MP_Ext_mainfile.m" from "Identification and Inference
under Narrative Restrictions" by Raffaella Giacomini, Toru Kitagawa, and
Matthew Read, 2021, arXiv:2102.06456
 
Yushi Yoshida modified the code for the
research paper, "Can exchange rate pass-throughs be perverse? 
A robust multiple-prior Bayesian SVAR approach," Yushi Yoshida and
Weiyang Zhai, available online Feb 2025, Journal of International Money and
Finance, 103312

The main contributions in these modified codes are the following.
(1) The original GKR code did not have zero restrictions, so we implanted the zero restriction part of the GK code. Of course, we had to make some adjustments.
(2) We added the ERPT range restrictions, allowing ERPT estimates between -1 and +1.
(3) Comparing graphically two alternative scenarios for narrative restrictions. (The code is original by Yoshida)

The main file is 'YZ_ERPT_mainfile.m'
After choosing the narrative restriction out of 7 alternatives by nr_yz = '', this code runs for more than a day.

To skip the long-hours estimation process, we provide the estimated results in the resultFiles folder.
With JPN_ERPT_YZ_MakingGraphOnly.m, you can enjoy the graph of two alternative narrative restrictions with one overlying another. 


The modified codes are indicated by the signature, @YZ2025JIMF
The following 8 codes in auxFunctions are also modified
(i)approximationBoundsYZ.m 
(ii)mainfileYZ.m 
(iii)intialComp_YZ.m
(iv)credibleRegionERPT.m 
(v)HighestPoseteriorDensityERPT.m
(vi)genVMA_YZ.m 
(vii)drawQ0eq.m 
(viii)drawQeq.m

Other 6 codes in auxFunctions are not modified.
