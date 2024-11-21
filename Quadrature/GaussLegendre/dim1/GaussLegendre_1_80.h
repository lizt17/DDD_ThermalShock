 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
#ifndef model_GAUSSLEGENDRE_1_80_H_ 
#define model_GAUSSLEGENDRE_1_80_H_ 

namespace model
{

   template<>
   struct GaussLegendre<1,80>
   {
       static Eigen::Matrix<double,2,80> abcsissasAndWeights()
       {
           Eigen::Matrix<double,80,2> aw;
           aw<<2.230886741844684e-04, 5.724750015985250e-04, 
               1.175067800880780e-03, 1.331766794751813e-03, 
               2.886229517155559e-03, 2.090156562348372e-03, 
               5.354348750122084e-03, 2.845461225701647e-03, 
               8.575713630685211e-03, 3.596452384058952e-03, 
               1.254542970713513e-02, 4.341972634629166e-03, 
               1.725745547810009e-02, 5.080883020553228e-03, 
               2.270461682818226e-02, 5.812057060374064e-03, 
               2.887861934506303e-02, 6.534380796926635e-03, 
               3.577006141377709e-02, 7.246754020056055e-03, 
               4.336844871412099e-02, 7.948091791764018e-03, 
               5.166221028061480e-02, 8.637326028070578e-03, 
               6.063871616089334e-02, 9.313407104105382e-03, 
               7.028429666844427e-02, 9.975305439036517e-03, 
               8.058426320987228e-02, 1.062201305786168e-02, 
               9.152293065926792e-02, 1.125254512313777e-02, 
               1.030836412476968e-01, 1.186594143293233e-02, 
               1.152487899324784e-01, 1.246126788201571e-02, 
               1.279998512082008e-01, 1.303761788373669e-02, 
               1.413174073189499e-01, 1.359411375021836e-02, 
               1.551811778289859e-01, 1.412990802862446e-02, 
               1.695700505069401e-01, 1.464418479206482e-02, 
               1.844621134765640e-01, 1.513616087973267e-02, 
               1.998346885851242e-01, 1.560508709360704e-02, 
               2.156643659386449e-01, 1.605024933654932e-02, 
               2.319270395514340e-01, 1.647096969882261e-02, 
               2.485979440556077e-01, 1.686660749230568e-02, 
               2.656516924147276e-01, 1.723656022587717e-02, 
               2.830623145841220e-01, 1.758026452247637e-02, 
               3.008032970590153e-01, 1.789719697681206e-02, 
               3.188476232502568e-01, 1.818687495292667e-02, 
               3.371678146261495e-01, 1.844885731913724e-02, 
               3.557359725577439e-01, 1.868274511936522e-02, 
               3.745238208038639e-01, 1.888818218100038e-02, 
               3.935027485711672e-01, 1.906485565723904e-02, 
               4.126438540836767e-01, 1.921249650348001e-02, 
               4.319179885954280e-01, 1.933087988703796e-02, 
               4.512958007792077e-01, 1.941982552952564e-02, 
               4.707477814237896e-01, 1.947919798138469e-02, 
               4.902443083716030e-01, 1.950890682815334e-02, 
               5.097556916283971e-01, 1.950890682815360e-02, 
               5.292522185762103e-01, 1.947919798138491e-02, 
               5.487041992207924e-01, 1.941982552952597e-02, 
               5.680820114045719e-01, 1.933087988703815e-02, 
               5.873561459163232e-01, 1.921249650347928e-02, 
               6.064972514288330e-01, 1.906485565723864e-02, 
               6.254761791961361e-01, 1.888818218100087e-02, 
               6.442640274422557e-01, 1.868274511936536e-02, 
               6.628321853738508e-01, 1.844885731913834e-02, 
               6.811523767497438e-01, 1.818687495291763e-02, 
               6.991967029409849e-01, 1.789719697670745e-02, 
               7.169376854158782e-01, 1.758026452237431e-02, 
               7.343483075852724e-01, 1.723656022587725e-02, 
               7.514020559443924e-01, 1.686660749230588e-02, 
               7.680729604485659e-01, 1.647096969882287e-02, 
               7.843356340613546e-01, 1.605024933674343e-02, 
               8.001653114148759e-01, 1.560508709405667e-02, 
               8.155378865234361e-01, 1.513616087977949e-02, 
               8.304299494930596e-01, 1.464418479163358e-02, 
               8.448188221710137e-01, 1.412990802863829e-02, 
               8.586825926810503e-01, 1.359411375024276e-02, 
               8.720001487917990e-01, 1.303761788378308e-02, 
               8.847512100675208e-01, 1.246126788205941e-02, 
               8.969163587523024e-01, 1.186594143296471e-02, 
               9.084770693407318e-01, 1.125254512316611e-02, 
               9.194157367901279e-01, 1.062201305789053e-02, 
               9.297157033315558e-01, 9.975305439071698e-03, 
               9.393612838391073e-01, 9.313407104150458e-03, 
               9.483377897193856e-01, 8.637326028133795e-03, 
               9.566315512858792e-01, 7.948091791861595e-03, 
               9.642299385862230e-01, 7.246754020254405e-03, 
               9.711213806549361e-01, 6.534380796201349e-03, 
               9.772953831718172e-01, 5.812057060399247e-03, 
               9.827425445219000e-01, 5.080883020551354e-03, 
               9.874545702928641e-01, 4.341972634629781e-03, 
               9.914242863693145e-01, 3.596452384056482e-03, 
               9.946456512498782e-01, 2.845461225702831e-03, 
               9.971137704828448e-01, 2.090156562346326e-03, 
               9.988249321991191e-01, 1.331766794757913e-03, 
               9.997769113258154e-01, 5.724750015945393e-04; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 
