#define WorkPath "/users/pep/mschulze/projects/ppttb_anomEW/Calc/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MT^2"
#define setP4P4 "MT^2"
#define setDST   "DST"
#define setDSTm4 "DSTm4"
#include /users/pep/mschulze/lib/FeynArtsToForm/header2.frm




* charged Goldstone boson exchanges

#include `WorkPath'TreeVirtGG2_input.frm



* removing also diags 2,5,6,8 because they are the u-channel diagrams that can be obtained from the remaining t-channel diagrams


*id DID(1) =0;
id DID(2) =0;
*id DID(3) =0;
*id DID(4) =0;
id DID(5) =0;
id DID(6) =0;
*id DID(7) =0;
id DID(8) =0;
*id DID(9) =0;

  
  
id MU=0;
id MB=0;
argument;
   id MB=0;
   id MU=0;
endargument;

id IndexDelta(Col1?,Col2?) = SUND(Col1,Col2);

id SumOver(?args) = 1;

         

#include `WorkPath'DiagGG2_calc.frm





id i_=cI;
id EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2 = PreFac;

id p1?.p2? = Dot(p1,p2);
argument;
  id p1?.p2? = Dot(p1,p2);
endargument;


Format mathematica;
Print;
Bracket PreFac,cI,EL,GS,ICol,SW,MW,Pi,PropDenom,TC,SI,voL;
.sort;


#write <`WorkPath'DiagGG2_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'DiagGG2_output.dat> "AmpList = \"  `AmpList'\";\n";


#write <`WorkPath'DiagGG2_output.dat> "M[0,0] = (%E);\n ", [10,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[1,1] = (%E);\n ", [1,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[2,1] = (%E);\n ", [2,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[3,1] = (%E);\n ", [3,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[4,1] = (%E);\n ", [4,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[5,1] = (%E);\n ", [5,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[6,1] = (%E);\n ", [6,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[7,1] = (%E);\n ", [7,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[8,1] = (%E);\n ", [8,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[9,1] = (%E);\n ", [9,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[11,1] = (%E);\n ", [11,1,1];
#write <`WorkPath'DiagGG2_output.dat> "M[12,1] = (%E);\n ", [12,1,1];





.end
