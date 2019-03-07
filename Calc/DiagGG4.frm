#define WorkPath "/users/pep/mschulze/projects/ppttb_anomEW/Calc/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MT^2"
#define setP4P4 "MT^2"
#define setDST   "DST"
#define setDSTm4 "DSTm4"
#include /users/pep/mschulze/lib/FeynArtsToForm/header2.frm


* Mass counter term insertions

#include `WorkPath'TreeVirtGG4_input.frm




id IndexDelta(Col1?,Col2?) = SUND(Col1,Col2);
id DID(Col1?)=1;
id SumOver(?args) = 1;

        

#include `WorkPath'DiagGG4_calc.frm


id i_=cI;

id EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2=PreFac;

id p1?.p2? = Dot(p1,p2);
argument;
  id p1?.p2? = Dot(p1,p2);
endargument;


Format mathematica;
Print;
Bracket PreFac,cI,EL,GS,ICol,SW,MW,Pi,PropDenom,TC;
.sort;


#write <`WorkPath'DiagGG4_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'DiagGG4_output.dat> "AmpList = \"  `AmpList'\";\n";

#write <`WorkPath'DiagGG4_output.dat> "M[1,1] = (%E);\n ", [1,1,1];
#write <`WorkPath'DiagGG4_output.dat> "M[2,1] = (%E);\n ", [2,1,1];


.end
