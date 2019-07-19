#define WorkPath "/users/pep/mschulze/projects/ppttb_anomEW/Calc/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MT^2"
#define setP4P4 "MT^2"
#define setDST   "DST"
#define setDSTm4 "DSTm4"
#include /users/pep/mschulze/lib/FeynArtsToForm/header2.frm



*  amps with closed quark triangle only

#include `WorkPath'TreeVirtGG3_input.frm



id MQDGen5=MB;
id MQUGen5=MT;
argument;
   id MQDGen5=MB;
   id MQUGen5=MT;
endargument;



id MU=0;
id MB=0;
argument;
   id MB=0;
   id MU=0;
endargument;

id IndexDelta(Col1?,Col2?) = SUND(Col1,Col2);
id DID(Col1?)=1;
id SumOver(?args) = 1;

        

#include `WorkPath'DiagGG3_calc.frm


id i_=cI;

id p1?.p2? = Dot(p1,p2);
argument;
  id p1?.p2? = Dot(p1,p2);
endargument;


Format mathematica;
Print;
Bracket PreFac,cI,EL,GS,ICol,SW,MW,Pi,PropDenom,TC;
.sort;


#write <`WorkPath'DiagGG3_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'DiagGG3_output.dat> "AmpList = \" `AmpList'\";\n";


#write <`WorkPath'DiagGG3_output.dat> "M[0,0] = (%E);\n ", [10,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[0,1] = (%E);\n ", [11,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[0,2] = (%E);\n ", [12,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[0,3] = (%E);\n ", [13,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[1,1] = (%E);\n ", [1,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[2,1] = (%E);\n ", [2,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[3,1] = (%E);\n ", [3,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[4,1] = (%E);\n ", [4,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[5,1] = (%E);\n ", [5,1,1];
#write <`WorkPath'DiagGG3_output.dat> "M[6,1] = (%E);\n ", [6,1,1];





.end
