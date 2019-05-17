#define WorkPath "/users/pep/mschulze/projects/ppttb_anomEW/Calc/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MT^2"
#define setP4P4 "MT^2"
#define setDST   "DST"
#define setDSTm4 "DSTm4"
#include /users/pep/mschulze/lib/FeynArtsToForm/header2.frm




#include `WorkPath'TreeVirtQQB1_input.frm

Local test =LoopDenom(q1,MZ,p3 + q1,MT, - p4 + q1,MT)*q1.q1;  
  

id MU=0;
id MB=0;
argument;
   id MB=0;
   id MU=0;
endargument;

id IndexDelta(Col1?,Col2?) = 1;
id SumOver(?args) = 1;

id SUNT(GluInt5,Col2,Col1)*SUNT(GluInt5,Col3,Col4) = ICol;
         

#include `WorkPath'DiagQQB1_calc.frm


id i_=cI;
id GS^4*ICol^2 = PreFac;

id p1?.p2? = Dot(p1,p2);
argument;
  id p1?.p2? = Dot(p1,p2);
endargument;


Format mathematica;
Print;
Bracket PreFac,cI,EL,GS,ICol,SW,MW,Pi;
.sort;


#write <`WorkPath'DiagQQB1_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'DiagQQB1_output.dat> "AmpList = \"  `AmpList'\";\n";


#write <`WorkPath'DiagQQB1_output.dat> "M[1,1] = (%E);\n ", [1,1,1];
#write <`WorkPath'DiagQQB1_output.dat> "M[2,1] = (%E);\n ", [2,1,1];
#write <`WorkPath'DiagQQB1_output.dat> "M[3,1] = (%E);\n ", [3,1,1];
#write <`WorkPath'DiagQQB1_output.dat> "M[4,1] = (%E);\n ", [4,1,1];





.end
