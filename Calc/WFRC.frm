#define WorkPath "/users/pep/mschulze/projects/ppttb_anomEW/Calc/"
#define setP1P1 "shat"
#define setP2P2 "mX^2"
#define setP3P3 "mX^2"
#define setP4P4 "mX^2"
#define setDST   "DST"
#define setDSTm4 "DSTm4"
#include /users/pep/mschulze/lib/FeynArtsToForm/header2.frm




#include `WorkPath'TreeWFRC_input.frm

  

id MB=0;
argument;
   id MB=0;
endargument;

id IndexDelta(Col1?,Col2?) = 1;
id SumOver(?args) = 1;
id DID(Col1?)=1;
        

#include `WorkPath'WFRC_calc.frm


id i_=cI;
id shat=p1.p1;


id p1?.p2? = Dot(p1,p2);
argument;
  id p1?.p2? = Dot(p1,p2);
endargument;



Format mathematica;
Print;
Bracket Ga,Chir;
.sort;



#write <`WorkPath'WFRC_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'WFRC_output.dat> "AmpList = \"  `AmpList'\";\n";

#write <`WorkPath'WFRC_output.dat> "M[1,1] = (%E);\n ", [1,1];
#write <`WorkPath'WFRC_output.dat> "M[2,1] = (%E);\n ", [2,1];
#write <`WorkPath'WFRC_output.dat> "M[3,1] = (%E);\n ", [3,1];





.end
