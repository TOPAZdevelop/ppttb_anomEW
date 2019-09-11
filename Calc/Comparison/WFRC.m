(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
Get[InitPath<>"EW_Renormalization.dat"];



convertNot={ Cpq1 -> -Cpq3, Cpq3 -> C33phiq3,   Cpu ->   C33phiu, A0[arg___] -> SI[1, arg],  B0[arg___] -> SI[2, arg], mt -> MT, Derivative[1, 0, 0][B0][args__] -> DB0[args], sw -> SW};


(*  Comparison of self energies  *)

SIGMASMEFTX0L=(alpha*(Cpq1 - Cpq3)*voL^2*((2*mt^2 - (Cpq1 - Cpq3)*(mt^2 + psq)*voL^2)*A0[mt^2] + (-2*mt^2 + (Cpq1 - Cpq3)*(mt^2 - psq)*voL^2)*A0[MZ^2] +   (2*mt^2*(-mt^2 + MZ^2 + psq) + (Cpq1 - Cpq3)*(mt^4 + psq*(-MZ^2 + psq) - mt^2*(MZ^2 + 2*psq))*voL^2)*B0[psq, mt^2, MZ^2]))/(32*MW^2*Pi*psq*sw^2)

SIGMASMEFTX0R=(alpha*Cpu*voL^2*(-((2*mt^2 + Cpu*(mt^2 + psq)*voL^2)*A0[mt^2]) + (2*mt^2 + Cpu*(mt^2 - psq)*voL^2)*A0[MZ^2] +   (2*mt^2*(mt^2 - MZ^2 - psq) + Cpu*(mt^4 + psq*(-MZ^2 + psq) - mt^2*(MZ^2 + 2*psq))*voL^2)*B0[psq, mt^2, MZ^2]))/(32*MW^2*Pi*psq*sw^2)

SIGMASMEFTX0S=(alpha*voL^2*((Cpq1 - Cpq3 - Cpu + 2*(Cpq1 - Cpq3)*Cpu*voL^2)*A0[mt^2] + (Cpq1 - Cpq3 - Cpu)*A0[MZ^2] + ((Cpq1 - Cpq3 - Cpu)*(mt^2 + MZ^2 - psq) + 2*(Cpq1 - Cpq3)*Cpu*MZ^2*voL^2)*    B0[psq, mt^2, MZ^2]))/(32*MW^2*Pi*sw^2)




TiSEChi0Left = SIGMASMEFTX0L //. convertNot // Expand;
TiSEChi0Right = SIGMASMEFTX0R //. convertNot // Expand;
TiSEChi0Scalar = SIGMASMEFTX0S //. convertNot // Expand;


MaSEChi0Left =sigmafL[psq,MT^2] //. EWSelfEnergies  //.  {ID[Chi] -> 1,  ID[__] -> 0} // Expand;
MaSEChi0LeftSM =MaSEChi0Left  //. voL->0  ;
MaSEChi0LeftBSM = MaSEChi0Left  - MaSEChi0LeftSM // Expand ;
MaSEChi0LeftBSM-TiSEChi0Left // Expand


MaSEChi0Right =sigmafR[psq,MT^2] //. EWSelfEnergies  //.  {ID[Chi] -> 1,  ID[__] -> 0} // Expand;
MaSEChi0RightSM =MaSEChi0Right  //. voL->0  ;
MaSEChi0RightBSM = MaSEChi0Right - MaSEChi0RightSM // Expand ;
MaSEChi0RightBSM-TiSEChi0Right // Expand



MaSEChi0Scalar =sigmafS[psq,MT^2] //. EWSelfEnergies  //.  {ID[Chi] -> 1,  ID[__] -> 0} // Expand;
MaSEChi0ScalarSM =MaSEChi0Scalar  //. voL->0  ;
MaSEChi0ScalarBSM = MaSEChi0Scalar - MaSEChi0ScalarSM // Expand ;
MaSEChi0ScalarBSM-TiSEChi0Scalar  // Expand


(*  Comparison of counter terms   *)

dZSMEFTX0L = ((alpha*voL^2*(-2*Cpu*A0[mt^2] + 2*Cpu*A0[MZ^2] - 2*Cpu*MZ^2*B0[mt^2, mt^2, MZ^2] -
    4*Cpq1*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] +
    4*Cpq3*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] +
    4*Cpu*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2]))/(32*MW^2*Pi*sw^2) +
 (alpha*voL^4*((Cpq1 - Cpq3)^2*A0[mt^2] - Cpu^2*A0[mt^2] + (Cpq1 - Cpq3)^2*A0[MZ^2] +
    Cpu^2*A0[MZ^2] + (Cpq1 - Cpq3)^2*MZ^2*B0[mt^2, mt^2, MZ^2] -
    Cpu^2*MZ^2*B0[mt^2, mt^2, MZ^2] + 2*Cpq1^2*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2,
      mt^2, MZ^2] - 4*Cpq1*Cpq3*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] +
    2*Cpq3^2*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] -
    4*Cpq1*Cpu*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] +
    4*Cpq3*Cpu*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2] +
    2*Cpu^2*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2]))/(32*MW^2*Pi*sw^2) );

dZSMEFTX0R = ((alpha*voL^2*(2*Cpq1*A0[mt^2] - 2*Cpq3*A0[mt^2] - 2*Cpq1*A0[MZ^2] + 2*Cpq3*A0[MZ^2] +
    2*Cpq1*MZ^2*B0[mt^2, mt^2, MZ^2] - 2*Cpq3*MZ^2*B0[mt^2, mt^2, MZ^2] +
    4*(-Cpq1 + Cpq3 + Cpu)*mt^2*MZ^2*Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2]))/
  (32*MW^2*Pi*sw^2) + (alpha*voL^4*((Cpq1 - Cpq3 + Cpu)*(-Cpq1 + Cpq3 + Cpu)*A0[mt^2] +
    ((Cpq1 - Cpq3)^2 + Cpu^2)*A0[MZ^2] + (Cpq1 - Cpq3 + Cpu)*(-Cpq1 + Cpq3 + Cpu)*MZ^2*
     B0[mt^2, mt^2, MZ^2] + 2*(-Cpq1 + Cpq3 + Cpu)^2*mt^2*MZ^2*
     Derivative[1, 0, 0][B0][mt^2, mt^2, MZ^2]))/(32*MW^2*Pi*sw^2) );




dZSMEFTPpmL=-(alpha*Cpq3*mt^2*voL^2*(B0[mt^2, 0, MW^2] + (mt^2 - MW^2)*Derivative[1, 0, 0][B0][mt^2, 0,MW^2]))/(8*MW^2*Pi*sw^2) + (alpha*Cpq3*voL^4*(-(Cpq3*(-2*A0[MW^2] + MW^2*B0[0, 0, MW^2] + (mt^2 - MW^2)*B0[mt^2, 0, MW^2])) - Cpq3*mt^2*(B0[mt^2, 0, MW^2] +   (mt^2 - MW^2)*Derivative[1, 0, 0][B0][mt^2, 0, MW^2])))/(16*MW^2*Pi*sw^2);


dZSMEFTPpmR=-(alpha*Cpq3*mt^2*voL^2*(B0[mt^2, 0, MW^2] + (mt^2 - MW^2)*Derivative[1, 0, 0][B0][mt^2, 0,MW^2]))/(8*MW^2*Pi*sw^2) - (alpha*Cpq3^2*mt^2*voL^4*(B0[mt^2, 0, MW^2] + (mt^2 - MW^2)*Derivative[1, 0, 0][B0][mt^2,0, MW^2]))/(16*MW^2*Pi*sw^2) ;





TiL = dZSMEFTX0L //.convertNot  // Expand;

MaL = delZiifL[MT^2] //. EWCounterterms //. {ID[Chi] -> 1,  ID[__] -> 0}  // Expand;
MaLSM = MaL //. voL -> 0 ;
MaLBSM = MaL - MaLSM // Expand;

MaLBSM-TiL // Expand // Factor


TiR = dZSMEFTX0R //.convertNot  // Expand;

MaR = delZiifR[MT^2] //. EWCounterterms //. {ID[Chi] -> 1,  ID[__] -> 0}  // Expand;
MaRSM = MaR //. voL -> 0 ;
MaRBSM = MaR - MaRSM // Expand;

MaRBSM-TiR // Expand // Factor


TiL = dZSMEFTPpmL //.convertNot //.  { SI[2,0,0,Msq_]-> 1/Msq SI[1,Msq] } // Expand;

MaL = delZiifL[MT^2] //. EWCounterterms //. {ID[phi] -> 1,  ID[__] -> 0} //. SU2Flip[MT^2]-> MB^2 //. MB->0 //.  { SI[1,0]->0 }  //Expand;
MaLSM = MaL //. voL -> 0 ;
MaLBSM = MaL - MaLSM // Expand;

MaLBSM- TiL // Expand // Factor


TiR = dZSMEFTPpmR //.convertNot  // Expand;

MaR = delZiifR[MT^2] //. EWCounterterms //. {ID[phi] -> 1,  ID[__] -> 0} //. SU2Flip[MT^2]-> MB^2 //. MB->0  //.  { SI[1,0]->0 }   //Expand;
MaRSM = MaR //. voL -> 0 ;
MaRBSM = MaR - MaRSM // Expand;

MaRBSM-TiR // Expand // Factor



