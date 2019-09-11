(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[InitPath<>"EW_Renormalization.dat"];
Get[ ProjectPath<>"DiagQQB1_output.dat" ];


Clear[TCReplList,SIList,TheRedAmpList];
TIReduction = Get[ InitPath<>"TIReduction.dat" ];


insertMasses = { 
  p1.p1 -> 0,
  p2.p2 -> 0,
  p3.p3 -> MT^2,
  p4.p4 -> MT^2,
  p1.p3 -> shat/4 (1 - beta z ),
  p1.p4 -> shat/4 (1 + beta z ),
  p2.p3 -> shat/4 (1 + beta z ),
  p2.p4 -> shat/4 (1 - beta z ),
  p3.p4 -> p1.p2 - MT^2,
  p1.p2-> shat/2
};

insertMT={ SI[args___]->SI[args], DB0[args___]->DB0[args], MT -> Sqrt[shat/4 (1-beta^2)]};

convertToMCFM = {shat->s , SW->Sqrt[SW2], voL^2-> vol2, voL^4-> vol4, voL^-2-> 1/vol2, 1/voL^-4-> 1/vol4, SI[1,args___]->xI1[args,musq,ep],SI[2,args___]->xI2[args,musq,ep],SI[3,args___]->xI3[args,musq,ep]};

InsertFinitePartSquared = {r_.* DSTm4* SI[1,m0sq_]-> r* (-2*m0sq), r_.* DSTm4* SI[2,p1_,m0sq_,m1sq_]-> r* (-2)}; 


(* for the Higgs exchange we remove the C33phibox operator and rewrite everything in terms of kappa and kappa_tilde *)
eq1 = kap == 1 - 2 SW MW/EL/MT 1/Sqrt[2] * voL^2 * (C33uphi+C33uphiS)/2;
eq2 = kapT ==  - 2 SW MW/EL/MT 1/Sqrt[2] * voL^2 * (C33uphi-C33uphiS)/(2*I);
InsertKappa = Solve[{eq1,eq2},{C33uphi,C33uphiS}][[1]] // FullSimplify;
InsertKappa = Append[InsertKappa, C33phibox->0 ]


TheAmpList    = StringReplace[AmpList,"["->"M["] // ToExpression;
TheRedAmpList = StringReplace[AmpList,"["->"redM["] // ToExpression;

(* generate list of all tensor integral coefficients *)
   TCList={};
   For[i=1,i<=Length[TheAmpList],i++,
         TCList = Union[TCList, Cases[ TheAmpList[[i]] ,TC[___],100] ];
   ];

   (* simplify this list *)
   TCReplList={};
   For[i=1, i<=Length[TCList],i++,
         Print[i,"/",Length[TCList],": simplifying ",TCList[[i]]];
         
         TheCoeff = TCList[[i]] //.TIReduction //Expand;
         TheCoeff = TheCoeff    //.insertMasses;
         TheCoeff = TheCoeff    // Expand;
         TheCoeff = SimplifyTC[ TheCoeff ];
         TheCoeff = Collect[ TheCoeff, SI[___], Simplify];
         TCReplList = Append[TCReplList,TCList[[i]]-> TheCoeff ];
   ];
 
 

  SIList={};
  For[i=1,i<=Length[TheAmpList],i++,
         Print["reducing diagram ",i,"/",Length[TheAmpList]];
         dummy = TheAmpList[[i]] //. TCReplList; 
         Print["expanding in D-4"];
         dummy = Expand[ dummy ] //. InsertFinitePart   /. DSTm4->0;
         Print["converting to invariants"];
         dummy = SIInvariants[ dummy,insertMasses];         
         dummy = SIOrder[ dummy ]  //. insertMasses;
         dummy = dummy //. { SI[1,0]->0 };
         Print["simplifying amplitude"];
         dummy = dummy //. PropDenom[x_]->1/x    // SimplAmp;
         Evaluate[ TheRedAmpList[[i]] ] = dummy;


         (* generate list of integrals in terms of invariants *)
         SIList = Union[SIList, Cases[ TheRedAmpList[[i]] ,SI[___],100] ] //SIOrder //Simplify;
   ];
SIList // TableForm









(* D-dimensional tree and WFR CT diagrams *)
TreeAmpCT = Expand[M[1,1]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmp = TreeAmpCT //. {  dZfL->0, dZfR->0  }
TreeAmpCT=TreeAmpCT-TreeAmp  // FullSimplify

(*  note: ICol^2 is 2, see Bachelor Thesis Bergner *)




(* Chi0 boson exchange *)
DiagChi0=TheRedAmpList[[2]]/(16 Pi^4)  //.{ PreFac-> gs^4*ICol^2 } //. insertMT;
% //. voL->0   // FullSimplify


(* Chi0 WFRC correction *)
CTDiagChi0 = TreeAmpCT //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Chi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)}   // Expand;
CTDiagChi0 = CTDiagChi0  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // FullSimplify


(* renormalized contribution *)
PreFacRen = EL^2 gs^4 ICol^2/(64 Pi^2 SW^2 MW^2 beta^2);

DiagChi0RENORM = 2*( DiagChi0 + CTDiagChi0 )/PreFacRen   // Collect[#,{SI[___],DB0[___]},FullSimplify]& 
UVcheck =  DiagChi0RENORM  //. InsertUVDiv  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // FullSimplify


(* convert to MCFM *)
qa3 = Expand[DiagChi0RENORM]  //. convertToMCFM // Collect[#,{SI[___],DB0[___]},FullSimplify]&  // FortranForm











(* Phi+/- boson exchange *)
DiagPhiPM0=TheRedAmpList[[3]]/(16 Pi^4)  //.{ PreFac-> gs^4*ICol^2, 1/Sqrt2^2-> 1/2 } //. insertMT;
% //. voL->0 // FullSimplify


(* Phi+/- WFRC correction *)
CTDiagPhiPM = TreeAmpCT //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[phi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi), SU2Flip[MT^2]->0, SI[1,0]->0  }  // Expand;
CTDiagPhiPM = CTDiagPhiPM  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // FullSimplify


(* renormalized contribution *)
PreFacRen = EL^2 gs^4 ICol^2/(64 Pi^2 SW^2 MW^2 beta^2);

DiagPhi0RENORM = 2*( DiagPhiPM0 + CTDiagPhiPM )/PreFacRen // Collect[#,{SI[___],DB0[___]},FullSimplify]& 
UVcheck =  DiagPhi0RENORM  //. InsertUVDiv  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // FullSimplify


(* convert to MCFM *)
(* remember to replace 0 --> 0d0 in arguments of integrals  *)
qa4 = Expand[DiagPhi0RENORM] //. convertToMCFM  // Collect[#,{SI[___],DB0[___]},FullSimplify]& // FortranForm

















(* Higgs boson exchange *)
DiagHiggs=TheRedAmpList[[4]]/(16 Pi^4)  //.{ PreFac-> gs^4*ICol^2 } //. insertMT //. Sqrt2->Sqrt[2];
% //. voL->0   // FullSimplify


(* Chi0 WFRC correction *)
CTDiagHiggs = TreeAmpCT //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Hig]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)}  // Expand;
CTDiagHiggs = CTDiagHiggs  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  //. Sqrt[EL^2]->EL //. 1/Sqrt[EL^2]->1/EL  // FullSimplify


(* renormalized contribution *)
PreFacRen = EL^2 gs^4 ICol^2/(64 Pi^2 SW^2 MW^2 beta^2);

DiagHiggsRENORM = 2*( DiagHiggs + CTDiagHiggs )/PreFacRen   //. InsertKappa  // Collect[#,{voL^2,voL^4,EL,SI[___],DB0[___]},FullSimplify]& 
UVcheck =  DiagHiggsRENORM  //. InsertUVDiv  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // FullSimplify 


(* convert to MCFM *)
(* remember to replace 0 --> 0d0 in arguments of integrals  *)
qa5 = Expand[DiagHiggsRENORM]  //. insertMT //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm










