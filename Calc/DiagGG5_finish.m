(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[InitPath<>"EW_Renormalization.dat"];
Get[ ProjectPath<>"DiagGG5_output.dat" ];

Clear[TCReplList,SIList];
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

convertToMCFM = {shat->s , SW->Sqrt[SW2], NCol->3,  voL^2-> vol2, voL^4-> vol4,  vol2^2->vol4, SI[1,args___]->xI1[args,musq,ep],SI[2,args___]->xI2[args,musq,ep],SI[3,args___]->xI3[args,musq,ep],SI[4,args___]->xI4[args,musq,ep]};

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
         TheRedAmpList[[i]] = TheRedAmpList[[i]]/(16 Pi^4); 
   ];

SIList // TableForm 





(* D-dimensional tree and WFR CT diagrams *)
inputTree = Expand[M[0,0]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmp = inputTree //. {  dZfL->0, dZfR->0  }
TreeAmpCT=inputTree-TreeAmp  // FullSimplify


 (* D-dimensional tree and WFR CT diagrams   S-CHANNEL CT ONLY*)
inputTree = Expand[M[13,1]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmpS = inputTree //. {  dZfL->0, dZfR->0  };
TreeAmpCTS=inputTree-TreeAmpS  // FullSimplify


 (* D-dimensional tree and WFR CT diagrams   T-CHANNEL CT ONLY*)
inputTree = Expand[M[14,1]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmpT = inputTree //. {  dZfL->0, dZfR->0  };
TreeAmpCTT=inputTree-TreeAmpT  // FullSimplify




TwoTimesRe = 2;



(* Higgs box diagram  t-channel *)
(* Note that TheRedAmpList[[9]] \[Equal] TheRedAmpList[[8]] with z\[Rule]-z  , i.e. t-channel vs u-channel *)

PreFacDiagHiggsBox=EL^2 GS^4/(256 MW^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm

DiagHiggsBox= TwoTimesRe * TheRedAmpList[[9]]/PreFacDiagHiggsBox  //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. InsertKappa //. insertMT //. Sqrt2->Sqrt[2]  // Collect[#,{voL^2,voL^4,EL,SI[___]},FullSimplify]&
% //. voL-> 0 ;
UVCheck=DiagHiggsBox  //. InsertUVDiv //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal // Factor
 % //. voL-> 0  // Factor





(* remember to replace arguments of xI function  0-->0d0  *)

bx5 = Expand[DiagHiggsBox] //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___]},FullSimplify]&  // FortranForm


(* s-channel Higgs vertex correction *)
DiagHiggsSVert = TheRedAmpList[[1]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT //. Sqrt2->Sqrt[2] //Simplify
UVCheck=DiagHiggsSVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor;


(* Higgs WFRC correction *)

CTDiagHiggsSVert = TreeAmpCTS  //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Hig]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)} // Expand;
CTDiagHiggsSVert = CTDiagHiggsSVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  //. {Sqrt[EL^2]->EL, 1/Sqrt[EL^2]->1/EL}  // Simplify
UVCheck=CTDiagHiggsSVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor;


(* renormalized contribution *)
PreFacDiagHiggsSVert = EL^2 GS^4 MZ^2 /(64 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


DiagHiggsSVertRENORM = TwoTimesRe*( DiagHiggsSVert + CTDiagHiggsSVert )/PreFacDiagHiggsSVert //. InsertKappa  // Collect[#,SI[___],FullSimplify]&
UVcheck =  DiagHiggsSVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor


(* convert to MCFM *)
vrts5 =Expand[ DiagHiggsSVertRENORM ] //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm


(* t-channel Higgs vertex corrections *)
DiagHiggsVert = (TheRedAmpList[[3]]+TheRedAmpList[[4]]) //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT //. Sqrt2->Sqrt[2] // Collect[#,SI[___],Simplify]&  
UVCheck=DiagHiggsVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* Higgs WFRC correction *)
CTDiagHiggsVert = TreeAmpCTT//. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Hig]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)}  // Expand;
CTDiagHiggsVert = CTDiagHiggsVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  //. {Sqrt[EL^2]->EL, 1/Sqrt[EL^2]->1/EL} // Simplify


(* renormalized contribution *)
PreFacDiagHiggsVert = EL^2 GS^4 /(256 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


DiagHiggsVertRENORM = TwoTimesRe*( DiagHiggsVert + CTDiagHiggsVert )/PreFacDiagHiggsVert  // Collect[#,SI[___],Simplify]&
UVcheck =  DiagHiggsVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor


(* convert to MCFM *)
vrt3 =Expand[ DiagHiggsVertRENORM ] //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm 










(* self energy corrections *)
DiagHiggsSelf = TheRedAmpList[[12]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT //. {Sqrt2->Sqrt[2], 1/Sqrt2->1/Sqrt[2]} // Collect[#,{SI[___]},Simplify]&  
UVCheck=DiagHiggsSelf  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor;


(* Higgs mass CT correction *)
Get[ ProjectPath<>"DiagGG4_output.dat" ];
MassCTHiggs = M[2,1]  // Simplify;
MassCTHiggs =MassCTHiggs   //. {dMf1[3,3]-> delM[MT^2], dZfL1[__]-> delZiifL[MT^2],  dZfR1[__]-> delZiifR[MT^2] }   //. EWCounterterms //. {ID[Hig]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi), Sqrt[MT^2]-> MT}  // Expand; 
MassCTHiggs = MassCTHiggs  //. InsertFinitePartSquared //. DSTm4->0  //. insertMT  //. {Sqrt[EL^2]->EL, 1/Sqrt[EL^2]->1/EL}  // Simplify;
UVCheck = MassCTHiggs  //. InsertUVDiv //. insertMT // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* renormalized contribution *)
PreFacDiagHiggsSelf = EL^2 GS^4/(256 MW^2 Pi^2 SW^2);

DiagHiggsSelfRENORM = TwoTimesRe*( DiagHiggsSelf + MassCTHiggs )/PreFacDiagHiggsSelf  // Simplify
UVcheck =  DiagHiggsSelfRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor


PreFacDiagHiggsSelf = EL^2 GS^4/(256 MW^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


(* convert to MCFM *)
slf3 = Expand[DiagHiggsSelfRENORM ] //. convertToMCFM  // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm  



























PreFacDiagHiggsVert
PreFacDiagHiggsSelf
PreFacDiagHiggsBox


UVCheckBox=DiagHiggsBox  //. InsertUVDiv //. insertMT //. NCol->3 // Series[#,{DSTm4,0,-1}]& // Normal// Factor
UVcheckTCh =  DiagHiggsVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  //. NCol->3// Series[#,{DSTm4,0,-1}]& // Normal  // Factor
UVcheckSelf =  DiagHiggsSelfRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT //. NCol->3 // Series[#,{DSTm4,0,-1}]& // Normal// Factor

UVCheckAll = DiagHiggsBox + DiagHiggsVertRENORM + DiagHiggsSelfRENORM //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT //. NCol->3 // Series[#,{DSTm4,0,-1}]& // Normal// Factor






