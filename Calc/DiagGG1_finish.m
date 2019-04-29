(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[InitPath<>"EW_Renormalization.dat"];
Get[ ProjectPath<>"DiagGG1_output.dat" ];

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

convertToMCFM = {shat->s , SW->Sqrt[SW2], NCol->3,  voL^2-> vol2, voL^4-> vol4, SI[1,args___]->xI1[args,musq,ep],SI[2,args___]->xI2[args,musq,ep],SI[3,args___]->xI3[args,musq,ep],SI[4,args___]->xI4[args,musq,ep]};

InsertFinitePartSquared = {r_.* DSTm4* SI[1,m0sq_]-> r* (-2*m0sq), r_.* DSTm4* SI[2,p1_,m0sq_,m1sq_]-> r* (-2)}; 





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








(* Chi0 box diagram  t-channel *)
(* Note that TheRedAmpList[[9]] \[Equal] TheRedAmpList[[8]] with z\[Rule]-z  , i.e. t-channel vs u-channel *)

DiagChi0Box=TheRedAmpList[[9]]  //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT   // Collect[#,SI[___],Simplify]&;
% //. voL-> 0 
UVCheck=DiagChi0Box  //. InsertUVDiv //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal // Factor
 % //. voL-> 0  // Factor


PreFacDiagChi0Box=EL^2 GS^4/(128 MZ^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


(* remember to replace arguments of xI function  0-->0d0  *)

bx3 = DiagChi0Box/PreFacDiagChi0Box  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___]},FullSimplify]&  // FortranForm






(* s-channel Chi0 vertex correction *)
DiagChi0SVert = TheRedAmpList[[1]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT //Simplify
UVCheck=DiagChi0SVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(******  OUTDATED 

(* Chi0 WFRC correction *)
Get[ProjectPath<>"Selfenergies_output.dat"];

SigL = SigmatChiL //.TIReduction // SIInvariants[ #,{p1.p1->psq}]& // Collect[#,{voL},FullSimplify]&;
SigR = SigmatChiR //.TIReduction // SIInvariants[ #,{p1.p1->psq}]& // Collect[#,{voL},FullSimplify]&;
SigS = SigmatChiS //.TIReduction // SIInvariants[ #,{p1.p1->psq}]& // Collect[#,{voL},FullSimplify]&;

deltaZChi0L = WFRCPreFactor*( - SigL - MT^2 ( D[SigL,psq]+ D[SigR,psq] + 2 *D[SigS,psq]) )  //. psq-> MT^2   //. Derivative[0,1,0,0][SI][2,psq_,m1sq_,m2sq_]-> DB0[psq,m1sq,m2sq];
deltaZChi0R = WFRCPreFactor*( - SigR - MT^2 ( D[SigL,psq]+ D[SigR,psq] + 2 *D[SigS,psq]) )  //. psq-> MT^2   //. Derivative[0,1,0,0][SI][2,psq_,m1sq_,m2sq_]-> DB0[psq,m1sq,m2sq];

CTDiagChi0SVert = TreeAmpCTS //. { dZfL-> deltaZChi0L, dZfR->deltaZChi0R, alphaO4Pi-> alpha/(4 Pi) , alpha-> EL^2/(4 Pi)  }   // Expand;
CTDiagChi0SVert = CTDiagChi0SVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Simplify

*********)


(* Chi0 WFRC correction *)

CTDiagChi0SVert = TreeAmpCTS  //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Chi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)}  // Expand;
CTDiagChi0SVert = CTDiagChi0SVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Simplify


(* renormalized contribution *)
PreFacDiagChi0SVert = EL^2 GS^4 MZ^2 /(64 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


DiagChi0SVertRENORM = 2*( DiagChi0SVert + CTDiagChi0SVert )/PreFacDiagChi0SVert  // Collect[#,SI[___],FullSimplify]&
UVcheck =  DiagChi0SVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL-> 0


(* convert to MCFM *)
vrts3 = DiagChi0SVertRENORM  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm








(* t-channel Chi0 vertex corrections *)
DiagChi0Vert = (TheRedAmpList[[3]]+TheRedAmpList[[4]]) //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT   // Collect[#,SI[___],Simplify]&  
UVCheck=DiagChi0Vert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(******  OUTDATED 
CTDiagChi0Vert = TreeAmpCTT //. { dZfL-> deltaZChi0L, dZfR->deltaZChi0R, alphaO4Pi-> alpha/(4 Pi) , alpha-> EL^2/(4 Pi)  }   // Expand;
CTDiagChi0Vert = CTDiagChi0Vert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Simplify

******)


(* Chi WFRC correction *)

CTDiagChi0Vert = TreeAmpCTT//. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[Chi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi)}  // Expand;
CTDiagChi0Vert = CTDiagChi0Vert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Simplify


(* renormalized contribution *)
PreFacDiagChi0Vert = EL^2 GS^4 /(256 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


DiagChi0VertRENORM = 2*( DiagChi0Vert + CTDiagChi0Vert )/PreFacDiagChi0Vert  // Collect[#,SI[___],Simplify]&
UVcheck =  DiagChi0VertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL->0 


(* convert to MCFM *)
vrt3 = DiagChi0VertRENORM  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm 










(* self energy corrections *)
DiagChi0Self = TheRedAmpList[[12]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2  } //. insertMT  // Collect[#,{SI[___]},Simplify]&  
UVCheck=DiagChi0Self  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* Chi0 mass CT correction *)
Get[ ProjectPath<>"DiagGG4_output.dat" ];
MassCTChi0 = M[2,1]  // Simplify;
MassCTChi0 =MassCTChi0   //. {dMf1[3,3]-> delM[MT^2], dZfL1[__]-> delZiifL[MT^2],  dZfR1[__]-> delZiifR[MT^2] }   //. EWCounterterms //. {ID[Chi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi), Sqrt[MT^2]-> MT} // Expand; 
MassCTChi0 = MassCTChi0  //. InsertFinitePartSquared //. DSTm4->0  //. insertMT // FullSimplify;
UVCheck = MassCTChi0  //. InsertUVDiv //. insertMT // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* renormalized contribution *)
PreFacDiagChi0Self = EL^2 GS^4/(256 MW^2 Pi^2 SW^2);

DiagChi0SelfRENORM = 2*( DiagChi0Self + MassCTChi0 )/PreFacDiagChi0Self  // Simplify
UVcheck =  DiagChi0SelfRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL->0 


PreFacDiagChi0Self = EL^2 GS^4/(256 MW^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm


(* convert to MCFM *)
slf3 = DiagChi0SelfRENORM  //. convertToMCFM  // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm  
