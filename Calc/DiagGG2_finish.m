(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[InitPath<>"EW_Renormalization.dat"];
Get[ ProjectPath<>"DiagGG2_output.dat" ];

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

convertToMCFM = {shat->s , SW->Sqrt[SW2], NCol->3,  voL^2-> vol2, voL^4-> vol4, vol2^2->vol4, SI[1,args___]->xI1[args,musq,ep],SI[2,args___]->xI2[args,musq,ep],SI[3,args___]->xI3[args,musq,ep],SI[4,args___]->xI4[args,musq,ep]};

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
         dummy = SIOrder[ dummy ];
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
inputTree = Expand[M[11,1]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmpS = inputTree //. {  dZfL->0, dZfR->0  };
TreeAmpCTS=inputTree-TreeAmpS  // FullSimplify


 (* D-dimensional tree and WFR CT diagrams   T-CHANNEL CT ONLY*)
inputTree = Expand[M[12,1]] //. insertMasses //. insertMT //. PropDenom[x_]->1/x   //.{ PreFac-> gs^4*ICol^2 } //. insertMT   // FullSimplify;
TreeAmpT = inputTree //. {  dZfL->0, dZfR->0  };
TreeAmpCTT=inputTree-TreeAmpT  // FullSimplify




TwoTimesRe = 2;






(* Phipm box diagram t-channel *)
(* Note that TheRedAmpList[[7]] \[Equal] TheRedAmpList[[6]] with z\[Rule]-z  , i.e. t-channel vs u-channel *)

PreFacDiagPhipmBox=EL^2 GS^4/(256 MW^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm

DiagPhipmBox= TwoTimesRe * TheRedAmpList[[7]]/PreFacDiagPhipmBox  //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2 , Sqrt2->Sqrt[2] } //. insertMT   // Collect[#,SI[___],Simplify]&;
% //. voL-> 0 
UVCheck=DiagPhipmBox  //. InsertUVDiv //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal // Factor
 % //. voL-> 0  // Factor





(* remember to replace arguments of xI function  0-->0d0  *)
bx4 = Expand[DiagPhipmBox]  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___]},FullSimplify]&  // FortranForm











(* s-channel Phipm vertex correction *)
DiagPhipmSVert = TheRedAmpList[[1]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2 , Sqrt2->Sqrt[2] } //. insertMT // Collect[#,SI[___],FullSimplify]&
UVCheck=DiagPhipmSVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Simplify



(* Phipm WFRC correction *)
CTDiagPhipmSVert = TreeAmpCTS  //. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[phi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi),SU2Flip[MT^2]-> MB^2 }  // Limit[#,MB->0 ]& // Limit[#,SI[1,0]->0 ]& // Expand;
CTDiagPhipmSVert = CTDiagPhipmSVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Collect[#,SI[___],FullSimplify]&


PreFacDiagPhipmSVert = EL^2 GS^4 /(256 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm

DiagPhipmSVertRENORM = TwoTimesRe*( DiagPhipmSVert + CTDiagPhipmSVert )/PreFacDiagPhipmSVert  // Collect[#,SI[___],FullSimplify]&
UVcheck =  DiagPhipmSVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL-> 0






(* convert to MCFM *)
vrts4 = Expand[DiagPhipmSVertRENORM]  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm











(* t-channel Phipm vertex corrections *)
DiagPhipmVert = (TheRedAmpList[[3]]+TheRedAmpList[[4]]) //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2 , Sqrt2->Sqrt[2] } //. insertMT   // Collect[#,SI[___],Simplify]&  
UVCheck=DiagPhipmVert  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* Phipm WFRC correction *)

CTDiagPhipmVert = TreeAmpCTT//. { dZfL->  delZiifL[MT^2] , dZfR->  delZiifR[MT^2]   }   //. EWCounterterms //. {ID[phi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi),SU2Flip[MT^2]-> MB^2 }  // Limit[#,MB->0 ]& // Limit[#,SI[1,0]->0 ]&  // Expand;
CTDiagPhipmVert = CTDiagPhipmVert  //. InsertFinitePartSquared //. DSTm4->0 //. insertMT  // Simplify


PreFacDiagPhipmVert = EL^2 GS^4 /(256 Pi^2 MW^2 SW^2 );
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm

DiagPhipmVertRENORM = TwoTimesRe*( DiagPhipmVert + CTDiagPhipmVert )/ PreFacDiagPhipmVert // Collect[#,SI[___],Simplify]&
UVcheck =  DiagPhipmVertRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL->0 


(* renormalized contribution *)


(* convert to MCFM *)
vrt4 =Expand[ DiagPhipmVertRENORM]  //. convertToMCFM // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm 












(* self energy corrections *)
DiagPhipmSelf = TheRedAmpList[[9]] //.{ PreFac-> EL^2*GS^4*SW^(-2)*MW^(-2)*Pi^2 , Sqrt2->Sqrt[2] } //. insertMT  // Collect[#,{SI[___]},Simplify]&  
UVCheck=DiagPhipmSelf  //. InsertUVDiv //. insertMT   // Series[#,{DSTm4,0,-1}]& // Normal // Factor


(* Phipm mass CT correction *)
Get[ ProjectPath<>"DiagGG4_output.dat" ];
MassCTPhipm = M[2,1]  // Simplify;
MassCTPhipm =MassCTPhipm   //. {dMf1[3,3]-> delM[MT^2] + MT EL^2/(128 MW^2 Pi^2 SW^2) voL^2 , dZfL1[__]-> delZiifL[MT^2],  dZfR1[__]-> delZiifR[MT^2] }   //. EWCounterterms //. {ID[phi]-> 1, ID[__]-> 0, alpha-> EL^2/(4 Pi), Sqrt[MT^2]-> MT, SU2Flip[MT^2]-> MB^2 } //. MB->0  // ReplaceAll[#,SI[1,0]->0 ]& // Expand; 
MassCTPhipm = MassCTPhipm  //. InsertFinitePartSquared //. DSTm4->0  //. insertMT  // ReplaceAll[#,SI[1,0]->0 ]&  // Simplify;
UVCheck = MassCTPhipm  //. InsertUVDiv //. insertMT // Series[#,{DSTm4,0,-1}]& // Normal // Factor


PreFacDiagPhipmSelf = EL^2 GS^4/(256 MW^2 Pi^2 SW^2);
% /. { EL^2-> alpha(4Pi), gs^4-> (alphas(4Pi))^2 } // FortranForm

DiagPhipmSelfRENORM = TwoTimesRe*( DiagPhipmSelf + MassCTPhipm )/PreFacDiagPhipmSelf  // Simplify
UVcheck =  DiagPhipmSelfRENORM  //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
% //. voL->0 


(* renormalized contribution *)


(* convert to MCFM *)
slf4 = Expand[DiagPhipmSelfRENORM]  //. convertToMCFM  // Collect[#,{xI1[___],xI2[___],xI3[___],xI4[___],DB0[___]},FullSimplify]&  // FortranForm  




























PreFacDiagPhipmBox
PreFacDiagPhipmVert
PreFacDiagPhipmSelf


UVCheckBox=DiagPhipmBox  //. InsertUVDiv //. insertMT  //. NCol->3 // Series[#,{DSTm4,0,-1}]& // Normal // Factor
UVcheckTCh =  DiagPhipmVertRENORM  //. InsertUVDiv //. NCol->3 //. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor
UVcheckSelf =  DiagPhipmSelfRENORM  //. InsertUVDiv  //. NCol->3//. DB0[MT^2,MT^2,__]->0  //. insertMT  // Series[#,{DSTm4,0,-1}]& // Normal  // Factor


UVCheckAll = DiagPhipmBox + DiagPhipmVertRENORM + DiagPhipmSelfRENORM //. InsertUVDiv //. DB0[MT^2,MT^2,__]->0  //. insertMT   //. NCol->3  // Series[#,{DSTm4,0,-1}]& // Normal   // Factor




