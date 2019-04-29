(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath <> "StandardInit.m" ];
Get[ InitPath <> "SIOrder.m" ];
Get[ ProjectPath <> "WFRC_output.dat" ];

Clear[TCReplList, SIList, TheRedAmpList];
TIReduction = Get[ InitPath <> "TIReduction.dat" ];

TheAmpList    = StringReplace[AmpList, "[" -> "M["] // ToExpression;
TheRedAmpList = StringReplace[AmpList, "[" -> "redM["] // ToExpression;
  
  







WFRCPreFactor = -alphaO4Pi/(4 MW^2 SW^2);

WFRCChi0=M[2, 1]/(16 Pi^4) /. EL^2-> alphaO4Pi (4 Pi)^2 // FullSimplify;
WFRCPhiPM=M[3, 1]/(16 Pi^4) /. EL^2-> alphaO4Pi (4 Pi)^2 //. 1/Sqrt2^2->1/2 // FullSimplify;

WFRCChi0=WFRCChi0/WFRCPreFactor  // Collect[#,{Ga[___],Chir[__],voL},FullSimplify]&;
WFRCPhiPM=WFRCPhiPM/WFRCPreFactor  // Collect[#,{Ga[___],Chir[__],voL},FullSimplify]&;


SigmatChiL = Coefficient[WFRCChi0,Ga[1,p1,-1]]/cI // Collect[#,{voL},FullSimplify]&
SigmatChiR = Coefficient[WFRCChi0,Ga[1,p1,+1]]/cI // Collect[#,{voL},FullSimplify]&
SigmatChiS = Coefficient[WFRCChi0,Chir[1,-1]]/cI/MT // Collect[#,{voL},FullSimplify]&








SigmatPhiL = Coefficient[WFRCPhiPM,Ga[1,p1,-1]]/cI // Collect[#,{voL},FullSimplify]&
SigmatPhiR = Coefficient[WFRCPhiPM,Ga[1,p1,+1]]/cI // Collect[#,{voL},FullSimplify]&
SigmatPhiS = Coefficient[WFRCPhiPM,Chir[1,-1]]/cI/MT // Collect[#,{voL},FullSimplify]&


(* ::InheritFromParent:: *)
(**)


(* this output file is no longer in use because the WFRC have been incoporated into /mnt/pep/mschulze/lib/MathLib/EW_Renormalization.dat


Save["Selfenergies_output.dat",{WFRCPreFactor,SigmatChiL,SigmatChiR,SigmatChiS,SigmatPhiL,SigmatPhiR,SigmatPhiS}];

*)
