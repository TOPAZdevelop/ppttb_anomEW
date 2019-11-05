(* ::Package:: *)

InitPath = "/users/pep/mschulze/lib/MathLib/";
FeynArtsPath = "/users/pep/mschulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/users/pep/mschulze/lib/FeynArtsToForm/";
ProjectPath = "/users/pep/mschulze/projects/ppttb_anomEW/Calc/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[ ProjectPath<>"DiagGG3_output.dat" ];




Clear[TCReplList,SIList,TheRedAmpList];
TIReduction = Get[ InitPath<>"TIReduction.dat" ];

insertMasses = { p1.p1 -> la2,
                 p2.p2 -> la2,
                 p3.p3 -> MT^2,
                 p4.p4 -> MT^2,
  p1.p3 -> shat/4 (1 - beta z ),
  p1.p4 -> shat/4 (1 + beta z ),
  p2.p3 -> shat/4 (1 + beta z ),
  p2.p4 -> shat/4 (1 - beta z ),
   p3.p4 -> p1.p2 - MT^2,
p1.p2-> shat/2
};

convertToMCFM = {shat->s , SW->Sqrt[SW2], CW->Sqrt[CW2], NCol->3,  voL^2-> vol2, voL^4-> vol4,  vol2^2->vol4, SI[1,args___]->xI1[args,musq,ep],SI[2,args___]->xI2[args,musq,ep],SI[3,args___]->xI3[args,musq,ep],SI[4,args___]->xI4[args,musq,ep]};



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
   ];

TheRedAmpList[[1]] = TheRedAmpList[[1]]/(16 Pi^4); 
TheRedAmpList[[2]] = TheRedAmpList[[2]]/(16 Pi^4); 
TheRedAmpList[[3]] = TheRedAmpList[[3]]/(16 Pi^4);
TheRedAmpList[[4]] = TheRedAmpList[[4]]/(16 Pi^4); 
TheRedAmpList[[5]] = TheRedAmpList[[5]]/(16 Pi^4); 
TheRedAmpList[[6]] = TheRedAmpList[[6]]/(16 Pi^4);
TheRedAmpList[[7]] = TheRedAmpList[[7]]/(16 Pi^4); 
TheRedAmpList[[8]] = TheRedAmpList[[8]]/(16 Pi^4); 

SIList // TableForm 




TwoTimesRe = 2;




DiagZChiTriangle = TwoTimesRe * (TheRedAmpList[[1]]+TheRedAmpList[[2]]+TheRedAmpList[[3]]+TheRedAmpList[[4]])  //.{ la2->0 , MW->CW MZ, Sqrt2->Sqrt[2], NCol->3, gat-> 1/(4 SW CW)} //. SI[2,0,MT^2,MT^2]-> SI[1,MT^2]/MT^2-1 //. {SI[arg___]->SI[arg], MT -> Sqrt[shat/4 (1-beta^2)] } // Collect[#,SI[___],FullSimplify]&


UVcheck =  DiagZChiTriangle  //. InsertUVDiv  //. {SI[arg___]->SI[arg], MT -> Sqrt[shat/4 (1-beta^2)] } //. Sqrt2->Sqrt[2]  // Series[#,{DSTm4,0,-1}]& // Normal // FullSimplify


trizx = Expand[DiagZChiTriangle] //. { EL^2 -> alpha (4 Pi), gs^4 -> (alphas (4 Pi))^2 }   //.    convertToMCFM //   Collect[#, {xI1[___], xI2[___], xI3[___], xI4[___], DB0[___]}, FullSimplify] &  // FortranForm  
Tilltrizx = (-16*alpha*alphas^2*(-1 + beta^2)^2*Pi*s^2*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 +   (2*Cpq3 + Cpu)^2*vol4)*xI3[0, 0, s, MT^2, MT^2, MT^2, musq, ep])/ (CW^2*MZ^2*SW2*(-1 + beta^2*z^2)) ;
Tilltrizx // Factor//FullSimplify

dum1=DiagZChiTriangle //. { EL^2 -> alpha (4 Pi), GS^4 -> (alphas (4 Pi))^2, gat-> 1/(4 SW CW) } //. shat->s // FullSimplify

dum2=Tilltrizx   //. {Cpq1->-Cpq3,  Cpq3->C33phiq3, Cpu->C33phiu , xI3[args__,musq,ep]-> SI[3,args],vol2->voL^2,vol4->voL^4, SW2->SW^2} // FullSimplify

dum1-dum2 // Expand 





(* for the Higgs exchange we remove the C33phibox operator and rewrite everything in terms of kappa and kappa_tilde *)
eq1 = kap == 1 - 2 SW MW/EL/MT 1/Sqrt[2] * voL^2 * (C33uphi+C33uphiS)/2;
eq2 = kapT ==  - 2 SW MW/EL/MT 1/Sqrt[2] * voL^2 * (C33uphi-C33uphiS)/(2*I);
InsertKappa = Solve[{eq1,eq2},{C33uphi,C33uphiS}][[1]] // FullSimplify;
InsertKappa = Append[InsertKappa, C33phibox->0 ]


DiagHiggsTriangle = TwoTimesRe * (TheRedAmpList[[7]]+TheRedAmpList[[8]])  //.{ la2->0 , MW->CW MZ, Sqrt2->Sqrt[2], NCol->3, EL^2 -> alpha (4 Pi), GS^4 -> (alphas (4 Pi))^2, EL-> 2 Sqrt[alpha Pi]} //. InsertKappa //. {SI[arg___]->SI[arg], MT -> Sqrt[shat/4 (1-beta^2)] } // FullSimplify


trih = DiagHiggsTriangle //. convertToMCFM// FullSimplify 
% // FortranForm  



