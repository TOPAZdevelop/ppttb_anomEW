  << HypothesisTesting`
(*--------------------------------------------------------------------------*)

SetCouplings[Cphiq_,Cphiu_]:=Module[{},
 
  gv = gvSM + (vev/Lambda)^2 * ( Cphiq - 0.5*Cphiu );
  ga = gaSM + (vev/Lambda)^2 * ( Cphiq + 0.5*Cphiu );
  gw = gwSM + (vev/Lambda)^2 * ( 0.5*Cphiq );
   
  (* normalizing to SM *)
  gv = gv/gvSM;
  ga = ga/gaSM;
  gw = gw/gwSM;
  grest = 1;
  gasq = ga^2;
  gvsq = gv^2;
  gwsq = gw^2; 
 
];

(*--------------------------------------------------------------------------*)


GenerateHistogram[InputHisto_]:=Module[{},

  InputHisto={}; 
  TotalXSec = 0;  TotalEvts = 0;
  Do[
     NewBin = ( InputCoupl[[iBin]][[2]] * grest 
               +InputCoupl[[iBin]][[3]] * gwsq
               +InputCoupl[[iBin]][[4]] * gasq
               +InputCoupl[[iBin]][[5]] * gvsq
               +InputCoupl[[iBin]][[6]] * ga  
               +InputCoupl[[iBin]][[7]] * gv  
               +InputLO[[iBin]][[2]] );
               
      NewBin = NewBin *BinSize * Lumi;
     
     InputHisto = Append[InputHisto,{InputCoupl[[iBin]][[1]],NewBin}];

     TotalEvts = TotalEvts + NewBin;     
     TotalXSec = TotalXSec + NewBin/Lumi;
   
   ,{iBin,1,NumBins} ];
   
    InputHisto=Prepend[InputHisto,{Cphiq,Cphiu,TotalXSec}];
   
];

(*--------------------------------------------------------------------------*)


CalculateChisq[ChiSq_,degeof_]:=Module[{},

  TotalXSecSM = Histo[0.00,0.00][[1]][[3]];
  TotalXSecBSM= Histo[Cphiq,Cphiu][[1]][[3]];

  (*Print["rescaling cross sections by ",DeltaN,"% uncertainty"];  *)
  If[ TotalXSecSM>TotalXSecBSM, 
     MaxAllowedDelta = 0.5*(TotalXSecSM/TotalXSecBSM-1) * 100;
     AvgTotalEvents  = 0.5*(TotalXSecSM+TotalXSecBSM);
     If[ DeltaN<MaxAllowedDelta,
       Delta1 = -DeltaN/100  (*make null hypothesis smaller*);
       Delta2 = +DeltaN/100  (*make alt. hypothesis larger*);
     ,
       Delta1 = AvgTotalEvents/TotalXSecSM  - 1; (* rescale such that total cross section is equal to the average between the two hypothesis *)
       Delta2 = AvgTotalEvents/TotalXSecBSM - 1;
     ];
  ,
     MaxAllowedDelta = 0.5*(TotalXSecBSM/TotalXSecSM-1) * 100;
     AvgTotalEvents  = 0.5*(TotalXSecSM+TotalXSecBSM);     
     Delta1 = +DeltaN/100  (*make null hypothesis larger *);
     Delta2 = -DeltaN/100  (*make alt. hypothesis smaller*);     
     If[ TotalXSecSM*(1+Delta1)>TotalXSecBSM*(1+Delta2),
       Delta1 = AvgTotalEvents/TotalXSecSM  - 1; (* rescale such that total cross section is equal to the average between the two hypothesis *)
       Delta2 = AvgTotalEvents/TotalXSecBSM - 1;
     ];  
  ];
  

  Unct1 = 1+Delta1;
  Unct2 = 1+Delta2;

 

  ChiSq=0;
  Do[
     If[ Histo[0.00,0.00][[iBin]][[2]]*Unct1<30 || Histo[Cphiq,Cphiu][[iBin]][[2]]*Unct2<30, 
           (*Print["Less than 30 events beyond bin ",iBin-1]; *)
           degeof = iBin-1-1;
           Return[]; 
     ,
           degeof = iBin-1;
     ];
     
     ChisqBin = (Histo[0.00,0.00][[iBin]][[2]]*Unct1-Histo[Cphiq,Cphiu][[iBin]][[2]]*Unct2)^2/(Histo[0.00,0.00][[iBin]][[2]]*Unct1);
     ChiSq = ChiSq + ChisqBin;
          
  ,{iBin,2,Length[Histo[0.00,0.00]]} ];  
   
];
  
  
(*--------------------------------------------------------------------------*)


CalculateSigni[Chi2_,degeof_,Sig_]:=Module[{omalpha},

   If[ Chi2<degeof, Sig=0; Return[]; ];

   omalpha = 1 - OneSidedPValue /. ChiSquarePValue[Chi2,degeof];
   Sig= Sqrt[2] * InverseErf[omalpha];
];


(*--------------------------------------------------------------------------*)
(*-----------------------      MAIN PROGRAM      ---------------------------*)
(*--------------------------------------------------------------------------*)


  MyPath = "/mnt/pep/mschulze/pool/projects/TOPAZ/data/ttbarEW/";


  InputCoupl =  ReadList[MyPath <> "input_couplings.dat", {Number, Number, Number, Number, Number, Number, Number}];
  InputLO    =  ReadList[MyPath <> "input_LO.dat", {Number, Number}];
  BinSize = 50;
  Lumi = 300;
  PreFactor = 1/3; (* ttbar semi-hadr. decays *)
  Lumi = Lumi * PreFactor;
  DeltaN = 4;  (* system. uncertainty in percent *)
  
  NumBins=Length[InputCoupl];
  

  
  sw2 = 0.2228972;
  cw2 = 1.0-sw2;
  gvSM = (0.5-4*sw2/3)/2/Sqrt[sw2*cw2];
  gaSM = 0.5/2/Sqrt[sw2*cw2];
  gwSM = 0.5/Sqrt[2*sw2];
  vev = 246;
  Lambda = 1000;
  


  Cphiq=+00.00;
  Cphiu=+00.00;
  SetCouplings[Cphiq,Cphiu];
  GenerateHistogram[Histo[Cphiq,Cphiu]];
  
  
  (*
  Cphiq=-02.00; 
  Cphiu=+00.00;
  SetCouplings[Cphiq,Cphiu];
  GenerateHistogram[Histo[Cphiq,Cphiu]];
  
  Clear[ChiSq,dof,Signif];
  CalculateChisq[ChiSq,dof]; 
  CalculateSigni[ChiSq,dof,Signif]; 
  
  Print["ChiSq = ",ChiSq," with ",dof," dof"];
  Print["Significance = ",Signif];
  
  *)
  
 
GetSignificance[xCphiq_,xCphiu_]:=Module[{Si},

  If[ xCphiq==0.00 && xCphiu==0.00, Return[]; ];

  Cphiq=xCphiq;
  Cphiu=xCphiu;
  SetCouplings[Cphiq,Cphiu];
  GenerateHistogram[Histo[Cphiq,Cphiu]];
  
  Clear[ChiSq,dof];
  CalculateChisq[ChiSq,dof]; 
  CalculateSigni[ChiSq,dof,Si]; 

  Si
];

   

pp=(vev/Lambda)^2 //N;
Fig1=ContourPlot[GetSignificance[CPhiQ/pp,CPhiU/pp] , {CPhiU,-0.4,2.0}, {CPhiQ, -0.8,0.4},Contours->{1.0,2.0,3.0},PlotLegends->Automatic,FrameLabel->Automatic]
