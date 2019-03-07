
* eval propagator functions
id PropDenom(pDumX?,mX?) = PropDenom(pDumX.pDumX-mX^2);

argument PropDenom;
#call simplify();
endargument;



* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);


* move Chir to the left
#call moveChirR();
#call simplify();


.sort






id p1.ep1 = 0;
id p1.cep1 = 0;
id p2.ep2 = 0;
id p2.cep2 = 0;
id p2.ep1 = 0;
id p2.cep1 = 0;
id p1.ep2 = 0;
id p1.cep2 = 0;






id Ga(DumStr1?,p2) = -Ga(DumStr1,p1) + Ga(DumStr1,p3) + Ga(DumStr1,p4);
id p2 = -p1+p3+p4;


* move Chir to the left
#call moveChirR();
#call simplify();


repeat;
   #call commuteToLeft(p2);
   #call commuteToRight(p1); 
#call simplify();
endrepeat;


* Dirac algebra
repeat;
   #call commuteToLeft(p3);
   #call commuteToRight(p4);
   #call applyDiracEq(1);  
endrepeat;
#call simplify();

id p1.ep1 = 0;
id p1.cep1 = 0;
id p2.ep2 = 0;
id p2.cep2 = 0;
id p2.ep1 = 0;
id p2.cep1 = 0;
id p1.ep2 = 0;
id p1.cep2 = 0;



.sort


#do i=1,1
   multiply SpiStr(`i');
   repeat;
      id SpiStr(`i',?args1)*Ga?{ASpi,Ga,Chir,Spi}(`i',?args2) = SpiStr(`i',?args1,Ga(`i',?args2));
   endrepeat;
   id SpiStr(`i',?args) = SpiStr(?args);
#enddo






* now we conjugate the LO

.store
Off statistics;

Global [cc3,1] = [3,1];
Global [cc4,1] = [4,1];
Global [cc5,1] = [5,1];

id i_ = -i_;
id cI = -cI;

id SpiStr(?args) = SpiStr(reverse_(?args));

argument SpiStr;
   id ASpi(?args) = SpiDum(?args);
   id Spi(?args) = ASpi(?args);
   id SpiDum(?args) = Spi(?args);
   id Chir(DumSymb1?,lambda?) = Chir(DumSymb1,-lambda);
endargument;

#call conjugateEp();

*rename open Lorentz indices
#do i = 1,20,1
id Lor{`i'} = LorP{`i'};
argument;
   id Lor{`i'} = LorP{`i'};
   argument;
      id Lor{`i'} = LorP{`i'};
   endargument;
endargument;
#enddo


*rename open Color indices
#do i = 1,20,1
id Col{`i'}    = ColP{`i'};
id ColInt{`i'} = ColIntP{`i'};
id ColDum{`i'} = ColDumP{`i'};
id Glu{`i'}    = GluP{`i'};
id GluInt{`i'} = GluIntP{`i'};
id GluDum{`i'} = GluDumP{`i'};
argument;
   id Col{`i'} = ColP{`i'};
   id ColInt{`i'} = ColIntP{`i'};
   id ColDum{`i'} = ColDumP{`i'};
   id Glu{`i'}    = GluP{`i'};
   id GluInt{`i'} = GluIntP{`i'};
   id GluDum{`i'} = GluDumP{`i'};   
endargument;
#enddo


* conjugate color matrices
id SUNT(Glu1?,Col1?,Col2?) = SUNT(Glu1,Col2,Col1);

.sort

.store
Off statistics;






* square matrix elements, missing the (2*Re) here

Global [1,1,1] =  [1,1]*([cc3,1]+[cc4,1]+[cc5,1]);
Global [2,1,1] =  [2,1]*([cc3,1]+[cc4,1]+[cc5,1]);
.sort


argument;
id GluP1=Glu1;
id GluP2=Glu2;
id ColP3=Col3;
id ColP4=Col4;
endargument;

#call colorAlgebra2();



repeat;
   id SpiStr(?args1,Spi(DumSymb1?,+p1?,mX?))*SpiStr(ASpi(DumSymb1?,+p1?,mX?),?args2) = SpiStr(?args1,Ga(DumSymb1,p1),?args2)+mX*SpiStr(?args1,?args2);
   id SpiStr(ASpi(DumSymb1?,+p1?,mX?),?args,Spi(DumSymb1?,+p1?,mX?)) = SpiStr(Ga(DumSymb1,p1),?args)+mX*SpiStr(?args);
endrepeat;


#do i=1,10
   id once SpiStr(?args) = SpiStr{`i'}(?args);
   argument SpiStr{`i'};
      id Ga?{ASpi,Spi,Ga,Chir}(DumSymb1?,?args) = Ga(`i',?args);
   endargument;

   repeat;
      id SpiStr{`i'}(Ga?{ASpi,Spi,Ga,Chir}(?args1),?args2) = Ga(?args1)*SpiStr{`i'}(?args2);
   endrepeat;
   id SpiStr{`i'} = 1;
#enddo



#call moveChir()
#call simplify;



id p1.ep1 = 0;
id p1.cep1 = 0;
id p2.ep2 = 0;
id p2.cep2 = 0;
id p2.ep1 = 0;
id p2.cep1 = 0;
id p1.ep2 = 0;
id p1.cep2 = 0;


***#call evalTrace(1);
#call FORMTrace(1);
id e_(LorX?,LorY?,LorZ?,LorW?)=LeviCiv(LorX,LorY,LorZ,LorW);
****#call FORMTraceD(1);

#call simplify();


id p1.ep1 = 0;
id p1.cep1 = 0;
id p2.ep2 = 0;
id p2.cep2 = 0;
id p2.ep1 = 0;
id p2.cep1 = 0;
id p1.ep2 = 0;
id p1.cep2 = 0;

.sort

id LeviCiv(p1,p2,p3,p4)=0;


#call pullOutPolVec();

id cep1(DumLor1?DumId)*ep1(DumLor2?DumId) = - MeT(DumLor1,DumLor2) +  (p1(DumLor1)*p2(DumLor2)+p2(DumLor1)*p1(DumLor2))/Dot(p1,p2) ;
id cep2(DumLor1?DumId)*ep2(DumLor2?DumId) = - MeT(DumLor1,DumLor2) +  (p2(DumLor1)*p1(DumLor2)+p1(DumLor1)*p2(DumLor2))/Dot(p1,p2) ;

** comment: summing over DumLor indices that don't appear doesn't introduce artificial factors
   sum DumLor1;
   sum DumLor2;
   sum DumLor3;
   sum DumLor4;
   sum DumLor5;
   sum DumLor6;
   sum DumLor7;
   sum DumLor8;
   sum DumLor9;
   sum DumLor10;
   
#call simplify();

id LeviCiv(p1,p2,p3,p4)=0;

id DSTm4=0;


id p3.p4 = p1.p2 - MT^2;
id p1.p2 = shat/2;
id 1/Dot(p1,p2) = 2/shat;
id MT^2 = shat/4*(1-beta^2);

id p1.p3 = shat/4*(1 - beta*z );
id p1.p4 = shat/4*(1 + beta*z );
id p2.p3 = shat/4*(1 + beta*z );
id p2.p4 = shat/4*(1 - beta*z );

argument PropDenom;
    id p3.p4 = p1.p2 - MT^2;
    id p1.p2 = shat/2;
    id MT^2 = shat/4*(1-beta^2);

    id p1.p3 = shat/4*(1 - beta*z );
    id p1.p4 = shat/4*(1 + beta*z );
    id p2.p3 = shat/4*(1 + beta*z );
    id p2.p4 = shat/4*(1 - beta*z );
endargument;

id PropDenom(shat?)=1/shat;


