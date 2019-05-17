

* eval propagator functions
id PropDenom(pDumX?,mX?) = PropDenom(pDumX.pDumX-mX^2);

argument PropDenom;
#call simplify();
endargument;



* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);


* move Chir to the right
#call moveChirR();
#call simplify();




repeat;
#call commuteToLeft(q1);
#call simplify();
endrepeat;



******************* higher-rank tensor reduction *******************************
************************* method: cancel D0 ************************************

** comment: implmement a higher-rank reduction and show that higher-ranks cancel out 


* cancelling l^2
repeat;
     id (q1.q1)  *LoopDenom(q1,m0?, pDumX?,m1?, pDumY?,m2? ) = LoopDenom(pDumX, m1, pDumY,m2) + m0^2*LoopDenom(q1,m0, pDumX, m1, pDumY,m2);
endrepeat;


*
* checking that LoopDenom(p3+q1,MT,-p4+q1,MT)*q1.q1 multiplies q1-indepdent terms
* and
* LoopDenom(p3+q1,MT,-p4+q1,MT) multiplies q1-terms only linearly 
* before doing the shifts
* 
*Print;
*Bracket LoopDenom,q1.q1;
*.sort
*.end


*
* doing the shift for the q1.q1 term  (CHECK THIS AGAIN!) HERE WAS A BUG
* 
id DID(2)*LoopDenom(p3+q1,MT, -p4+q1,MT)*q1.q1 = DID(2)*( LoopDenom(q1,MT) + (MT^2-2*q1.p3+p3.p3)*LoopDenom(q1,MT, -p3-p4+q1,MT) );  
id DID(3)*LoopDenom(p3+q1, 0, -p4+q1, 0)*q1.q1 = DID(3)*( LoopDenom(q1, 0) + (0000-2*q1.p3+p3.p3)*LoopDenom(q1, 0, -p3-p4+q1, 0) );



* detecting other quadratic terms 
* 
id q1(LorX?)*Ga(DumStr1?,q1)*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = MaydayMaydayMayday;
id q1(LorX?)*q1(LorY?)*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = MaydayMaydayMayday;
id q1.pDumX?*Ga(DumStr1?,q1)*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = MaydayMaydayMayday;
id q1.pDumX?*q1.pDumY?*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = MaydayMaydayMayday;
id Ga(DumStr1?,q1)*Ga(DumStr2?,q1)*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = MaydayMaydayMayday;

* doing the shift for the term linear in q1
id Ga(DumStr1?,q1)*LoopDenom(p3+q1,MT?, -p4+q1,MT?) = Ga(DumStr1,q1-p3)*LoopDenom(q1,MT, -p3-p4+q1,MT);
*
id LoopDenom(p3+q1,MT?, -p4+q1,MT?)*(q1.pDumX?) = (q1.pDumX-p3.pDumX)*LoopDenom(q1,MT, -p3-p4+q1,MT);
*
id LoopDenom(p3+q1,MT?, -p4+q1,MT?) = LoopDenom(q1,MT, -p3-p4+q1,MT);






**************************************************


******************* higher-rank tensor reduction *******************************
************************* method: cancel D1 ************************************



*id q1.q1*LoopDenom(q1,m0?, pDumX?,m1?, pDumY?,m2? ) =  LoopDenom(q1,m0,pDumY,m2) - ( (pDumX.pDumX-q1.q1-2*(pDumX.q1-q1.q1)) +  2*(pDumX.q1-q1.q1)-m1^2)*LoopDenom(q1,m0,pDumX,m1,pDumY,m2);
*id q1.q1*LoopDenom(q1,m0?, pDumX?,m1?) =  LoopDenom(q1,m0) - ((pDumX.pDumX-q1.q1-2*(pDumX.q1-q1.q1))+2*(pDumX.q1-q1.q1)-m1^2)*LoopDenom(q1,m0,pDumX,m1);

*id q1.q1*LoopDenom(q1,m0?, pDumX?,m1?, pDumY?,m2? ) =  LoopDenom(q1,m0,pDumY,m2) - ( (pDumX.pDumX-q1.q1-2*(pDumX.q1-q1.q1)) +  2*(pDumX.q1-q1.q1)-m1^2)*LoopDenom(q1,m0,pDumX,m1,pDumY,m2);
*id q1.q1*LoopDenom(q1,m0?, pDumX?,m1?) =  LoopDenom(q1,m0) - ((pDumX.pDumX-q1.q1-2*(pDumX.q1-q1.q1))+2*(pDumX.q1-q1.q1)-m1^2)*LoopDenom(q1,m0,pDumX,m1);


**************************************************
* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);
id DID(Col1?)=1;

*Print;
*Bracket LoopDenom,q1.q1;
*.sort
*.end




id Ga(DumStr1?,p3) = -Ga(DumStr1,p4) + Ga(DumStr1,p1) + Ga(DumStr1,p2);
id p3 = -p4+p1+p2;


* Dirac algebra
repeat;
   #call commuteToLeft(p2);
   #call commuteToRight(p1);
   #call applyDiracEq(1);  
endrepeat;
#call simplify();

repeat;
   #call commuteToLeft(p3);
   #call commuteToRight(p4);
   #call applyDiracEq(2);  
endrepeat;
#call simplify();

.sort




* pull out loop momenta q1
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor1)*q1(DumLor1);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor2)*q1(DumLor2);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor3)*q1(DumLor3);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor4)*q1(DumLor4);
id once pDumX?FVec.q1 = pDumX(DumLor5)*q1(DumLor5);
id once pDumX?FVec.q1 = pDumX(DumLor6)*q1(DumLor6);
id once pDumX?FVec.q1 = pDumX(DumLor7)*q1(DumLor7);
id once pDumX?FVec.q1 = pDumX(DumLor8)*q1(DumLor8);
id once LeviCiv(q1,LorX?,LorY?,LorZ?) = LeviCiv(DumLor9,LorX,LorY,LorZ)*q1(DumLor9);


* identify tensor integrals (Denner's conventions)
multiply SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)*q1(LorY?DumId)*q1(LorZ?DumId) = LoopMom(LorW,LorX,LorY,LorZ)/SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)*q1(LorY?DumId)                = LoopMom(LorW,LorX,LorY)     /SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)                               = LoopMom(LorW,LorX)          /SIntDummy;
id q1(LorW?DumId)                                              = LoopMom(LorW)               /SIntDummy;
id SIntDummy                                                   = LoopMom(0);


id LoopDenom( q1,m0? ) * LoopMom(?args)                   
   = i_*Pi^2 * TI(1,?args,m0);
id LoopDenom( q1,m0?, pDumW?,m1? ) * LoopMom(?args)                   
   = i_*Pi^2 * TI(2,?args,pDumW-q1,m0,m1);
id LoopDenom( q1,m0?, pDumW?,m1?, pDumX?,m2? ) * LoopMom(?args) 
   = i_*Pi^2 * TI(3,?args,pDumW-q1,pDumX-q1,m0,m1,m2);
id LoopDenom( q1,m0?, pDumW?,m1?, pDumX?,m2?, pDumY?,m3? ) * LoopMom(?args) 
   = i_*Pi^2 * TI(4,?args,pDumW-q1,pDumX-q1,pDumY-q1,m0,m1,m2,m3);

   
   
id TI(1,0,m0?) 
 = SI(1,m0);
id TI(2,0,pDumW?,m0?,m1?) 
 = SI(2,pDumW,m0,m1);
id TI(3,0,pDumW?,pDumX?,m0?,m1?,m2?) 
 = SI(3,pDumW,pDumX,m0,m1,m2);
id TI(4,0,pDumW?,pDumX?,pDumY?,m0?,m1?,m2?,m3?) 
 = SI(4,pDumW,pDumX,pDumY,m0,m1,m2,m3);

* remove LoopMom(0) from the LO contribution
id LoopMom(0)=1;



* insert tensor decomposition ( scalar integrals must be identified first )
#include /mnt/pep/mschulze/lib/FeynArtsToForm/LorDec.frm




sum DumLor1;
sum DumLor2;
sum DumLor3;
sum DumLor4;
sum DumLor5;
sum DumLor6;
sum DumLor7;
sum DumLor8;
sum DumLor9;

#call simplify();
.sort 


* move Chir to the left
#call moveChirR();
#call simplify();



id Ga(DumStr1?,p3) = -Ga(DumStr1,p4) + Ga(DumStr1,p1) + Ga(DumStr1,p2);
id p3 = -p4+p1+p2;


* Dirac algebra
repeat;
   #call commuteToLeft(p2);
   #call commuteToRight(p1);
   #call applyDiracEq(1);  
endrepeat;
#call simplify();

repeat;
   #call commuteToLeft(p3);
   #call commuteToRight(p4);
   #call applyDiracEq(2);  
endrepeat;
#call simplify();

.sort





* now we conjugate the LO


#do i=1,2
   multiply SpiStr(`i');
   repeat;
      id SpiStr(`i',?args1)*Ga?{ASpi,Ga,Chir,Spi}(`i',?args2) = SpiStr(`i',?args1,Ga(`i',?args2));
   endrepeat;
   id SpiStr(`i',?args) = SpiStr(?args);
#enddo





.store
Off statistics;


Global [cc1,1] = [1,1];

id i_ = -i_;
id cI = -cI;

* removing the WFRC for this piece (higher order)
id dZfL=1;
id dZfR=1;

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



.store
Off statistics;








* square matrix elements

Global [1,1,1] = [1,1]*[cc1,1];
Global [2,1,1] = [2,1]*[cc1,1];
Global [3,1,1] = [3,1]*[cc1,1];
Global [4,1,1] = [4,1]*[cc1,1];


repeat;
   id SpiStr(?args1,Spi(DumSymb1?,+p1?,mX?))*SpiStr(ASpi(DumSymb1?,+p1?,mX?),?args2) = SpiStr(?args1,Ga(DumSymb1,p1),?args2)+mX*SpiStr(?args1,?args2);
endrepeat;
id SpiStr(ASpi(DumSymb1?,+p1?,mX?),?args,Spi(DumSymb1?,+p1?,mX?)) = SpiStr(Ga(DumSymb1,p1),?args)+mX*SpiStr(?args);




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



#call evalTrace(1);
#call evalTrace(2);

*#call FORMTraceD(1);
*#call FORMTraceD(2);

#call simplify();

id LeviCiv(p1,p2,p3,p4)=0;

