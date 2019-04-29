

* eval propagator functions
id PropDenom(pDumX?,mX?) = PropDenom(pDumX.pDumX-mX^2);

argument PropDenom;
#call simplify();
endargument;



* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);


id p2=p1;
argument;
  id p2=p1;
endargument;
id Spi(1,p1,MT) = 1;
id ASpi(1,p1,MT) =1;
#call simplify();


#call moveChirR();
#call commuteToLeft(p1);
#call commuteToLeft(q1);
#call simplify();





******************* higher-rank tensor reduction *******************************




* cancelling l^2
repeat;
     id (q1.q1)  *LoopDenom(q1,m0?, pDumX?,m1? ) = LoopDenom(pDumX, m1) + m0^2*LoopDenom(q1,m0, pDumX, m1);
endrepeat;


*Print;
*Bracket LoopDenom,q1.q1;
*.sort
*.end

* doing the shift for the term linear in q1: Ga(1,q1) 
id Ga(DumStr1?,q1)*LoopDenom(-p1+q1,MT?) = Ga(DumStr1,q1+p1)*LoopDenom(q1,MT);
id LoopDenom(-p1+q1,MT?) = LoopDenom(q1,MT);


* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);

* removing A_mu(M^2) integral because it is zero
id Ga(DumStr1?,q1)*LoopDenom(q1,MZ?) = 0;


*Print;
*Bracket LoopDenom,q1.q1;
*.sort
*.end


**************************************************



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


#call moveChirR();
#call commuteToLeft(p1);
#call simplify();



id Ga(1,p1)*Chir(1,+1) = Ga(1,p1,+1);
id Ga(1,p1)*Chir(1,-1) = Ga(1,p1,-1);

