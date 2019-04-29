#define NumAmps "3" 


** 1: higgs
** 2: chi0
** 3: phi +/-


Global [1,1] = (-(DID(1)*IndexDelta(Col1, Col2)*LoopDenom(q1, MT, -p2 + q1, MH)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*
(-i_*EL/(2*MW*SW))*((MT*Chir(1, -1)) + (MT*Chir(1, 1)))*
(MT + Ga(1, q1))*
(-i_*EL/(2*MW*SW))*((MT*Chir(1, -1)) + (MT*Chir(1, 1)))*
Spi(1, p1, MT));
#define NumSM1 "1" 

Global [2,1] = (-(DID(2)*IndexDelta(Col1, Col2)*LoopDenom(q1, MT, -p2 + q1, MZ)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*
(-EL/(2*MW*SW))*((-MT*Chir(1, -1)) + (MT*Chir(1, 1)) + voL^2*(-2*C33phiq3*Ga(1,p1-q1)*Chir(1,-1)+C33phiu*Ga(1,p1-q1)*Chir(1,+1)) )*
(MT + Ga(1, q1))*
(-EL/(2*MW*SW))*((-MT*Chir(1, -1)) + (MT*Chir(1, 1)) + voL^2*(+2*C33phiq3*Ga(1,p1-q1)*Chir(1,-1)-C33phiu*Ga(1,p1-q1)*Chir(1,+1)) )*
Spi(1, p1, MT));
#define NumSM2 "1" 

Global [3,1] = (-(DID(3)*IndexDelta(Col1, Col2)*LoopDenom(q1, MB, -p2 + q1, MW)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*
(+i_*EL/(MW*Sqrt2*SW))*(+MT*Chir(1, -1) - MB*Chir(1, 1) + voL^2*Ga(1,p1-q1)*Chir(1,-1)*C33phiq3)*
(MB + Ga(1, q1))*
(-i_*EL/(MW*Sqrt2*SW))*(+MB*Chir(1, -1) - MT*Chir(1, 1) + voL^2*Ga(1,q1-p1)*Chir(1,-1)*C33phiq3)*
Spi(1, p1, MT));
#define NumSM3 "1" 
#define NumSM3 "1" 

#define AmpList "{[1\,1],[2\,1],[3\,1]}"

