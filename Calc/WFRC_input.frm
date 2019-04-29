#define NumAmps "3" 

Global [1,1] = (-(DID(1)*IndexDelta(Col1, Col2)*LoopDenom(q1, MT, -p2 + q1, MH)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*(-(i_*EL*MT*Chir(1, -1))/(2*MW*SW) - (i_*EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, q1))*(-(i_*EL*MT*Chir(1, -1))/(2*MW*SW) - (i_*EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p1, MT));
#define NumSM1 "1" 

Global [2,1] = (-(DID(2)*IndexDelta(Col1, Col2)*LoopDenom(q1, MT, -p2 + q1, MZ)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, q1))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p1, MT));
#define NumSM2 "1" 

Global [3,1] = (-(DID(3)*IndexDelta(Col1, Col2)*LoopDenom(q1, MB, -p2 + q1, MW)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)))*(ASpi(1, p2, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p1, MT));
#define NumSM3 "1" 

#define AmpList "{[1\,1],[2\,1],[3\,1]}"

