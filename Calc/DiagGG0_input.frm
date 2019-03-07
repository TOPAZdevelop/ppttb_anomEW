#define NumAmps "3" 

Global [1,1] = (i_*GS*DID(1)*ep1(Lor1)*ep2(Lor2)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*Spi(1, p4, -MT));
#define NumSM1 "1" 

Global [2,1] = (-(i_*DID(2)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM2 "1" 

Global [3,1] = (-(i_*DID(3)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM3 "1" 

#define AmpList "{[1\,1],[2\,1],[3\,1]}"

