#define NumAmps "5" 




* mass counter term insertions



Global [1,1] = (DID(10)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*(-(i_*GS*(Ga(1, Lor1)*Chir(1, -1) + Ga(1, Lor1)*Chir(1, 1))*SUNT(Glu1, Col3, ColInt5)))*(MT + Ga(1, p2 - p4))*(-(i_*(2*dMf1(3, 3) + MT*(dZfL1(3, 3, 3) + dZfR1(3, 3, 3)) - 2*dZfR1(3, 3, 3)*Ga(1, p2 - p4)*Chir(1, 1) + 2*dZfL1(3, 3, 3)*Ga(1, -p2 + p4)*Chir(1, -1)))/2)*(MT + Ga(1, p2 - p4))*(-(i_*GS*(Ga(1, Lor2)*Chir(1, -1) + Ga(1, Lor2)*Chir(1, 1))*SUNT(Glu2, ColInt5, Col4)))*Spi(1, p4, -MT));
#define NumSM10 "1" 

Global [2,1] = (DID(12)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*(-(i_*GS*(Ga(1, Lor2)*Chir(1, -1) + Ga(1, Lor2)*Chir(1, 1))*SUNT(Glu2, Col3, ColInt5)))*(MT + Ga(1, -p2 + p3))*(-(i_*(2*dMf1(3, 3) + MT*(dZfL1(3, 3, 3) + dZfR1(3, 3, 3)) + 2*dZfL1(3, 3, 3)*Ga(1, p2 - p3)*Chir(1, -1) - 2*dZfR1(3, 3, 3)*Ga(1, -p2 + p3)*Chir(1, 1)))/2)*(MT + Ga(1, -p2 + p3))*(-(i_*GS*(Ga(1, Lor1)*Chir(1, -1) + Ga(1, Lor1)*Chir(1, 1))*SUNT(Glu1, ColInt5, Col4)))*Spi(1, p4, -MT));
#define NumSM12 "1" 


*** tree amplitudes


Global [3,1] = (i_*GS*DID(3)*ep1(Lor1)*ep2(Lor2)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*Spi(1, p4, -MT));
#define NumSM10 "1" 

Global [4,1] = (-(i_*DID(4)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM11 "1" 

Global [5,1] = (-(i_*DID(5)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM12 "1"  



#define AmpList "{[1\,1],[2\,1],[3\,1],[4\,1],[5\,1]}"

