#define NumAmps "12" 




*  loop amps with charged goldstone Gpm

Global [1,1] = (GS*DID(1)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MW, p3 + q1, MB, -p4 + q1, MB)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, p3 + q1))*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*(MB + Ga(1, -p4 + q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p4, -MT));
#define NumSM1 "1" 

Global [2,1] = (-(DID(2)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, -p2 + q1, MB, -p4 + q1, MW)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, p2 - q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*(MB + Ga(1, -q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p4, -MT));
#define NumSM2 "1" 

Global [3,1] = (-(DID(3)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, -p4 + q1, MW, p2 - p3 - p4 + q1, MB)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, -p2 + p3 + p4 - q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*(MB + Ga(1, -q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p4, -MT));
#define NumSM3 "1" 

Global [4,1] = (-(DID(4)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, -p2 + q1, MB, -p3 + q1, MW)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MB + Ga(1, -p2 + q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM4 "1" 

Global [5,1] = (-(DID(5)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, -p3 + q1, MW, p2 - p3 - p4 + q1, MB)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MB + Ga(1, p2 - p3 - p4 + q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM5 "1" 

Global [6,1] = (-(DID(6)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, p2 + q1, MB, p2 - p4 + q1, MW, p2 - p3 - p4 + q1, MB)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, -p2 + p3 + p4 - q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MB + Ga(1, -q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*(MB + Ga(1, -p2 - q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p4, -MT));
#define NumSM6 "1" 

Global [7,1] = (-(DID(7)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, p2 + q1, MB, p2 - p3 + q1, MW, p2 - p3 - p4 + q1, MB)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, p2 + q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MB + Ga(1, q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*(MB + Ga(1, p2 - p3 - p4 + q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*Spi(1, p4, -MT));
#define NumSM7 "1" 

Global [8,1] = (-(DID(8)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, p2 - p4 + q1, MW)*PropDenom(-p2 + p4, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, -q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM8 "1" 

Global [9,1] = (-(DID(9)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MB, p2 - p3 + q1, MW)*PropDenom(p2 - p3, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*((i_*EL*MT*Chir(1, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(1, 1))/(MW*Sqrt2*SW))*(MB + Ga(1, q1))*(-((i_*EL*MB*Chir(1, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(1, 1))/(MW*Sqrt2*SW))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM9 "1" 




*** tree amplitudes


Global [10,1] = (i_*GS*DID(10)*ep1(Lor1)*ep2(Lor2)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*Spi(1, p4, -MT));
#define NumSM10 "1" 

Global [11,1] = (-(i_*DID(11)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM11 "1" 

Global [12,1] = (-(i_*DID(12)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM12 "1"  



#define AmpList "{[1\,1],[2\,1],[3\,1],[4\,1],[5\,1],[6\,1],[7\,1],[8\,1],[9\,1],[10\,1],[11\,1],[12\,1]}"

