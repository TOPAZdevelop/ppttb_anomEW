#define NumAmps "9" 




*  loop amps with closed quark triangle



Global [1,1] = (DID(1)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*((EL*MQUGen5*Chir(2, -1))/(2*MW*SW) - (EL*MQUGen5*Chir(2, 1))/(2*MW*SW))*(MQUGen5 + Ga(2, p2 + q1))*(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt5, ColInt6) + i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt5, ColInt6))*(MQUGen5 + Ga(2, q1))*(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt6, ColInt5) + i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt6, ColInt5)));
#define NumSM1 "2" 

Global [2,1] = (DID(2)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*((EL*MQUGen5*Chir(2, -1))/(2*MW*SW) - (EL*MQUGen5*Chir(2, 1))/(2*MW*SW))*(MQUGen5 + Ga(2, p2 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt6, ColInt5)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt6, ColInt5))*(MQUGen5 + Ga(2, q1))*(-(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt5, ColInt6)) - i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt5, ColInt6)));
#define NumSM2 "2" 

Global [3,1] = (-(DID(3)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*MeT(Lor3, Lor4)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(i_*EL*IZ(3, -1)*Ga(1, Lor3)*Chir(1, -1) + i_*EL*IZ(3, 1)*Ga(1, Lor3)*Chir(1, 1))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*(-(i_*EL*IZ(3, 1)*Ga(2, Lor4)*Chir(2, -1)) - i_*EL*IZ(3, -1)*Ga(2, Lor4)*Chir(2, 1))*(MQUGen5 + Ga(2, p2 + q1))*(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt5, ColInt6) + i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt5, ColInt6))*(MQUGen5 + Ga(2, q1))*(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt6, ColInt5) + i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt6, ColInt5)));
#define NumSM3 "2" 

Global [4,1] = (-(DID(4)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*MeT(Lor3, Lor4)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(i_*EL*IZ(3, -1)*Ga(1, Lor3)*Chir(1, -1) + i_*EL*IZ(3, 1)*Ga(1, Lor3)*Chir(1, 1))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*(i_*EL*IZ(3, -1)*Ga(2, Lor4)*Chir(2, -1) + i_*EL*IZ(3, 1)*Ga(2, Lor4)*Chir(2, 1))*(MQUGen5 + Ga(2, p2 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt6, ColInt5)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt6, ColInt5))*(MQUGen5 + Ga(2, q1))*(-(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt5, ColInt6)) - i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt5, ColInt6)));
#define NumSM4 "2" 

Global [5,1] = (-(DID(5)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQDGen5, p2 + q1, MQDGen5, p2 - p3 - p4 + q1, MQDGen5)*MeT(Lor3, Lor4)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(i_*EL*IZ(3, -1)*Ga(1, Lor3)*Chir(1, -1) + i_*EL*IZ(3, 1)*Ga(1, Lor3)*Chir(1, 1))*Spi(1, p4, -MT))*((MQDGen5 + Ga(2, p2 - p3 - p4 + q1))*(-(i_*EL*IZ(4, 1)*Ga(2, Lor4)*Chir(2, -1)) - i_*EL*IZ(4, -1)*Ga(2, Lor4)*Chir(2, 1))*(MQDGen5 + Ga(2, p2 + q1))*(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt5, ColInt6) + i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt5, ColInt6))*(MQDGen5 + Ga(2, q1))*(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt6, ColInt5) + i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt6, ColInt5)));
#define NumSM5 "2" 

Global [6,1] = (-(DID(6)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQDGen5, p2 + q1, MQDGen5, p2 - p3 - p4 + q1, MQDGen5)*MeT(Lor3, Lor4)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(i_*EL*IZ(3, -1)*Ga(1, Lor3)*Chir(1, -1) + i_*EL*IZ(3, 1)*Ga(1, Lor3)*Chir(1, 1))*Spi(1, p4, -MT))*((MQDGen5 + Ga(2, p2 - p3 - p4 + q1))*(i_*EL*IZ(4, -1)*Ga(2, Lor4)*Chir(2, -1) + i_*EL*IZ(4, 1)*Ga(2, Lor4)*Chir(2, 1))*(MQDGen5 + Ga(2, p2 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt6, ColInt5)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt6, ColInt5))*(MQDGen5 + Ga(2, q1))*(-(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt5, ColInt6)) - i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt5, ColInt6)));
#define NumSM6 "2" 


*** tree amplitudes


Global [10,1] = (i_*GS*DID(10)*ep1(Lor1)*ep2(Lor2)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*Spi(1, p4, -MT));
#define NumSM10 "1" 

Global [11,1] = (-(i_*DID(11)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM11 "1" 

Global [12,1] = (-(i_*DID(12)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM12 "1"  



#define AmpList "{[1\,1],[2\,1],[3\,1],[4\,1],[5\,1],[6\,1],[10\,1],[11\,1],[12\,1]}"

