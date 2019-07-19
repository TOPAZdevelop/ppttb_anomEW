#define NumAmps "15" 




*  loop amps with neutral goldstone G0


Global [1,1] = (GS*DID(1)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MZ, p3 + q1, MT, -p4 + q1, MT)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*
     ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,-q1)*Chir(1,-1) + C33phiu*Ga(1,-q1)*Chir(1,+1)  ) )*
(MT + Ga(1, p3 + q1))*(-(i_*GS*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*(MT + Ga(1, -p4 + q1))*
    ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,+q1)*Chir(1,-1) + C33phiu*Ga(1,+q1)*Chir(1,+1)  ) )*
Spi(1, p4, -MT));
#define NumSM1 "1" 

Global [2,1] = (-(DID(2)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, -p2 + q1, MT, -p4 + q1, MZ)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, p2 - q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*(MT + Ga(1, -q1))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT));
#define NumSM2 "1" 

Global [3,1] = (-(DID(3)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, -p4 + q1, MZ, -p1 + q1, MT)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*
       ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,q1-p4)*Chir(1,-1) + C33phiu*Ga(1,q1-p4)*Chir(1,+1)  ) )*
(MT + Ga(1, -p2 + p3 + p4 - q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*(MT + Ga(1, -q1))*
       ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,p4-q1)*Chir(1,-1) + C33phiu*Ga(1,p4-q1)*Chir(1,+1)  ) )*
Spi(1, p4, -MT));
#define NumSM3 "1" 

Global [4,1] = (-(DID(4)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, -p2 + q1, MT, -p3 + q1, MZ)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*
       ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,-q1+p3)*Chir(1,-1) + C33phiu*Ga(1,-q1+p3)*Chir(1,+1)  ) )*
(MT + Ga(1, q1))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + q1))*
       ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,-p3+q1)*Chir(1,-1) + C33phiu*Ga(1,-p3+q1)*Chir(1,+1)  ) )*       
(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM4 "1" 

Global [5,1] = (-(DID(5)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, -p3 + q1, MZ, p2 - p3 - p4 + q1, MT)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, q1))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p3 - p4 + q1))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM5 "1" 

Global [6,1] = (DID(6)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*((EL*MQUGen5*Chir(2, -1))/(2*MW*SW) - (EL*MQUGen5*Chir(2, 1))/(2*MW*SW))*(MQUGen5 + Ga(2, p2 + q1))*(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt5, ColInt6) + i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt5, ColInt6))*(MQUGen5 + Ga(2, q1))*(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt6, ColInt5) + i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt6, ColInt5)));
#define NumSM6 "2" 

Global [7,1] = (DID(7)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, p2 + q1, MQUGen5, p2 - p3 - p4 + q1, MQUGen5)*PropDenom(-p3 - p4, MZ)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(ColInt6, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, p2 - p3 - p4 + q1))*((EL*MQUGen5*Chir(2, -1))/(2*MW*SW) - (EL*MQUGen5*Chir(2, 1))/(2*MW*SW))*(MQUGen5 + Ga(2, p2 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(Glu2, ColInt6, ColInt5)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(Glu2, ColInt6, ColInt5))*(MQUGen5 + Ga(2, q1))*(-(i_*GS*Ga(2, Lor1)*Chir(2, -1)*SUNT(Glu1, ColInt5, ColInt6)) - i_*GS*Ga(2, Lor1)*Chir(2, 1)*SUNT(Glu1, ColInt5, ColInt6)));
#define NumSM7 "2" 


* removed chiral gluon coupling
Global [8,1] = (-(DID(8)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, p2 + q1, MT, p2 - p4 + q1, MZ, p2 - p3 - p4 + q1, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, -p2 + p3 + p4 - q1))*(-(i_*GS*Ga(1, Lor1)*SUNT(Glu1, Col3, ColInt5)) )*(MT + Ga(1, -q1))*(-(i_*GS*Ga(1, Lor2)*SUNT(Glu2, ColInt5, Col4)) )*(MT + Ga(1, -p2 - q1))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT));
#define NumSM8 "1" 

Global [9,1] = (-(DID(9)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, p2 + q1, MT, p2 - p3 + q1, MZ, -p1 + q1, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*
     ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,-q1-p2+p3)*Chir(1,-1) + C33phiu*Ga(1,-q1-p2+p3)*Chir(1,+1)  ) )*
(MT + Ga(1, p2 + q1))*(-(i_*GS*Ga(1, Lor2)*SUNT(Glu2, Col3, ColInt5)) )*(MT + Ga(1, q1))*(-(i_*GS*Ga(1, Lor1)*SUNT(Glu1, ColInt5, Col4)) )*(MT + Ga(1, -p1 + q1))*
     ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,+q1+p2-p3)*Chir(1,-1) + C33phiu*Ga(1,+q1+p2-p3)*Chir(1,+1)  ) )*     
Spi(1, p4, -MT));
#define NumSM9 "1" 

Global [10,1] = (-(GS*DID(10)*ep1(Lor1)*ep2(Lor2)*IndexDelta(Col3, Col4)*LoopDenom(q1, MQUGen5, -p3 - p4 + q1, MQUGen5)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(-p3 - p4, MZ)*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Gen5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5)))*(ASpi(1, p3, MT)*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*Spi(1, p4, -MT))*((MQUGen5 + Ga(2, -q1))*((EL*MQUGen5*Chir(2, -1))/(2*MW*SW) - (EL*MQUGen5*Chir(2, 1))/(2*MW*SW))*(MQUGen5 + Ga(2, p3 + p4 - q1))*(-(i_*GS*Ga(2, Lor4)*Chir(2, -1)*SUNT(GluInt5, ColInt5, ColInt5)) - i_*GS*Ga(2, Lor4)*Chir(2, 1)*SUNT(GluInt5, ColInt5, ColInt5)));
#define NumSM10 "2" 

Global [11,1] = (-(DID(11)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, p2 - p4 + q1, MZ)*PropDenom(-p2 + p4, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, -q1))*((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW))*(MT + Ga(1, p2 - p4))*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM11 "1" 

Global [12,1] = (-(DID(12)*ep1(Lor1)*ep2(Lor2)*LoopDenom(q1, MT, p2 - p3 + q1, MZ)*PropDenom(p2 - p3, MT)^2*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*
      ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,-q1-p2+p3)*Chir(1,-1) + C33phiu*Ga(1,-q1-p2+p3)*Chir(1,+1)  ) )*
(MT + Ga(1, q1))*
      ((EL*MT*Chir(1, -1))/(2*MW*SW) - (EL*MT*Chir(1, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(1,+q1+p2-p3)*Chir(1,-1) + C33phiu*Ga(1,+q1+p2-p3)*Chir(1,+1)  ) )*     
(MT + Ga(1, -p2 + p3))*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM12 "1" 



*** tree amplitudes


Global [13,1] = (i_*GS*DID(13)*ep1(Lor1)*ep2(Lor2)*MeT(Lor3, Lor4)*(MeT(Lor1, Lor2)*(-p1(Lor3) + p2(Lor3)) + MeT(Lor2, Lor3)*(-p2(Lor1) - p3(Lor1) - p4(Lor1)) + MeT(Lor1, Lor3)*(p1(Lor2) + p3(Lor2) + p4(Lor2)))*PropDenom(p3 + p4, 0)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)*SumOver(GluInt5, 8, Internal)*SUNF(Glu1, Glu2, GluInt5))*(ASpi(1, p3, MT)*(-(i_*GS*(1+dZfL)*Ga(1, Lor4)*Chir(1, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*(1+dZfR)*Ga(1, Lor4)*Chir(1, 1)*SUNT(GluInt5, Col3, Col4))*Spi(1, p4, -MT));
#define NumSM10 "1" 

Global [14,1] = (-(i_*DID(14)*ep1(Lor1)*ep2(Lor2)*PropDenom(-p2 + p4, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*(1+dZfL)*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, Col3, ColInt5)) - i_*GS*(1+dZfR)*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, Col3, ColInt5))*(MT + Ga(1, p2 - p4))*(-(i_*GS*(1+dZfL)*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, ColInt5, Col4)) - i_*GS*(1+dZfR)*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM11 "1" 

Global [15,1] = (-(i_*DID(15)*ep1(Lor1)*ep2(Lor2)*PropDenom(p2 - p3, MT)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(ColInt5, 3, Internal)*SumOver(Glu1, 8, External)*SumOver(Glu2, 8, External)))*(ASpi(1, p3, MT)*(-(i_*GS*(1+dZfL)*Ga(1, Lor2)*Chir(1, -1)*SUNT(Glu2, Col3, ColInt5)) - i_*GS*(1+dZfR)*Ga(1, Lor2)*Chir(1, 1)*SUNT(Glu2, Col3, ColInt5))*(MT + Ga(1, -p2 + p3))*(-(i_*GS*(1+dZfL)*Ga(1, Lor1)*Chir(1, -1)*SUNT(Glu1, ColInt5, Col4)) - i_*GS*(1+dZfR)*Ga(1, Lor1)*Chir(1, 1)*SUNT(Glu1, ColInt5, Col4))*Spi(1, p4, -MT));
#define NumSM12 "1"  



#define AmpList "{[1\,1],[2\,1],[3\,1],[4\,1],[5\,1],[6\,1],[7\,1],[8\,1],[9\,1],[10\,1],[11\,1],[12\,1],[13\,1],[14\,1],[15\,1]}"

