#define NumAmps "3" 


* 1:tree, 2:G0-vertex, 3:Gpm-vertex



Global [1,1] = (i_*DID(1)*MeT(Lor1, Lor2)*PropDenom(p3 + p4, 0)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(GluInt5, 8, Internal))*(ASpi(1, p2, -MU)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(GluInt5, Col2, Col1)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(GluInt5, Col2, Col1))*Spi(1, p1, MU))*(ASpi(2, p3, MT)*(-(i_*GS*(1+dZfL)*Ga(2, Lor2)*Chir(2, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*(1+dZfR)*Ga(2, Lor2)*Chir(2, 1)*SUNT(GluInt5, Col3, Col4))*Spi(2, p4, -MT));
#define NumSM1 "2" 


Global [2,1] = (DID(1)*LoopDenom(q1, MZ, p3 + q1, MT, -p4 + q1, MT)*MeT(Lor1, Lor2)*PropDenom(p3 + p4, 0)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(GluInt5, 8, Internal))*
(ASpi(1, p2, -MU)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(GluInt5, Col2, Col1)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(GluInt5, Col2, Col1))*Spi(1, p1, MU))*
(ASpi(2, p3, MT)*(   
     (EL*MT*Chir(2, -1))/(2*MW*SW) - (EL*MT*Chir(2, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(2,-q1)*Chir(2,-1) + C33phiu*Ga(2,-q1)*Chir(2,+1)  ) 
)*(MT + Ga(2, p3 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(GluInt5, Col3, Col4))*(MT + Ga(2, -p4 + q1))*(    
     (EL*MT*Chir(2, -1))/(2*MW*SW) - (EL*MT*Chir(2, 1))/(2*MW*SW) - EL/(2*SW*MW)*voL^2*( -2*C33phiq3*Ga(2,+q1)*Chir(2,-1) + C33phiu*Ga(2,+q1)*Chir(2,+1)  )   
)*Spi(2, p4, -MT));
#define NumSM2 "2" 

Global [3,1] = (DID(2)*LoopDenom(q1, MW, p3 + q1, MB, -p4 + q1, MB)*MeT(Lor1, Lor2)*PropDenom(p3 + p4, 0)*SumOver(Col1, 3, External)*SumOver(Col2, 3, External)*SumOver(Col3, 3, External)*SumOver(Col4, 3, External)*SumOver(GluInt5, 8, Internal))*
(ASpi(1, p2, -MU)*(-(i_*GS*Ga(1, Lor1)*Chir(1, -1)*SUNT(GluInt5, Col2, Col1)) - i_*GS*Ga(1, Lor1)*Chir(1, 1)*SUNT(GluInt5, Col2, Col1))*Spi(1, p1, MU))*
(ASpi(2, p3, MT)*(
      (i_*EL*MT*Chir(2, -1))/(MW*Sqrt2*SW) - (i_*EL*MB*Chir(2, 1))/(MW*Sqrt2*SW) + i_*EL/(MW*Sqrt2*SW)*voL^2*Ga(2,-q1)*Chir(2,-1)*C33phiq3
)*(MB + Ga(2, p3 + q1))*(-(i_*GS*Ga(2, Lor2)*Chir(2, -1)*SUNT(GluInt5, Col3, Col4)) - i_*GS*Ga(2, Lor2)*Chir(2, 1)*SUNT(GluInt5, Col3, Col4))*(MB + Ga(2, -p4 + q1))*(
     -((i_*EL*MB*Chir(2, -1))/(MW*Sqrt2*SW)) + (i_*EL*MT*Chir(2, 1))/(MW*Sqrt2*SW) -i_*EL/(MW*Sqrt2*SW)*voL^2*Ga(2,q1)*Chir(2, -1)*C33phiq3
)*Spi(2, p4, -MT));
#define NumSM3 "2" 

#define AmpList "{[1\,1]\,[2\,1],[3\,1]}"
