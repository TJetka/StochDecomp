function R = JAKSTAT_rates(x, par, t, stimulus)
R =[
par(1)*stimulus*x(1)*exp(x(17)-par(8)/par(9));
par(2)*x(2)*x(2);
par(3)*x(3);
par(4)*x(4);
par(4)*x(5);
par(4)*x(6);
par(4)*x(7);
par(4)*x(8);
par(4)*x(9);
par(4)*x(10);
par(4)*x(11);
par(4)*x(12);
par(4)*x(13);
par(8);
par(9)*x(17);
];
 end
