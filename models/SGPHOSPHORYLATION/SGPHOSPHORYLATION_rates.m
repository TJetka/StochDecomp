function R = PHOSPHORYLATION_rates(y,par,t,stimulus)
R=[par(1);
par(2)*y(1);
par(3)*y(1);
par(4)*y(2);
par(5)*y(2);
par(6)*y(3);
par(7)*y(3);];
end

