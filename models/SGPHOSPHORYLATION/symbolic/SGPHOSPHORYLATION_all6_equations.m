function R = f(t,y,p)

R = [
     [ (p(1) - p(2)*y(1)) ];
     [ (p(3)*y(1) - p(4)*y(2) - p(5)*y(2) + p(6)*y(3)) ];
     [ (p(5)*y(2) - p(6)*y(3) - p(7)*y(3)) ];
     [ (-2*p(2)*y(4)) ];
     [ (p(6)*y(3) - 2*y(5)*(p(4) + p(5)) + 2*p(3)*y(7) + 2*p(6)*y(9)) ];
     [ (p(6)*y(3) - 2*y(6)*(p(6) + p(7)) + 2*p(5)*y(9)) ];
     [ (p(3)*y(4) - y(7)*(p(4) + p(5)) - p(2)*y(7) + p(6)*y(8)) ];
     [ (p(5)*y(7) - p(2)*y(8) - y(8)*(p(6) + p(7))) ];
     [ (p(5)*y(5) - y(9)*(p(6) + p(7)) - p(6)*y(3) - y(9)*(p(4) + p(5)) + p(3)*y(8) + p(6)*y(6)) ];
];
