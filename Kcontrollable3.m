function [Ak,Bk,Ck] = Kcontrollable3(p)

% Third-order system parameterized in controllable canonical form
Ak = [0 1 0; 0 0 1; p(1) p(2) p(3)];
Bk = [0; 0; 1];
Ck = [p(4) p(5) p(6)];

