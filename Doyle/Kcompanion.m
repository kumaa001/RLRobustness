function [Ak,Bk,Ck] = Kcompanion(p)

% Second-order system parameterized in companion form with 4 parameters.
Ak = [0 p(1); 1 p(2)];
Bk = [1; 0];
Ck = [p(3) p(4)];
