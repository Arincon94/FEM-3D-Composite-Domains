function C=orthocl(E1,E2,E3,G4,G5,G6,v23,v13,v12);
% ORTHOCL - Constructs Cij constitutive matrix for Orthotropic materials.
%
%	Syntaxe: C=orthocl(E1,E2,E3,G23,G13,G12,v23,v13,v12)
%

% Copyright 2000 (c) Marcelo A. Trindade
% Last modified 12/05/2000

v32=v23/E2*E3; v31=v13/E1*E3; v21=v12/E1*E2; D=1-v12*v21-v23*v32-v31*v13-2*v21*v32*v13;

C=zeros(6);

C(1,1)=E1*(1  -v23*v32)/D;
C(1,2)=E1*(v21+v31*v23)/D; C(2,1)=C(1,2);
C(1,3)=E1*(v31+v21*v32)/D; C(3,1)=C(1,3);

C(2,2)=E2*(1  -v13*v31)/D;
C(2,3)=E2*(v32+v12*v31)/D; C(3,2)=C(2,3);

C(3,3)=E3*(1  -v21*v12)/D;

C(4,4)=G4;
C(5,5)=G5;
C(6,6)=G6;
