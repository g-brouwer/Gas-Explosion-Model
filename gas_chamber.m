function gas_chamber_height = gas_chamber(r21, rhogz, rhos, N)

%this function (called in func_explode) calculates the gas chamber size using cap thickness and
%gas mass fraction as inputs. 

%polynomial coefficients
a=r21^3;
b=3*r21^2;
c=3*r21;
d=rhogz/(rhos*N);

%cubic polynomial
poly=[a b c -d]; 
%calculates roots of the cubic polynomial
r=roots(poly);
%find the real root
real_root=real(r(find(imag(r)==0)));

%gas chamber height: feeds into func_explode
gas_chamber_height=1/real_root;

end

