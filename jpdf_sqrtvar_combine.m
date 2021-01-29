
% Find mass posn and var jpdf's given 3 or more modes
% jpdf is >= 3D but slice through different choices of sv (sqrt of var).

%INPUTS:
    % nmodes       - number of modes
    % df           - fractional frequency jumps
    % sv           - square root of variance (2nd linear mass moment)
    % SIGMA        - noise covariance matrix
%OUTPUTS:
    % jpdf_ma - 2D matrix containting jpdf value versus mass and position

function jpdf_mu = jpdf_sqrtvar_combine(nmodes,df,sv,SIGMA)

global m Ph_alpha PH2 D1PH2 D2PH2 D3PH2

if isnan(df)% || nmodes < 3
    jpdf_mu = NaN(length(m));
    return;
end

SigmaN = SIGMA(1:nmodes,1:nmodes);
Si = inv(SigmaN);
vv = sv^2;

%constructing the mode shapes
h_mu1 = cell(nmodes,1);
ZZ = 0;

for i=1:nmodes
    h_mu1{i} = -m*(PH2(i,:)+.5*D2PH2(i,:)*vv)/(2*Ph_alpha(i));
end

for i=1:nmodes
    for j=1:nmodes
        h_mu2 = (h_mu1{i}-df(i)).*(h_mu1{j}-df(j));
        ZZ = ZZ + Si(i,j)*h_mu2;
    end
end


J1 = ( (PH2(1,:)+.5*D2PH2(1,:)*vv) .* ...
          ( ( (D1PH2(2,:)+.5*D3PH2(2,:)*vv) .* (.5*D2PH2(3,:)*sv) ) -...
            ( (D1PH2(3,:)+.5*D3PH2(3,:)*vv) .* (.5*D2PH2(2,:)*sv) ) ) );
J2 = ( (PH2(2,:)+.5*D2PH2(2,:)*vv) .* ...
          ( ( (D1PH2(1,:)+.5*D3PH2(1,:)*vv) .* (.5*D2PH2(3,:)*sv) ) -...
            ( (D1PH2(3,:)+.5*D3PH2(3,:)*vv) .* (.5*D2PH2(1,:)*sv) ) ) );
J3 = ( (PH2(3,:)+.5*D2PH2(3,:)*vv) .* ...
          ( ( (D1PH2(1,:)+.5*D3PH2(1,:)*vv) .* (.5*D2PH2(2,:)*sv) ) -...
            ( (D1PH2(2,:)+.5*D3PH2(2,:)*vv) .* (.5*D2PH2(1,:)*sv) ) ) );
        
J = (J1-J2+J3) / (Ph_alpha(1)*Ph_alpha(2)*Ph_alpha(3));

for i=4:nmodes
    J = J .* D2PH2(i,:) * sv / (2*Ph_alpha(i));
end

J = m.^(nmodes-1) .* abs(J);

jpdf_mu = (2*pi)^(-nmodes/2)*J.*exp(-ZZ/2);
  % answer off by constant factor due to x scale from derivatives

end






