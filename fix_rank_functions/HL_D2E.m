%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Hess L(x) eta_i by Algorithm 7 in Huang et al (2017)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HL_A, HL_B, HL_C] = HL_D2E(P, egradL, x, vU, sU, vV, sV, m, n, p, dim)
    
%% Compute HL_A,HL_B,HL_C with Hess L(x)[eta_i] = [U U_perp][A B'; C 0][V V_perp]'

p2 = p^2; mp = m*p; 

HL_A = cell(dim, 1); HL_B = cell(dim, 1); HL_C = cell(dim, 1);
HessL = @(t)func_HessL(P, egradL, x, vU, sU, vV, sV, m, n, p, mp, p2, t);
for i = 1:dim
    [tmpA, tmpB, tmpC] = HessL(i);
    HL_A{i} = tmpA;
    HL_B{i} = tmpB;
    HL_C{i} = tmpC;
end

%% Compute A,B,C with Hess L(x)[eta] = [U U_perp][A B'; C 0][V V_perp]' by combining Algorithm 7 in Huang et al (2017)
function [A, B, C] = func_HessL(P, egradL, x, vU, sU, vV, sV, m, n, p, mp, p2, t)
if t <= p2
    [row, col] = ind2sub([p,p], t);
    reta = x.U(:,row)*x.V(:,col)'; % compute reta x.U*Sdot*x.V'\in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    UU_perpPreta = fun_alpha(Preta, vU, sU, m, p);
    U_perpPreta = UU_perpPreta(p+1:m,:);
    C = U_perpPreta*x.V;

    VV_perpPretat = fun_alpha(Preta', vV, sV, n, p);
    V_perpPretat = VV_perpPretat(p+1:n,:);
    B = V_perpPretat*x.U;
elseif t > mp
    W = zeros(n-p,p);
    W(t-mp) = 1;
    W = [zeros(p,p);W];

    Vp = fun_beta(W, vV, sV, n, p); 
    reta = x.U*Vp'; %reta \in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    VV_perpPretat = fun_alpha(Preta', vV, sV, n, p);
    V_perpPretat = VV_perpPretat(p+1:n,:);
    B = V_perpPretat*x.U;

    Cbar = Preta + egradL*Vp*(x.invS.*x.V');
    UU_perpCbar = fun_alpha(Cbar, vU, sU, m, p);
    U_perpCbar = UU_perpCbar(p+1:m,:);
    C = U_perpCbar*x.V;
else
    K = zeros(m-p,p);
    K(t-p2) = 1;
    K = [zeros(p,p); K];

    Up = fun_beta(K, vU, sU, m, p);
    reta = Up*x.V'; %reta \in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    UU_perpPreta = fun_alpha(Preta, vU, sU, m, p);
    U_perpPreta = UU_perpPreta(p+1:m,:);
    C = U_perpPreta*x.V;

    Bbar = Preta + x.U*(x.invS.*Up')*egradL;
    VV_perpBbar = fun_alpha(Bbar', vV, sV, n, p);
    V_perpBbar = VV_perpBbar(p+1:n,:);
    B = V_perpBbar*x.U;
end
end

end