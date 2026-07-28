%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Hess L(x) eta_i by Algorithm 7 in Huang et al (2017)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HL_A, HL_B, HL_C] = HL_D2E_mex(P, egradL, x, qr_U, tau_U, qr_V, tau_V, m, n, p, dim)
    
%% Compute HL_A,HL_B,HL_C with Hess L(x)[eta_i] = [U U_perp][A B'; C 0][V V_perp]'

p2 = p^2; mp = m*p; 

HL_A = cell(dim, 1); HL_B = cell(dim, 1); HL_C = cell(dim, 1);
HessL = @(t)func_HessL(P, egradL, x, qr_U, tau_U, qr_V, tau_V, m, n, p, mp, p2, t);
for i = 1:dim
    [tmpA, tmpB, tmpC] = HessL(i);
    HL_A{i} = tmpA;
    HL_B{i} = tmpB;
    HL_C{i} = tmpC;
end

%% Compute A,B,C with Hess L(x)[eta] = [U U_perp][A B'; C 0][V V_perp]' by combining Algorithm 7 in Huang et al (2017)
function [A, B, C] = func_HessL(P, egradL, x, qr_U, tau_U, qr_V, tau_V, m, n, p, mp, p2, t)
if t <= p2
    [row, col] = ind2sub([p,p], t);
    reta = x.U(:,row)*x.V(:,col)'; % compute reta x.U*Sdot*x.V'\in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    UU_perpPreta = apply_q_mex(Preta,qr_U, tau_U,  'T');
    U_perpPreta = UU_perpPreta(p+1:m,:);
    C = U_perpPreta*x.V;
    VV_perpPretat = apply_q_mex(Preta', qr_V, tau_V, 'T'); 
    V_perpPretat = VV_perpPretat(p+1:n,:);
    B = V_perpPretat*x.U;
elseif t > mp
    W = zeros(n-p,p);
    W(t-mp) = 1;
    W = [zeros(p,p);W];

    Vp = apply_q_mex(W,qr_V, tau_V, 'N');
    reta = x.U*Vp'; %reta \in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    VV_perpPretat = apply_q_mex(Preta',qr_V, tau_V,  'T');
    V_perpPretat = VV_perpPretat(p+1:n,:);
    B = V_perpPretat*x.U;

    Cbar = Preta + egradL*Vp*(x.invS.*x.V');
    UU_perpCbar = apply_q_mex(Cbar,qr_U, tau_U,  'T');
    U_perpCbar = UU_perpCbar(p+1:m,:);
    C = U_perpCbar*x.V;
else
    K = zeros(m-p,p);
    K(t-p2) = 1;
    K = [zeros(p,p); K];

    Up = apply_q_mex(K,qr_U, tau_U,  'N'); 
    reta = Up*x.V'; %reta \in T_xM
    Preta = P.*reta;
    A = x.U'*Preta*x.V;

    UU_perpPreta = apply_q_mex(Preta,qr_U, tau_U,  'T'); 
    U_perpPreta = UU_perpPreta(p+1:m,:);
    C = U_perpPreta*x.V;

    Bbar = Preta + x.U*(x.invS.*Up')*egradL;
    VV_perpBbar = apply_q_mex(Bbar',qr_V, tau_V,  'T'); 
    V_perpBbar = VV_perpBbar(p+1:n,:);
    B = V_perpBbar*x.U;
end
end

end