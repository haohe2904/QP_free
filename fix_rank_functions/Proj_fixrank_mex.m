%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the projection from R^{m X n} \to T_x M.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PZ = Proj_fixrank_mex(Z, x, qr_U, tau_U, qr_V, tau_V, m, n, p)
    % Compute U'*Z*V
    ZV = Z*x.V;
    PZ.A = x.U'*ZV;

    % Compute (U'*Z*V_perp)' 
    ZU = Z'*x.U;
    ZB = apply_q_mex(ZU, qr_V, tau_V, 'T'); 
    PZ.B = ZB(p+1:n,:);

    % Compute U_perp'*Z*V             
    ZC = apply_q_mex(ZV, qr_U, tau_U, 'T'); 
    PZ.C = ZC(p+1:m,:);
 end