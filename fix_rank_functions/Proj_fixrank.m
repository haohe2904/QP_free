%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the projection from R^{m X n} \to T_x M.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PZ = Proj_fixrank(Z, x, vU, sU, vV, sV, m, n, p)
    % Compute U'*Z*V
    ZV = Z*x.V;
    PZ.A = x.U'*ZV;

    % Compute (U'*Z*V_perp)' 
    ZU = Z'*x.U;
    ZB = fun_alpha(ZU, vV, sV, n, p);
    PZ.B = ZB(p+1:n,:);

    % Compute U_perp'*Z*V             
    ZC = fun_alpha(ZV, vU, sU, m, p);
    PZ.C = ZC(p+1:m,:);
 end