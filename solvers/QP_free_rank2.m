function [xfinal, costfinal, residual,stopreason] = QP_free_rank2(m, n, k, P, A, eqindices, ineqindices, x0, options)

    if isfield(options,'tolKKTres');             tol                  = options.tolKKTres;              end
    if isfield(options,'verbosity');             printyes             = options.verbosity;              end
    if isfield(options,'maxiter');               maxiter              = options.maxiter;                end
    if isfield(options,'mineigval_correction');  mineigval_correction = options.mineigval_correction;   end
    
    %% Initialization
    Time = tic;
    stopreason = 0;
    mn = m*n; mk = m*k;
    x = x0; dS = diag(x.S); 
    x.invS = 1./dS;
    X = x.U*(x.V'.*dS);  % compute x.U*x.S*x,V';
    
    % Papameters setting
    nu = 2.3; tau = 0.75; varrho = 2.4; kappa = 0.55;
    sigma = 0.45; varsigma = 0.5; gamma_max = 50;
    rho = 2; rho_u = 1.5; c1 = 0.5; c2 = 0.5; c3 = 0.5;
    cor_tol = 1e-5; ls_max_steps = 50;
    
    % Problem informations
    dim = (m+n-k)*k;
    ineq_num = length(ineqindices); eq_num = length(eqindices);
    constraint_number = ineq_num + eq_num;
    dimsystem = dim + constraint_number;
    
    % Initial Lagrange multiplier vectors and gammas
    mus = ones(ineq_num, 1); 
    lambdas = ones(eq_num, 1);
    gammas = 0.1*ones(constraint_number, 1);
    
    % Store the vector expression of Riemannian gradients of constriant
    Ineqgrads = zeros(dim,ineq_num); 
    Eqgrads = zeros(dim,eq_num);
    
    % Store cx which is the constriant function values
    cx = zeros(constraint_number,1);
    ls_cx = zeros(constraint_number,1);
    
    %% Make alpha and beta and delta
    for i = 1:constraint_number
        if i <= ineq_num
            cx(i) = -X(ineqindices(i));
        else 
            cx(i) = X(eqindices(i-ineq_num)) - A(eqindices(i-ineq_num));
        end
    end
    
    % Compute objective value
    PXA = P.*(X-A);
    xcost = 0.5*norm(PXA, 'fro')^2;
    
    % Update penalty parameter flag
    rho_flag = 0; ls_tag =0;

   % base = cell(dim,1);

    %M = fixedrankembeddedfactory(m, n, k); % m: # of rows, n: # of cols, k: rank
    
    for iter = 1:maxiter
        

        % x.Upe = null(x.U'); x.Vpe = null(x.V');
        % for i = 1:k^2
        %     [row, col] = ind2sub([k,k], i);
        %     base{i} = x.U(:,row)*x.V(:,col)';
        %     base{i} = x.U(:,row)*x.Vpe(:,col)';
        %     base{i} = x.Upe(:,row)*x.V(:,col)';
        % end
        % 
        % for i = 1:mk-k^2
        %     %[d1,d2] = size(x.Upe)
        %     [row, col] = ind2sub([m-k,k], i);
        %     base{i+k^2} = x.Upe(:,row)*x.V(:,col)';
        % end
        % 
        % for i = 1:n*k-k^2
        %     [row, col] = ind2sub([k,n-k], i);
        %     base{i+mk} = x.U(:,row)*x.Vpe(:,col)';
        % end
        
        %% Check whether the iteration point information needs to be updated
        if rho_flag == 0
            cx_gam = (cx.^2+gammas.^2).^0.5;
            alpha = cx./cx_gam + 1;
            beta = (1-gammas./cx_gam).^0.5;
            %delta = -beta./cx;
            delta = (1./(cx_gam.*(gammas+cx_gam))).^.5;
       
            % Compute v,s at Q by Algorithm 3 in Huang et al (2017), where Q is Orthogonal matrix.
            [vU, sU] = fun_vs(x.U, m, k);
            [vV, sV] = fun_vs(x.V, n, k);
        
            %% Make the matrix representation for the Riemannian gradients of constriants
            for i = 1:ineq_num
                [row, col] = ind2sub([m,n], ineqindices(i));    
        
                % Compute x.U'*egrad c_ineqs(x)*x.V
                ineqA = x.U(row,:)'*x.V(col,:);

                % G0 = zeros(m,n);
                % G0(ineqindices(i))=1;
                % G1 = x.U*x.U'*G0*x.V*x.V'+(eye(m)-x.U*x.U')*G0*x.V*x.V'+x.U*x.U'*G0*(eye(n)-x.V*x.V')
        
                % Compute (x.U'*egrad c_ineqs(x)*x.V_perp)'
                B0 = zeros(n,k);  
                B0(col,:) = x.U(row,:);             %   egrad c_ineqs(x)'*x.U
                ineqB = fun_alpha(B0, vV, sV, n, k);
                ineqB = ineqB(k+1:n,:);
        
                % Compute x.U_perp'*egrad c_ineqs(x)*x.V
                C0 = zeros(m,k);  
                C0(row,:) = x.V(col,:);               
                ineqC = fun_alpha(C0, vU, sU, m, k);
                ineqC = ineqC(k+1:m,:);

                % ABC = [ineqA ineqB'; ineqC zeros(m-k,n-k)];
                % 
                % G22 = fun_beta(ABC, vU, sU, m, k);
                % G2 = fun_beta(G22', vV, sV, n, k);
                % G2 = G2'
                
                % Make ineqgrads
                Ineqgrads(:,i) = -[ineqA(:); ineqC(:); ineqB(:)];

                % egardc = zeros(m,n);
                % egardc(row,col) = -1;
                % gradc = zeros(dim,1);
                % for l = 1:dim
                %     gradc(l)=base{l}(:)'*egardc(:);
                % end
                % 
                % Ineqgrads(:,i) = gradc;

                % G1 = zeros(m,n);
                % for m = 1:dim
                %     G1 = G1 + gradc(m)*base{m};
                % end
                % G1

            end
        
            for i = 1:eq_num
                [row, col] = ind2sub([m,n], eqindices(i));
        
                % Compute x.U'*egrad c_eqs(x)*x.V
                eqA = x.U(row,:)'*x.V(col,:);
        
                % Compute (x.U'*egrad c_eqs(x)*x.V_perp)'
                B0 = zeros(n,k);  
                B0(col,:) = x.U(row,:);               
                eqB = fun_alpha(B0, vV, sV, n, k);
                eqB = eqB(k+1:n,:);
        
                % Compute x.U_perp'*egrad c_eqs(x)*x.V
                C0 = zeros(m,k);  
                C0(row,:) = x.V(col,:);               
                eqC = fun_alpha(C0, vU, sU, m, k);
                eqC = eqC(k+1:m,:);
                
                % Make eqgrads
                Eqgrads(:,i) = [eqA(:); eqC(:); eqB(:)];


                % egardc = zeros(m,n);
                % egardc(row,col) = 1;
                % gradc = zeros(dim,1);
                % for l = 1:dim
                %     gradc(l)=base{l}(:)'*egardc(:);
                % end
                % 
                % Eqgrads(:,i) = gradc;
            end
            
            Gradc = [Ineqgrads Eqgrads];
            aGradct = alpha.*(Gradc');
            Dbeta = diag(-sqrt(2).*beta);
            
            %% Compute Euclidean and Riemannian gradient (in vector representation) of Lagrangian function
            egradf = PXA;
            gradf = Proj_fixrank(egradf, x, vU, sU, vV, sV, m, n, k);
            vgradf = [gradf.A(:); gradf.C(:); gradf.B(:)];
            egradLag = gradLagrangian(egradf, mus, lambdas, ineqindices, eqindices, ineq_num, eq_num, m, n);
            gradLag = Proj_fixrank(egradLag, x, vU, sU, vV, sV, m, n, k);
            vgradLag = [gradLag.A(:); gradLag.C(:); gradLag.B(:)];

            % ABC = [gradLag.A gradLag.B'; gradLag.C zeros(m-k,n-k)];
            % 
            % G22 = fun_beta(ABC, vU, sU, m, k);
            % G2 = fun_beta(G22', vV, sV, n, k);
            % G2 = G2'
            % 
            % G1 = M.egrad2rgrad(x, egradLag)


            % main_e_dir = vgradLag;
            % Sdot_main = reshape(main_e_dir(1:k^2), [k,k]); 
            % K_main = reshape(main_e_dir(k^2+1:mk), [m-k,k]);
            % W_mian = reshape(main_e_dir(mk+1:end), [n-k,k]);
            % 
            % main_dir.S = Sdot_main;
            % if m <= n
            %     M1_main = fun_beta([Sdot_main; K_main], vU, sU, m, k);
            %     M2_mian = fun_beta([zeros(k); W_mian], vV, sV, n, k);
            %     main_dir.Up = M1_main - x.U*Sdot_main;
            %     main_dir.Vp = M2_mian;
            % else
            %     M1_main = fun_beta([zeros(k); K_main], vU, sU, m, k);
            %     M2_mian = fun_beta([Sdot_main'; W_mian], vV, sV, n, k);
            %     main_dir.Up = M1_main;
            %     main_dir.Vp = M2_mian - x.V*Sdot_main';
            % end
            % G1 = x.U*main_dir.S*x.V'+main_dir.Up*x.V'+x.U*main_dir.Vp'
        else
            d_eqs = zeros(mn, 1);
            d_eqs(eqindices) = diff_rho;
            egradLag = egradLag - reshape(d_eqs, [m, n]);
            vgradLag = vgradLag - diff_rho*sum(Eqgrads, 2);
        end
    
        gradLagnorm = norm(vgradLag, 2);
    
        % Compute KKT residual and check the stopping conditions
        KKT_res = KKT_residual(X, A, mus, gradLagnorm, eqindices, ineqindices, ineq_num, eq_num);
        if KKT_res <= tol
            %fprintf('KKT Residual tolerance reached\n');
            stopreason = 1;
            break;
        end
    
        ttime = toc(Time);
        if printyes >= 3
            if iter == 1
                fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, xcost, KKT_res, ttime);
                kkt_flag0 = floor(log10(KKT_res));
            else
                kkt_flag1 = floor(log10(KKT_res));
                if kkt_flag1 < kkt_flag0
                    fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, xcost, KKT_res, ttime);
                    kkt_flag0 = kkt_flag1;
                end
            end
        elseif printyes >= 2
            fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, xcost, KKT_res, ttime);
        elseif printyes >= 1
            if mod(iter, 100) == 0 && iter ~= 0
                fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, xcost, KKT_res, ttime);
            end
        end
    
        %% Compute A,B,C with Hess L(x)[eta_i] = [U U_perp][A B'; C 0][V V_perp]', where {\eta_i}_{i\in [dim]} is the basis of T_xM
        [HL_A, HL_B, HL_C] = HL_D2E(P, egradLag, x, vU, sU, vV, sV, m, n, k, dim); 
        HL = HL_matrix(HL_A, HL_B, HL_C, dim, m, k);
        [V,T] = schur(HL);
        for i = 1 : dim
            if T(i,i) < 1e-5  
                T(i,i) = mineigval_correction;
            end
        end
        HL = V*T*V'; HL = 0.5*(HL.'+HL);
    
        % Make coefficient matrix B
       % HL = eye(dim);
        B = [HL Gradc; aGradct Dbeta];
     
        % Make the vector expression of gradpen
        vgradpen = vgradf - rho*sum(Eqgrads,2);

        % gradeqc = zeros(m,n);
        % gradeqc(eqindices) = 1;
        % gradpen = egradf - rho*gradeqc;
        % vgradpen = zeros(dim,1);
        % for l = 1:dim
        %     vgradpen(l)=base{l}(:)'*gradpen(:);
        % end
    
        %% Compute the direction and Lagrange multipliers by solving three linear system with LU decomposition.
        [L, U, Q] = lu(B); 

        % Solve the first linear system.
        b1 = [-vgradpen; zeros(constraint_number,1)];
        %dir_multipliers0 = B\b1;
        b1 = Q*b1;
        dir_multipliers1 = linearsolve(L, U, b1, dimsystem); 
        dir1 = dir_multipliers1(1:dim);
        multipliers1 = dir_multipliers1(dim+1:dimsystem);
        m_ineq1 = multipliers1(1:ineq_num);
        m_eq1 = multipliers1(ineq_num+1:end);
    
        min_ineq1 = min(m_ineq1); min_eq1 = min(m_eq1); 
        dir1_norm = norm(dir1,2);
    
        % Update rho, a penalty parameter, if necessary.
        if  dir1_norm <= c1 && min_eq1 < c2 && min_eq1 >= -c3 && min_ineq1 >= -c3
            diff_rho = (rho_u-1)*rho;
            rho = rho + diff_rho;
            rho_flag = 1;
            continue
        else
            rho_flag = 0;
        end
        
        % Solve the second linear system.
        al_lam3 = alpha.*(min(multipliers1,0).^3);
        b2 = [-vgradpen; al_lam3];
        b2 = Q*b2;
        dir_multipliers2 = linearsolve(L, U, b2, dimsystem); 
        dir2 = dir_multipliers2(1:dim);
        multipliers2 = dir_multipliers2(dim+1:dimsystem);
        dir2_norm_nu = norm(dir2,2)^nu;
    
        % If the norm of direction2 is too smaller, then we don't need to 
        % solve the thrid linear system. And using the direction2 as a master
        % search direction.
        if dir2_norm_nu > 1e-16
    
            % Solve the thrid linear system.
            b3 = [-vgradpen; al_lam3-dir2_norm_nu*alpha];
            b3 = Q*b3;
            dir_multipliers3 = linearsolve(L, U, b3, dimsystem);
            dir3 = dir_multipliers3(1:dim);
            multipliers3 = dir_multipliers3(dim+1:dimsystem);
    
            % Computing the main line serach direction.
            vgptd2 = vgradpen'*dir2; vgptd3 = vgradpen'*dir3;
            if vgptd3 <= tau*vgptd2
                theta = 1;
            else
                theta = (1-tau)*vgptd2/(vgptd2-vgptd3);
            end
            main_e_dir = (1-theta)*dir2 + theta*dir3;
            multipliers = (1-theta)*multipliers2 + theta*multipliers3;
         else
            main_e_dir = dir2;
            multipliers = multipliers2;
        end

        %main_e_dir
    
        %% Converting the main direction to the tangent space at current iteration point by Algorithm 7 in Huang et al (2017).
        main_dir_norm = norm(main_e_dir,2);
        Sdot_main = reshape(main_e_dir(1:k^2), [k,k]); 
        K_main = reshape(main_e_dir(k^2+1:mk), [m-k,k]);
        W_mian = reshape(main_e_dir(mk+1:end), [n-k,k]);

        main_dir.S = Sdot_main;
        % if m <= n
        %     M1_main = fun_beta([Sdot_main; K_main], vU, sU, m, k);
        %     M2_mian = fun_beta([zeros(k); W_mian], vV, sV, n, k);
        %     main_dir.Up = M1_main - x.U*Sdot_main;
        %     main_dir.Vp = M2_mian;
        % else
        %     M1_main = fun_beta([zeros(k); K_main], vU, sU, m, k);
        %     M2_mian = fun_beta([Sdot_main'; W_mian], vV, sV, n, k);
        %     main_dir.Up = M1_main;
        %     main_dir.Vp = M2_mian - x.V*Sdot_main';
        % end

        M1_main = fun_beta([zeros(k); K_main], vU, sU, m, k);
        M2_mian = fun_beta([zeros(k); W_mian], vV, sV, n, k);
        main_dir.Up = M1_main;
        main_dir.Vp = M2_mian;
        
        % Sdot_main
        % M1_main

        %VV = x.U*main_dir.S*x.V'+main_dir.Up*x.V'+x.U*(main_dir.Vp)'
        
        % G1 = x.U*main_dir.S*x.V' + x.U*main_dir.Vp' + main_dir.Up*x.V';
        % G1 = M.egrad2rgrad(x, G1)
        %
        % G1 = zeros(m,n);
        % for l = 1:dim
        %     G1 = G1 + main_e_dir(l)*base{l};
        % end
        % main_dir.S = x.U'*G1*x.V;
        % main_dir.Up = x.Upe*x.Upe'*G1*x.V;
        % main_dir.Vp = x.Vpe*x.Vpe'*G1'*x.U;
                
         %% Computing the correction direction if current point sufficiently close an KKT point.
        if KKT_res < cor_tol
            active_index = (cx+multipliers >= 0);
            active_Gradc = Gradc(:, active_index);
            al_dl_lam = abs(alpha./(sqrt(2)*delta.*multipliers)-1);
            omega = max(main_dir_norm^varrho,max(al_dl_lam(active_index))^kappa*main_dir_norm^2);
            
            % % Security guarantee for algorithm.
            % if omega == inf
            %     omega = main_dir_norm^varrho;
            % end
             
            % Compute newx_cor and c_i(newx_cor) for i \in active set
            newx_cor = retr(x, main_dir, k); % retraction
            newX_cor = newx_cor.U*(newx_cor.V'.*diag(newx_cor.S));
            
            active_newx_cor = zeros(nnz(active_index),1); 
            i_act = 1;
            for i = 1:constraint_number
                if i <= ineq_num && active_index(i)
                    active_newx_cor(i_act) = -newX_cor(ineqindices(i));
                    i_act = i_act + 1;
                elseif i > ineq_num && active_index(i)
                    active_newx_cor(i_act) = newX_cor(eqindices(i-ineq_num)) - A(eqindices(i-ineq_num));
                    i_act = i_act + 1;
                end
            end
            
            % quadprog, a matlab solver for QP.
            [coeff, ~, ~, ~, ~] = quadprog(HL, [], [], [], active_Gradc',...
                -active_newx_cor-omega, [], [], [], optimset('Display','off'));
            
            % Computing the correction direction on tanget space
            if norm(coeff,2) <= main_dir_norm
                Sdot_cor = reshape(coeff(1:k^2), [k,k]); 
                K_cor = reshape(coeff(k^2+1:mk), [m-k,k]);
                W_cor = reshape(coeff(mk+1:end), [n-k,k]);
        
                cor_dir.S = Sdot_cor;
                % if m <= n
                %     M1_cor = fun_beta([Sdot_cor; K_cor], vU, sU, m, k);
                %     M2_cor = fun_beta([zeros(k); W_cor], vV, sV, n, k);
                %     cor_dir.Up = M1_cor - x.U*Sdot_cor;
                %     cor_dir.Vp = M2_cor;
                % else
                %     M1_cor = fun_beta([zeros(k); K_cor], vU, sU, m, k);
                %     M2_cor = fun_beta([Sdot_cor'; W_cor], vV, sV, n, k);
                %     cor_dir.Up = M1_cor;
                %     cor_dir.Vp = M2_cor - x.V*Sdot_cor';
                % end

                M1_cor = fun_beta([zeros(k); K_cor], vU, sU, m, k);
                M2_cor = fun_beta([zeros(k); W_cor], vV, sV, n, k);
                cor_dir.Up = M1_cor;
                cor_dir.Vp = M2_cor;
            else
                cor_dir.S = 0; cor_dir.Up = 0; cor_dir.Vp = 0;
            end
        end
    
        %% curvilinear search
        t = 1; ls = 1; 
        xpencost = xcost - rho*sum(cx(ineq_num+1:end));
        %rho_cnewx = 0; 
        descent = sigma*t*(main_e_dir'*vgradpen);
        
        while ls <= ls_max_steps
            if KKT_res < cor_tol
                Dir.S = t*main_dir.S + t^2*cor_dir.S;
                Dir.Up = t*main_dir.Up + t^2*cor_dir.Up;
                Dir.Vp = t*main_dir.Vp + t^2*cor_dir.Vp;
            else
                Dir.S = t*main_dir.S;
                Dir.Up = t*main_dir.Up;
                Dir.Vp = t*main_dir.Vp;
                %Dir.M = Dir.S;
            end
            newx = retr(x, Dir, k); % retractionn
            %newx = M.retr(x, Dir, 1);
            newX = newx.U*(newx.V'.*diag(newx.S));

            PnewXA = P.*(newX-A);
            newxcost = 0.5*norm(PnewXA, 'fro')^2;
            
            % For constriants
            for i = 1:constraint_number
                if i <= ineq_num
                    cnewx = -newX(ineqindices(i));
                else
                    cnewx = newX(eqindices(i-ineq_num)) - A(eqindices(i-ineq_num));
                end
                if cnewx >= 0
                    ls_tag = 1;
                    break
                else
                    ls_cx(i) = cnewx;
                end
            end
            if ls_tag == 1
                t = t*varsigma;
                descent = varsigma*descent;
                ls_tag = 0;
                ls = ls +1;
                continue
            end

            % For penalty objective
            newxpencost = newxcost - rho*sum(ls_cx(ineq_num+1:end));
            if newxpencost - xpencost > descent
                t = t*varsigma;
                descent = varsigma*descent;
                ls_tag = 0;
                ls = ls +1;
                continue 
            end

            if ls_tag == 0
                %ls
                break
            end

        end

        % Update variables to new iterate
        x = newx; X = newX; cx = ls_cx;
        dS = diag(x.S); x.invS = 1./dS;
        PXA = PnewXA; xcost = newxcost;

        % Update gammas, mus and lambdas
        gammas = min(max(multipliers1,main_dir_norm),gamma_max);
        mus = m_ineq1; lambdas = m_eq1 - rho; 

    end
    xfinal = x; costfinal = xcost; residual = KKT_res;
    
      
    %% Support functions
    function egradLag = gradLagrangian(gf, mus, lambdas, ineqindices, eqindices, ineq_num, eq_num, m, n)
        Lams = zeros(m,n); Mus = zeros(m,n);
        for numeq = 1:eq_num
            keq = eqindices(numeq);
            Lams(keq) = lambdas(numeq);
        end
        for numineq = 1:ineq_num
            kineq = ineqindices(numineq);
            Mus(kineq) = -mus(numineq);
        end
        egradLag = gf + Lams + Mus;
    end
    
    % Calculating the KKT residual at (X, mus, Xgradnorm, eqindices, ineqindices)
    function val = KKT_residual(X, A, mus, Xgradnorm, eqindices, ineqindices, ineq_num, eq_num)
%        KKT_residual(X, mus, gradLagnorm, eqindices, ineqindices, ineq_num, eq_num);
        val = Xgradnorm^2;        
        compowvio = complementaryPowerViolation(X, mus, ineqindices, ineq_num);
        muspowvio = musposiPowerViolation(mus, ineq_num);
        
        val = val + compowvio;
        val = val + muspowvio;
   
        for numineq = 1: ineq_num
            index_ineq = ineqindices(numineq);
            cost_at_x = -X(index_ineq);
            violation = max(0, cost_at_x);
            val = val + violation^2;
        end
        
        for numeq = 1: eq_num
            index_eq = eqindices(numeq);
            cost_at_x = abs(X(index_eq)-A(index_eq));
            val = val + cost_at_x^2;
        end
                
        val = sqrt(val);   
    end
    
    function compowvio = complementaryPowerViolation(X, mus, ineqindices, ineq_num)
        compowvio = 0; 
        for numineq = 1: ineq_num
            index_ineq = ineqindices(numineq);
            violation = -mus(numineq)*X(index_ineq);
            compowvio = compowvio + violation^2;
        end
    end
    
    function musvio = musposiPowerViolation(mus, ineq_num)
        musvio = 0;
        for numineq = 1: ineq_num
            violation = max(-mus(numineq), 0);
            musvio = musvio + violation^2;
        end
    end

    function Y = retr(x, Z, p)

        [Qu, Ru] = qr([x.U, Z.Up], 0);
        [Qv, Rv] = qr([x.V, Z.Vp], 0);

        [U1, S1, V1] = svd(Ru*[x.S + Z.S, eye(p); eye(p), zeros(p)]*Rv');
    
        Y.U = Qu*U1(:, 1:p); 
        Y.V = Qv*V1(:, 1:p); 
        Y.S = S1(1:p, 1:p);
    end
    
    function val = linearsolve(L, U, b, n)  
        y = ltri(L, b, n);
        val = utri(U, y, n);
        function y = ltri(L, b, n)
            y = zeros(n,1);
            for j = 1:n-1
                y(j) = b(j)/L(j,j);
                b(j+1:n) = b(j+1:n) - y(j)*L(j+1:n,j);
            end
            y(n) = b(n)/L(n,n);
        end
        function x = utri(U, y, n)
            x=zeros(n,1);
            for j = n:-1:2
                x(j) = y(j)/U(j,j);
                y(1:j-1) = y(1:j-1) - x(j)*U(1:j-1,j);
            end
            x(1) = y(1)/U(1,1);
      end
    end
    
end