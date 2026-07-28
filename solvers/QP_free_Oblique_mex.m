function [xfinal, costfinal, residual, stopreason] = QP_free_Oblique_mex(d, s, A, C, X0, options)

    if isfield(options,'tolKKTres');             tol                  = options.tolKKTres;              end
    if isfield(options,'verbosity');             printyes             = options.verbosity;              end
    if isfield(options,'maxiter');               maxiter              = options.maxiter;                end
    if isfield(options,'mineigval_correction');  mineigval_correction = options.mineigval_correction;   end
    
    %% Initialization
    Time = tic;
    stopreason = 0;
    X = X0; I_d = eye(d);

   % M = obliquefactory(d,s);
    
    % Algoritmh papameters setting
    nu = 2.3; tau = 0.75; varrho = 2.4; kappa = 0.55;
    sigma = 0.45; varsigma = 0.5; gamma_max = 50;
    cor_tol = 1e-5; ls_max_steps = 5;
    
    % Problem informations
    dim = (d-1)*s; ineq_num = d*s; d_1 = d - 1;
    constraint_number = ineq_num;
    dimsystem = dim + constraint_number;
    Upsilon = 0.5;
    
    % Initial Lagrange multiplier vectors and gammas
    mus = zeros(ineq_num, 1); 
    gammas = 1e-1*ones(constraint_number, 1);
    %gammas = mus;
    
    % Store the vector expression of Riemannian gradients of constriant
    Ineqgrads = zeros(dim,ineq_num); 
    zerosdim = zeros(dim,1);
    
    % Store cx which is the constriant function values
    cx = zeros(constraint_number,1);
    ls_cx = zeros(constraint_number,1);
    
    %% Make cx
    for i = 1:constraint_number
      cx(i) = -X(i);
    end
    
    % Compute objective value
    XtA = X'*A; sX = sum(X,2);
    Xcost = -norm(XtA, 'fro')^2 + Upsilon*(norm(sX,2)^2/s-1);
    
    % Update line search parameter flag
    ls_tag = 0;

    % Store reflector
    qr_V = cell(s,1); tau_V = cell(s,1);
    zerosV = zeros(d,d-1);
    
    %% Main loop
    for iter = 1:maxiter
        
        cx_gam = (cx.^2+gammas.^2).^0.5;
        alpha = cx./cx_gam + 1;
        beta = (1-gammas./cx_gam).^0.5;
        % delta = -beta./cx;
        delta = (1./(cx_gam.*(gammas+cx_gam))).^.5;
   
        % Compute v,s at Q by Algorithm 3 in Huang et al (2017), where Q is Orthogonal matrix.
        for i = 1:s
            [qr_V{i}, tau_V{i}] = dgeqrf_mex([X(:,i), zerosV]);
        end
    
        %% Make the matrix representation for the Riemannian gradients of constriants
        for i = 1:ineq_num
            gradc0 = zeros(d,1);
            [row, col] = ind2sub([d,s], i); 
            gradc0(row) = -1;

            % Compute the vector representation of gradc by Algorithm 4 in Huang et al (2017).
            vgradc0 = apply_q_mex(gradc0,qr_V{col}, tau_V{col}, 'T'); 
            vgradc0(1) = [];
            vgradc = zerosdim;
            vgradc(1+(col-1)*d_1:col*d_1) = vgradc0;
            
            % Make ineqgrads
            Ineqgrads(:,i) = vgradc;
        end
       
        Gradc = Ineqgrads;
        aGradct = alpha.*(Gradc');
        Dbeta = diag(-sqrt(2).*beta);
        
        %% Compute Euclidean and Riemannian gradient (in vector representation) of Lagrangian function
        CX = -A*XtA';
        egradf = 2*CX + (2*Upsilon/s)*(sX*ones(1,s));
        vgradm = zeros(d,s);
        for i = 1:s
            vgradm(:,i) = apply_q_mex(egradf(:,i),qr_V{i}, tau_V{i}, 'T');
        end
        vgradm(1,:) = [];
        vgradf = vgradm(:);

        %%
        Mus = reshape(mus, [d,s]);
        mMus = zeros(d,s);
        for i = 1:s
            mMus(:,i) = apply_q_mex(Mus(:,i),qr_V{i}, tau_V{i}, 'T');
        end
        mMus(1,:) = [];
        vMus = mMus(:);

        vgradLag = vgradf - vMus;
        gradLagnorm = norm(vgradLag, 2);
       % gradLagnorm=0;
    
        % Compute KKT residual and check the stopping conditions
        KKT_res = KKT_residual(X, mus, gradLagnorm, ineq_num);
        ttime = toc(Time);
        if printyes >= 3
            if iter == 1
                fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, Xcost, KKT_res, ttime);
                kkt_flag0 = floor(log10(KKT_res));
            else
                kkt_flag1 = floor(log10(KKT_res));
                if kkt_flag1 < kkt_flag0
                    fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, Xcost, KKT_res, ttime);
                    kkt_flag0 = kkt_flag1;
                end
            end
        elseif printyes >= 2
            fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, Xcost, KKT_res, ttime);
        elseif printyes >= 1
            if mod(iter, 100) == 0 && iter ~= 0
                fprintf('Iter: %d, Cost: %f, KKT residual: %.4e, Time: %.4e \n', iter, Xcost, KKT_res, ttime);
            end
        end

        if KKT_res <= tol
            fprintf('KKT Residual tolerance reached\n');
            stopreason = 1;
            break;
        end
    
        %% Compute Riemannian Hessian (in matrix representation) of Lagrangian function
        MH = cell(s^2,1);
        %sX = sum(X,2);
        for i = 1:s
            egradL = 2*CX(:,i) - Mus(:,i) + (2*Upsilon/s)*sX;
            for k = i:s
                if k == i
                    H0 = 2*C + (2*Upsilon/s - X(:,i)'*egradL)*I_d;
                    BtH = apply_q_mex(H0, qr_V{i}, tau_V{i}, 'T');
                    BtH(1,:) = [];
                    HtB = BtH';
                    BtHtB =  apply_q_mex(HtB, qr_V{i}, tau_V{i}, 'T');
                    BtHtB(1,:) = [];
                    MH{(i-1)*s + k} = BtHtB;
                else
                    H0 = 2*Upsilon/s*I_d;
                    BtH = apply_q_mex(H0, qr_V{k}, tau_V{k}, 'T'); 
                    BtH(1,:) = [];
                    HtB = BtH';
                    BtHtB =  apply_q_mex(HtB, qr_V{i}, tau_V{i}, 'T');
                    BtHtB(1,:) = [];
                    MH{(i-1)*s + k} = BtHtB';
                    MH{(k-1)*s + i} = BtHtB;
                end
            end
        end

        MH = reshape(MH, s, s);
        HL = cell2mat(MH);

        % Ensure that matrix BtHtB is positive definite.
        [V,T] = schur(HL);
        for k = 1 : dim
            if T(k,k) < 1e-5  
                T(k,k) = mineigval_correction;
            end
        end
        HL = V*T*V'; HL = 0.5*(HL.'+HL);
        %HL = eye(dim);
    
        % Make coefficient matrix B
        B = [HL Gradc; aGradct Dbeta];
    
        %% Compute the direction and Lagrange multipliers by solving three linear system with LU decomposition.
        [L, U, Q] = lu(B); 

        % Solve the first linear system.
        b1 = [-vgradf; zeros(constraint_number,1)];
        %dir_multipliers0 = B\b1;
        b1 = Q*b1;
        dir_multipliers1 = linearsolve(L, U, b1, dimsystem); 
        %dir1 = dir_multipliers1(1:dim);
        m_ineq1 = dir_multipliers1(dim+1:dimsystem);
        
        % Solve the second linear system.
        al_lam3 = alpha.*(min(m_ineq1,0).^3);
        b2 = [-vgradf; al_lam3];
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
            b3 = [-vgradf; al_lam3-dir2_norm_nu*alpha];
            b3 = Q*b3;
            dir_multipliers3 = linearsolve(L, U, b3, dimsystem);
            dir3 = dir_multipliers3(1:dim);
            multipliers3 = dir_multipliers3(dim+1:dimsystem);
    
            % Computing the main line serach direction.
            vgptd2 = vgradf'*dir2; vgptd3 = vgradf'*dir3;
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

        %main_e_dir = dir1;
    
        %% Converting the main direction to the tangent space at current iteration point by Algorithm 7 in Huang et al (2017).
        main_dir_norm = norm(main_e_dir,2);
        main_e_dirm = reshape(main_e_dir, [d-1,s]);
        main_dir = zeros(d,s);
        for i = 1:s
            main_dir(:,i) = apply_q_mex([0; main_e_dirm(:,i)], qr_V{i}, tau_V{i}, 'N');
        end
                
         %% Computing the correction direction if current point sufficiently close an KKT point.
        if KKT_res < cor_tol
            active_index = (cx+multipliers >= 0);
            active_Gradc = Gradc(:, active_index);
            al_dl_lam = abs(alpha./(sqrt(2)*delta.*multipliers)-1);

            num_act = nnz(active_index);
            if num_act == 0
                omega = main_dir_norm^varrho;
            else
                omega = max(main_dir_norm^varrho,max(al_dl_lam(active_index))^kappa*main_dir_norm^2);         
            end
            % Compute newX_cor and c_i(newx_cor) for i \in active set
            newX_cor = normalize(X+main_dir, 'norm');
            
            active_newx_cor = zeros(num_act,1); 
            i_act = 1;
            for i = 1:constraint_number
                if active_index(i)
                    active_newx_cor(i_act) = -newX_cor(i);
                    i_act = i_act + 1;
                end
            end
            % quadprog, a matlab solver for QP.
            [coeff, ~, ~, ~, ~] = quadprog(HL, [], [], [], active_Gradc',...
                -active_newx_cor-omega, [], [], [], optimset('Display','off'));

            %[coeff, ~, ~, ~, ~] = quadprog(HL, [], [], [], [],...
             %   [], [], [], [], optimset('Display','off'));
            
            % Computing the correction direction on tanget space
            if norm(coeff,2) <= main_dir_norm
                cor_e_dir = reshape(coeff, [d-1,s]);
                cor_dir = zeros(d,s);
                for i = 1:s
                    cor_dir(:,i) = apply_q_mex([0; cor_e_dir(:,i)], qr_V{i}, tau_V{i}, 'N');
                end
            else
                cor_dir = 0;
            end
        end
    
        %% curvilinear search
        t = 1; ls = 1; 
        descent = sigma*t*(main_e_dir'*vgradf);
        
        while ls <= ls_max_steps
            if KKT_res < cor_tol
                Dir = t*main_dir + t^2*cor_dir;
            else
                Dir = t*main_dir;
                %Dir.M = Dir.S;
            end
            newX = normalize(X+Dir, 'norm'); % retractionn
            

            newXtA = newX'*A; newsX = sum(newX,2);
            newXcost = -norm(newXtA, 'fro')^2 + Upsilon*(norm(newsX,2)^2/s-1);
            
            % For constriants
            for i = 1:constraint_number
                cnewx = -newX(i);
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

            % For objective
            if newXcost - Xcost > descent
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
        X = newX; sX = newsX; 
        XtA = newXtA; Xcost = newXcost;

        if ls <= ls_max_steps
            cx = ls_cx;
        else
            for i = 1:constraint_number
                cx(i) = -X(i);
            end
        end

        % Update gammas, mus and lambdas
        gammas = min(max(m_ineq1,main_dir_norm),gamma_max);
        mus = m_ineq1;  

    end
    xfinal = X; costfinal = Xcost; residual = KKT_res;
    
      
    %% Support functions

    % Calculating the KKT residual at (X, mus)
    function val = KKT_residual(X, mus, Xgradnorm, ineq_num)
        val = Xgradnorm^2;        
        compowvio = complementaryPowerViolation(X, mus, ineq_num);
        muspowvio = musposiPowerViolation(mus, ineq_num);
        
        val = val + compowvio;
        val = val + muspowvio;
   
        for numineq = 1: ineq_num
            cost_at_x = -X(numineq);
            violation = max(0, cost_at_x);
            val = val + violation^2;
        end    
        val = sqrt(val);   
    end
    
    function compowvio = complementaryPowerViolation(X, mus, ineq_num)
        compowvio = 0; 
        for numineq = 1: ineq_num
            violation = -mus(numineq)*X(numineq);
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