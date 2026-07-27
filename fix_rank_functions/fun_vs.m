%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute v,s at Q by Algorithm 3 in Huang et al (2017), where Q is Orthogonal matrix.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v, s] = fun_vs(Q, n, p)
v = cell(p,1); s = zeros(p,1);
for i = 1:p
    k = n - i + 1;
    e1 = zeros(k,1); e1(1) = 1;
    z = Q(:,i);
    z = z(i:end);
    sgn = -sign(z(1));
    s(i) = sgn;
    a = sgn*norm(z,2);
    z_ae = z - a*e1;
    vi = z_ae/norm(z_ae,2);
    v{i} = vi;

    if i == 1
        Q = Q - 2*vi*(vi'*Q);
    else
        Q2 = Q(i:end,:);
        vivitQ2 = 2*vi*(vi'*Q2);
        Q = [Q(1:i-1,:); Q2-vivitQ2];
    end
end
end