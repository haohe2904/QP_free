%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% alpha(Q,A) that Algorithm 4 in Huang et al (2017), where Q is Orthogonal matrix.
%% i.e., to compute [Q Q_perp]^T*A
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = fun_alpha(A0,v,s,n,p)

A = A0;
for i = 1:p
    if i == 1
        vvtA = 2*v{i}*(v{i}'*A);
        A = A - vvtA;
    else
        A1 = A(1:i-1,:);
        A2 = A(i:n,:);
        vvtA2 = 2*v{i}*(v{i}'*A2);
        A2 = A2 - vvtA2;
        A = [A1;A2];
    end
end

%% Compute A <-- diag(s1,s2,...,sp, I_{n−p})A
for i = 1:p
    A(i,:) = s(i)*A(i,:);
end
end