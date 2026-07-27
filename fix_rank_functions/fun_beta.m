%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beta(Q,A) that Algorithm 5 in Huang et al (2017), where Q is Orthogonal matrix.
%% i.e., to compute [Q Q_perp]*A
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = fun_beta(A0,v,s,n,p)

A = A0;
%% Compute A <-- diag(s1,s2,...,sp, I_{n−p})A
for i = 1:p
    A(i,:) = s(i)*A(i,:);
end

for i = p:-1:1
    if i > 1
        A1 = A(1:i-1,:);
        A2 = A(i:n,:);
        vvtA2 = 2*v{i}*(v{i}'*A2);
        A2 = A2 - vvtA2;
        A = [A1;A2];
    else
        vvtA = 2*v{i}*(v{i}'*A);
        A = A - vvtA;
    end
end
end