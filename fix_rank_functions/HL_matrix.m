%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the matix representation HL by HL_A, HL_B and HL_C.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HL = HL_matrix(HL_A,HL_B,HL_C,dim,m,p)

p2 = p^2; mp = m*p;
%% Generate a zero matrix
HL = zeros(dim);
for i = 1:dim
    for j = i:dim
        if j == i
            if j <= p2
               HL(i,j) = HL_A{i}(j);
            elseif j > mp
                HL(i,j) = HL_B{i}(j-mp);
            else
                HL(i,j) = HL_C{i}(j-p2);
            end
        else
            if j <= p2
                HL(i,j) = HL_A{i}(j);
                HL(j,i) = HL(i,j);
            elseif j > mp
                HL(i,j) = HL_B{i}(j-mp);
                HL(j,i) = HL(i,j);
            else
                HL(i,j) = HL_C{i}(j-p2);
                HL(j,i) = HL(i,j);
            end
        end 
    end
end
end