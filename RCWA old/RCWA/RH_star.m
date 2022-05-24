function S_AB = RH_star(S_A,S_B)
% Redheffer Star product
% Input is cell array of matrices.
sz=size(S_A{1,1});
D=S_A{1,2}/(eye(sz)-S_B{1,1}*S_A{2,2});
F=S_B{2,1}/(eye(sz)-S_A{2,2}*S_B{1,1});

S_AB{1,1}=S_A{1,1}+D*S_B{1,1}*S_A{2,1};
S_AB{1,2}=D*S_B{1,2};
S_AB{2,1}=F*S_A{2,1};
S_AB{2,2}=S_B{2,2}+F*S_A{2,2}*S_B{1,2};
end

