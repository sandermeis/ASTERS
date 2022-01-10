function c = get_mode(c_ref,c_inc,Sg,W_0,V_0,W,V)
c_mn=(Sg{1,2})\(c_ref-Sg{1,1}*c_inc);
c_pl=Sg{2,1}*c_inc+Sg{2,2}*c_mn;
A=W\W_0+V\V_0;
B=W\W_0-V\V_0;
c{1}=0.5*A*c_pl+0.5*B*c_mn;
c{2}=0.5*B*c_pl+0.5*A*c_mn;
end

