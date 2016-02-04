function sol=Fitfkt_T2(param,c_0,data)

global bestparam besterr

c1=0.1;
c2=10;
for i=1:length(param)
    bm{1,i}=param(i)*c1;
    bm{2,i}=param(i)*c2;
end

lb=[bm{1,1} bm{1,2} bm{1,3} bm{1,4} bm{1,5} bm{1,6} bm{1,7} bm{1,8} bm{1,9} bm{1,10} bm{1,11} bm{1,12} bm{1,13} bm{1,14} bm{1,15} bm{1,16}];
ub=[bm{2,1} bm{2,2} bm{2,3} bm{2,4} bm{2,5} bm{2,6} bm{2,7} bm{2,8} bm{2,9} bm{2,10} bm{2,11} bm{2,12} bm{2,13} bm{2,14} bm{2,15} bm{2,16}];

options=optimset('MaxFunEvals',inf,'MaxIter',500);
[nparam,fval,exitflag,output]=fmincon(@(param)ode23s_solver_T2(param,c_0,data),bestparam,[],[],[],[],lb,ub,[],options);
output
