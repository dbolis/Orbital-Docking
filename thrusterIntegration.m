function [pulses,out] = thrusterIntegration(X0,options,tsim_Actu,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mt,mu,beta,p,q,eta,q_d,rho_d,drho_d)


global logicMat2 Tcontrol Fcontrol 

delt=1 % delta t integration time step
tstep=1 % continuous time step length
thrustMag=2 % constant magnitude thrust value
t0=0 % simulation begin time

xout=[transpose(X0),t0,0,0,0,0,0,0] % initialize output
pulseTableOut=[] % initialize pulse table output

for z=1:1:tsim_Actu 
    %find continuous control thrust
    logicMat2=[]
    Tcontrol=[]
    Fcontrol=[]
    [t,x] = ode45(@integrationLi,t0:delt:t0+tstep,X0,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mt,mu,beta,p,q,eta,q_d,rho_d,drho_d)
  
% PWM
t=[1:1:size(logicMat2,2)]
% thrusterArea=trapz(transpose(logicMat2(:,1:length(t))))*1/(length(t)-1) 

thrusterAreaFirst=trapz(transpose(logicMat2(:,1:floor(length(t)/2))))*((floor(length(t)/2)-1)/(length(t)-1))/(floor(length(t)/2)-1) % first half of 1 second step
thrusterAreaSecond=trapz(transpose(logicMat2(:,floor(length(t)/2):length(t))))*((length(t)-floor(length(t)/2))/(length(t)-1))/((length(t)-floor(length(t)/2))) % second half of 1 second step


eqTimeFirst=thrusterAreaFirst/thrustMag % equivalent time first half 
eqTimeSecond=thrusterAreaSecond/thrustMag % equivalent time second half 


dataFirst=round(eqTimeFirst,1)*10 %round first half to number of pulses 
dataSecond=round(eqTimeSecond,1)*10 %round second half to number of pulses
filterFirst=dataFirst>5 % find values larger than 5 pulses
filterSecond=dataSecond>5 % find values larger than 5 pulses

indecesFirst=find(filterFirst) % find indeces of larger than 5
indecesSecond=find(filterSecond) % find indeces of larger than 5
for j=1:length(indecesFirst)
    dataFirst(indecesFirst(j))=5 % replace indeces larger than 5 with 5
end
for j=1:length(indecesSecond)
    dataSecond(indecesSecond(j))=5 % replace indeces larger than 5 with 5
end

%create pulse table for 12 thrusters
pulseTable=[zeros(1,5-dataFirst(1)),thrustMag*ones(1,dataFirst(1)),thrustMag*ones(1,dataSecond(1)),zeros(1,5-dataSecond(1));
        zeros(1,5-dataFirst(2)),thrustMag*ones(1,dataFirst(2)),thrustMag*ones(1,dataSecond(2)),zeros(1,5-dataSecond(2));
        zeros(1,5-dataFirst(3)),thrustMag*ones(1,dataFirst(3)),thrustMag*ones(1,dataSecond(3)),zeros(1,5-dataSecond(3));
        zeros(1,5-dataFirst(4)),thrustMag*ones(1,dataFirst(4)),thrustMag*ones(1,dataSecond(4)),zeros(1,5-dataSecond(4));
        zeros(1,5-dataFirst(5)),thrustMag*ones(1,dataFirst(5)),thrustMag*ones(1,dataSecond(5)),zeros(1,5-dataSecond(5));
        zeros(1,5-dataFirst(6)),thrustMag*ones(1,dataFirst(6)),thrustMag*ones(1,dataSecond(6)),zeros(1,5-dataSecond(6));
        zeros(1,5-dataFirst(7)),thrustMag*ones(1,dataFirst(7)),thrustMag*ones(1,dataSecond(7)),zeros(1,5-dataSecond(7));
        zeros(1,5-dataFirst(8)),thrustMag*ones(1,dataFirst(8)),thrustMag*ones(1,dataSecond(8)),zeros(1,5-dataSecond(8));
        zeros(1,5-dataFirst(9)),thrustMag*ones(1,dataFirst(9)),thrustMag*ones(1,dataSecond(9)),zeros(1,5-dataSecond(9));
        zeros(1,5-dataFirst(10)),thrustMag*ones(1,dataFirst(10)),thrustMag*ones(1,dataSecond(10)),zeros(1,5-dataSecond(10));
        zeros(1,5-dataFirst(11)),thrustMag*ones(1,dataFirst(11)),thrustMag*ones(1,dataSecond(11)),zeros(1,5-dataSecond(11));
        zeros(1,5-dataFirst(12)),thrustMag*ones(1,dataFirst(12)),thrustMag*ones(1,dataSecond(12)),zeros(1,5-dataSecond(12))]
        

    
pulseAndTime=[pulseTable;[t0:0.1:t0+0.9]] % output pulse table
pulseTableOut=[pulseTableOut, pulseAndTime];



X0forLoop=X0 % updated initial conditions

% continuous integration over each of 10 constant thrust pulses
for k=0:1:9
    
    [tforLoop,xforLoop] = ode45(@pulseIntegration,(t0+k/10):0.1:(t0+(k+1)/10),X0forLoop,options,pulseTable(1:12,k+1),e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mt,mu,beta,p,q,eta,q_d,rho_d,drho_d)
    X0forLoop=transpose(xforLoop(end,:))

end

t0=t0+1 % increment outerloop time step
X0=X0forLoop % update initial conditions
xout=[xout; transpose(X0), t0, mean(transpose(Fcontrol)), mean(transpose(Tcontrol))]; % appent output, take mean of control force and contrl torque over each outer loop timestep
end

out=xout
pulses=pulseTableOut
