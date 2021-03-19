function [pulses,out] = thrusterIntegration(X0,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)


global logicMat2 Tcontrol Fcontrol 

delt=1
tstep=1
thrustMag=2
t0=0

xout=[transpose(X0),t0,0,0,0,0,0,0]
pulseTableOut=[]

for z=1:1:20
    logicMat2=[]
    Tcontrol=[]
    Fcontrol=[]
    [t,x] = ode45(@integrationLi,t0:delt:t0+tstep,X0,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
  

t=[1:1:size(logicMat2,2)]
thrusterArea=trapz(transpose(logicMat2(:,1:length(t))))*1/(length(t)-1) %% check lengths t to make sure trapz integration is right (tried to fix already)

thrusterAreaFirst=trapz(transpose(logicMat2(:,1:floor(length(t)/2))))*((floor(length(t)/2)-1)/(length(t)-1))/(floor(length(t)/2)-1)
thrusterAreaSecond=trapz(transpose(logicMat2(:,floor(length(t)/2):length(t))))*((length(t)-floor(length(t)/2))/(length(t)-1))/((length(t)-floor(length(t)/2)))

test=thrusterAreaFirst+thrusterAreaSecond


eqTimeFirst=thrusterAreaFirst/thrustMag
eqTimeSecond=thrusterAreaSecond/thrustMag


dataFirst=round(eqTimeFirst,1)*10
dataSecond=round(eqTimeSecond,1)*10
filterFirst=dataFirst>5
filterSecond=dataSecond>5

indecesFirst=find(filterFirst)
indecesSecond=find(filterSecond)
for j=1:length(indecesFirst)
    dataFirst(indecesFirst(j))=5
end
for j=1:length(indecesSecond)
    dataSecond(indecesSecond(j))=5
end

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
        

    
pulseAndTime=[pulseTable;[t0:0.1:t0+0.9]]
pulseTableOut=[pulseTableOut, pulseAndTime];



X0forLoop=X0
for k=0:1:9
    
    [tforLoop,xforLoop] = ode45(@pulseIntegration,(t0+k/10):0.1:(t0+(k+1)/10),X0forLoop,options,pulseTable(1:12,k+1),e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
    X0forLoop=transpose(xforLoop(end,:))

end
t0=t0+1
X0=X0forLoop
xout=[xout; transpose(X0), t0, mean(transpose(Fcontrol)), mean(transpose(Tcontrol))];
end

out=xout
pulses=pulseTableOut
