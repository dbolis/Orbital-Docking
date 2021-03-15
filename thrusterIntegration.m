function [pulses,out] = thrusterIntegration(delt,X0,options,res,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)


global logicMat2 Tcontrol Fcontrol 


tstep=1
thrustMag=2
t0=0

xout=[transpose(X0),t0,0,0,0,0,0,0]
pulseTableOut=[]

for z=1:1:1500
    logicMat2=[]
    Tcontrol=[]
    Fcontrol=[]
    [t,x] = ode45(@integrationLi,t0:delt:t0+tstep,X0,options,e_v,p_v,i_v,omega_v,Re,J2,omega_E,AerS,c_D,e_c,n_c,e_t,n_t,Ic,It,mc,mu,beta,p,q,eta,q_d,rho_d,drho_d)
%     timeAll=size(t);
%     q_t=transpose(x(1:timeAll(1),5:8));
%     w_c=transpose(x(1:timeAll(1),9:11));
%     w_t=transpose(x(1:timeAll(1),2:4));
%     rho_c=transpose(x(1:timeAll(1),16:18));
%     drho_c=transpose(x(1:timeAll(1),19:21));
%     rho_t=transpose(x(1:timeAll(1),22:24));
%     drho_t=transpose(x(1:timeAll(1),25:27));
%     rho=rho_c-rho_t
%     drho=drho_c-drho_t
%     q_c=transpose(x(1:timeAll(1),12:15));
%     q_r=quatProdMat(quatRecipMat(q_t),q_c); 
%     logicMat2=[]
%     for i=1:1:size(t)
% 
%     % Relative position
%         Ch=0.5
%         rho_obs=rho_d
%         d_0=1.9
%         delta_0=2
%         alpha=0.5
%         etaPos=5e-4
% 
% 
%         if norm(rho(1:3,i)-rho_obs)>=delta_0
%             k=0
%         else
%             k=0
%         end
% 
%         gradUatt = Ch*(rho(1:3,i)-rho_d)/sqrt((norm(rho(1:3,i)-rho_d))^2+1)
% 
%         gradUrep = 4*k*(rho(1:3,i)-rho_obs)*((norm(rho(1:3,i)-rho_obs))^2-delta_0^2)/((norm(rho(1:3,i)-rho_obs))^2-d_0^2)^3
% 
%         sRel=drho(1:3,i)+alpha*(gradUatt+gradUrep)
%         sRelout(1:3,i)=sRel
% 
%         if abs(sRel(1))>=etaPos
%             s1Pos=tanh(sRel(1))
%         else
%             s1Pos=sRel(1)/etaPos
%         end
% 
%         if abs(sRel(2))>=etaPos
%             s2Pos=tanh(sRel(2))
%         else
%             s2Pos=sRel(2)/etaPos
%         end
% 
%         if abs(sRel(3))>=etaPos
%             s3Pos=tanh(sRel(3))
%         else
%             s3Pos=sRel(3)/etaPos
%         end
% 
%         satPos=[s1Pos; s2Pos; s3Pos]
% 
% 
%         grad2Uatt=Ch*(((norm(rho(1:3,i)-rho_d))^2+1)*eye(3)-(rho(1:3,i)-rho_d)*transpose(rho(1:3,i)-rho_d))/((norm(rho(1:3,i)-rho_d))^2+1)^(3/2)
% 
% 
% 
% 
% 
%         gamma1=0.02 ;
%         gamma2=0.005;
%         gamma3=0.05;
%         pdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
%         mdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
%         % spos=drho; % s position control, excluding APF term  Eq 20 Li
%         % satt=derror(2:4); % s att control, excluding APF term Eq 36 
%         Fclog=-mc*alpha*grad2Uatt*drho(1:3,i)-(gamma1+mc*gamma3)*satPos-pdiag*satPos;
%         Fc(1:3,i)=-mc*alpha*grad2Uatt*drho(1:3,i)-(gamma1+mc*gamma3)*satPos-pdiag*satPos;  %Eq 27 Li excluding APF terms
% 
%         Ar=(q_r(1,i)^2-transpose(q_r(2:4,i))*q_r(2:4,i))*eye(3)+2*(q_r(2:4,i))*transpose(q_r(2:4,i))-2*q_r(1,i)*[0, -q_r(4,i), q_r(3,i); %EQ 10 Liu
%                                                                                                     q_r(4,i), 0, -q_r(2,i);
%                                                                                                     -q_r(3,i), q_r(2,i), 0]; 
% 
%         w_r=w_c(1:3,i)-Ar*w_t(1:3,i);
%         w_rOut(1:3,i)=w_r;
% 
%         w_r_tilde = [0, -w_r(3), w_r(2); %EQ 12 Liu
%                      w_r(3), 0, -w_r(1);
%                       -w_r(2), w_r(1), 0];
% 
%         G_r=[0, -transpose(w_r);
%             w_r, -w_r_tilde];
% 
%         dq_r=0.5*G_r*q_r(1:4,i);
% 
% 
%         gamma1=0.02 ;
%         gamma2=0.005;
%         gamma3=0.05;
%         pdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
%         mdiag=[0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5];
% 
%         Chrel=0.5;
%         etaAtt=5e-4
%         beta=0.6;
%         gradUquat=Chrel*(q_r(2:4,i)-[0;0;0])/sqrt((norm(q_r(2:4,i)-[0;0;0]))^2+1);
% 
%         sAtt=dq_r(2:4)+beta*gradUquat
%         sAttout(1:3,i)=sAtt
% 
%         if abs(sAtt(1))>=etaAtt
%             s1Att=sign(sAtt(1))
%         else
%             s1Att=sAtt(1)/etaAtt
% 
%         end
% 
%         if abs(sAtt(2))>=etaAtt
%             s2Att=sign(sAtt(2))
%         else
%             s2Att=sAtt(2)/etaAtt
%         end
% 
%         if abs(sAtt(3))>=etaAtt
%             s3Att=sign(sAtt(3))
%         else
%             s3Att=sAtt(3)/etaAtt
%         end
% 
%         satAtt=[s1Att; s2Att; s3Att]
% 
%         grad2Uquat=Chrel*(((norm(q_r(2:4,i)-[0;0;0]))^2+1)*eye(3)-(q_r(2:4,i)-[0;0;0])*transpose(q_r(2:4,i)-[0;0;0]))/((norm(q_r(2:4,i)-[0;0;0]))^2+1)^(3/2);
% 
%         T=[0, -q_r(4,i), q_r(3,i); q_r(4,i), 0, -q_r(2,i);-q_r(3,i), q_r(2,i), 0] + q_r(1,i)*eye(3);
%         dT=[0, -dq_r(4), dq_r(3); dq_r(4), 0, -dq_r(2);-dq_r(3), dq_r(2), 0] + dq_r(1)*eye(3);
%         P=inv(T);
%         dP=inv(dT);
%         Istar=transpose(P)*Ic*P;
% 
%         IcPdq_r=Ic*P*dq_r(2:4);
%         Pdq_r=P*dq_r(2:4);
%         Arw_t=Ar*w_t(1:3,i);
%         twoPdq_r=2*P*dq_r(2:4);
%         w_t_tilde = [0, -w_t(3,i), w_t(2,i); %EQ 8 Liu
%                      w_t(3,i), 0, -w_t(1,i);
%                       -w_t(2,i), w_t(1,i), 0];
%         dw_t=It\(-w_t_tilde*It*w_t(1:3,i));
% 
%         Cstar=-Istar/dP*P-2*transpose(P)*[0, -IcPdq_r(3), IcPdq_r(2); IcPdq_r(3), 0, -IcPdq_r(1); -IcPdq_r(2), IcPdq_r(1), 0]*P;
%         Nstar=transpose(P)*([0, -Pdq_r(3), Pdq_r(2); Pdq_r(3), 0, -Pdq_r(1); -Pdq_r(2), Pdq_r(1), 0]*...
%             Ic*Ar*w_t(1:3,i))+transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1);-Arw_t(2), Arw_t(1), 0]*Ic*P*dq_r(2:4))+...
%             0.5*transpose(P)*([0, -Arw_t(3), Arw_t(2); Arw_t(3), 0, -Arw_t(1); -Arw_t(2), Arw_t(1), 0]*Ic*Ar*w_t(1:3,i))-...
%             0.5*transpose(P)*Ic*([0, -twoPdq_r(3), twoPdq_r(2); twoPdq_r(3), 0, -twoPdq_r(1);-twoPdq_r(2), twoPdq_r(1), 0]*Ar*w_t(1:3,i)-Ar*dw_t);
% 
%         Tc=Nstar-Cstar*beta*gradUquat-Istar*beta*grad2Uquat*dq_r(2:4)-gamma2*satAtt-mdiag*satAtt;
%         Tcprimelog = 2*transpose(T)*Tc
%         logicMat2=[logicMat2, thrusterLogic2(Fclog,Tcprimelog,q_c(1:4,i), x(i,1), i_v, omega_v)]
%         
%     end
  

t=[1:1:size(logicMat2,2)]
thrusterArea=trapz(transpose(logicMat2(:,1:length(t))))*1/(length(t)-1) %% check lengths t to make sure trapz integration is right (tried to fix already)
% 
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

% eqTime=thrusterArea/thrustMag
% eqTimeHalf=eqTime/2
% 
% 
% data=round(eqTimeHalf,1)*10
% filter=data>5
% indeces=find(filter)
% for j=1:length(indeces)
%     data(indeces(j))=5
% end
% pulseTable=[zeros(1,5-data(1)), thrustMag*ones(1,data(1)*2), zeros(1,5-data(1)); %% having equal time around center of pulse neccesitates rounding to even number of pulses... 
%       zeros(1,5-data(2)), thrustMag*ones(1,data(2)*2), zeros(1,5-data(2)); %% maybe can try with odd pulses if this doesntw work
%       zeros(1,5-data(3)), thrustMag*ones(1,data(3)*2), zeros(1,5-data(3));
%       zeros(1,5-data(4)), thrustMag*ones(1,data(4)*2), zeros(1,5-data(4));
%       zeros(1,5-data(5)), thrustMag*ones(1,data(5)*2), zeros(1,5-data(5));
%       zeros(1,5-data(6)), thrustMag*ones(1,data(6)*2), zeros(1,5-data(6));
%       zeros(1,5-data(7)), thrustMag*ones(1,data(7)*2), zeros(1,5-data(7));
%       zeros(1,5-data(8)), thrustMag*ones(1,data(8)*2), zeros(1,5-data(8));
%       zeros(1,5-data(9)), thrustMag*ones(1,data(9)*2), zeros(1,5-data(9));
%       zeros(1,5-data(10)), thrustMag*ones(1,data(10)*2), zeros(1,5-data(10));
%       zeros(1,5-data(11)), thrustMag*ones(1,data(11)*2), zeros(1,5-data(11));
%       zeros(1,5-data(12)), thrustMag*ones(1,data(12)*2), zeros(1,5-data(12))]

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
