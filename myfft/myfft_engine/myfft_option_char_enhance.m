function val=myfft_option_char_enhance(model,inputVar, x) 
%{
[a, b] =xlsread(strcat(input_path,input_fileName),model);
[l,w]=size(a);
if scenario > l
    error('invalid scenario %d',  scenario);
end
inputVar = containers.Map(b,a(scenario,:));
NInter = inputVar('NInter');
Delta = inputVar('TT')/NInter;
%DeltaT = TT/NInter;
l = inputVar('l');
NCoef = 2*inputVar('NCoeff')-1;
%}

NInter = inputVar('NInter');
Delta = inputVar('TT')/NInter;
Rate = inputVar('r');
if strcmp(model, 'Gauss') %BS model
    Sigma = inputVar('sigma');
    val=exp(complex(0,1)*x*(Rate-Sigma^2/2)*Delta-Sigma^2*x.^2*Delta/2);
elseif  strcmp(model, 'Merton') %Merton's model
    Sigma = inputVar('sigma');
    Mlambda = inputVar('lambda');
    Mdelta = inputVar('delta');  
    Mq = inputVar('q');
    Mmu = inputVar('mu');
    Mgamma=Rate-Mq-0.5*Sigma^2-Mlambda*(exp(0.5*Mdelta^2+Mmu)-1);
    val=exp(complex(0,1)*Mgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+Mlambda*Delta*(exp(-0.5*x.^2*Mdelta^2+complex(0,1)*x*Mmu)-1));
elseif strcmp(model, 'Kou') %Kou's model
    Sigma = inputVar('sigma');
    Klamm = inputVar('lamm');
    Klamp = inputVar('lamp');
    Klam = inputVar('lambda');
    Kp = inputVar('p');
    Kgamma=Rate-0.5*Sigma^2-Klam*(Kp/(Klamp-1)-(1-Kp)/(Klamm+1));
    val=exp(complex(0,1)*Kgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+complex(0,1)*x*Klam*Delta.*(Kp./(Klamp-complex(0,1)*x)-(1-Kp)./(Klamm+complex(0,1)*x)));
    %Kgamma=Rate-0.5*Sigma^2-Klam*(Kp/(Klamp-1)-(1-Kp)/(Klamm+1));
    %val=exp(complex(0,1)*Kgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+complex(0,1)*x*Klam*Delta.*(Kp./(Klamp-complex(0,1)*x)-(1-Kp)./(Klamm+complex(0,1)*x)));
    
    
    %error('not implementted yet, do not know char function of Kou model');
elseif strcmp(model, 'NIG') %NIG
    NIGalpha = inputVar('alphaNIG');
    NIGbeta = inputVar('beta');
    NIGdelta = inputVar('delta');    
    NIGmu = inputVar('mu');
    w=-NIGdelta*( (NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+1)^2)^0.5 )-NIGmu;
    val=exp(complex(0,1)*x*Delta*(Rate+w+NIGmu)+NIGdelta*Delta*((NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+complex(0,1)*x).^2).^0.5));
    %temp test for density
    %val=exp(NIGdelta*((NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+complex(0,1)*x)^2)^0.5)+complex(0,1)*NIGmu*x);
elseif strcmp(model, 'VG') %VG model
    VGtheta = inputVar('theta');
    Sigma = inputVar('sigma');
    VGnu = inputVar('nu');
    w=log(1-VGtheta*VGnu-0.5*Sigma^2*VGnu)/VGnu;
    val=exp(complex(0,1)*Delta*x*(Rate+w) - (Delta/VGnu) * log(1-complex(0,1)*x*VGnu*VGtheta+0.5*Sigma^2*VGnu*x.^2) );
    %k=0.25;
    %sigma=0.3;
    %theta=0.03;
    %phi=-1/k * log(1+x^2*sigma^2*k/2-complex(0,1)*theta*k*x);
    %val=exp(phi);
elseif strcmp(model, 'TS') %TS model
    TSalphap = inputVar('alphap');
    TSalpham = inputVar('alpham');
    TScp = inputVar('Cp');
    TScm = inputVar('Cm');
    TSlambdap = inputVar('Lamp');
    TSlambdam = inputVar('Lamm');
    w=-gamma(-TSalphap)*TScp*((TSlambdap-1)^TSalphap-TSlambdap^TSalphap) - gamma(-TSalpham)*TScm*( (TSlambdam+1)^TSalpham-TSlambdam^TSalpham);
    val=exp( complex(0,1)*x*Delta*(Rate+w)+Delta*gamma(-TSalphap)*TScp*( (TSlambdap-complex(0,1)*x).^TSalphap-TSlambdap^TSalphap ) + Delta*gamma(-TSalpham)*TScm*(  (TSlambdam+complex(0,1)*x).^TSalpham - TSlambdam^TSalpham) );
else
    error('wrong model %s',model)
%elseif modelType==7
%    val=exp(-IGSdelta* ( (-2*complex(0,1)*x+IGSgamma^2).^0.5-IGSgamma));
end