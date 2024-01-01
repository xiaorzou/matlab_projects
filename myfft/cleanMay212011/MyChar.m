function val=MyChar(modelType,x)
global Delta Sigma Rate  %BS

global Mq Mlambda Mmu Mgamma Mdelta %Merton's Model
global Kp Klam Klamp Klamm %Kou's model
global NIGmu NIGdelta NIGalpha NIGbeta %Normal Inverse Gaussian Model
global VGnu VGmu VGtheta %VG model
global TSgamma TSalpha TSalphap TSalpham TSlambdap TSlambdam TScp TScm %Tempered Stable Process

global IGSdelta IGSgamma %Inverse Gaussian subordinating 

if modelType==1 %BS model
    val=exp(complex(0,1)*x*(Rate-Sigma^2/2)*Delta-Sigma^2*x.^2*Delta/2);
elseif  modelType==2 %Merton's model
    val=exp(complex(0,1)*Mgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+Mlambda*Delta*(exp(-0.5*x.^2*Mdelta^2+complex(0,1)*x*Mmu)-1));
elseif modelType==3 %Kou's model
    Kgamma=Rate-0.5*Sigma^2-Klam*(Kp/(Klamp-1)-(1-Kp)/(Klamm+1));
    val=exp(complex(0,1)*Kgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+complex(0,1)*x*Klam*Delta.*(Kp./(Klamp-complex(0,1)*x)-(1-Kp)./(Klamm+complex(0,1)*x)));
    
    
    %Kgamma=Rate-0.5*Sigma^2-Klam*(Kp/(Klamp-1)-(1-Kp)/(Klamm+1));
    %val=exp(complex(0,1)*Kgamma*x*Delta-0.5*x.^2*Sigma^2*Delta+complex(0,1)*x*Klam*Delta.*(Kp./(Klamp-complex(0,1)*x)-(1-Kp)./(Klamm+complex(0,1)*x)));
    
    %error('not implementted yet, do not know char function of Kou model');
elseif modelType==4 %NIG
    w=-NIGdelta*( (NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+1)^2)^0.5 )-NIGmu;
    val=exp(complex(0,1)*x*Delta*(Rate+w+NIGmu)+NIGdelta*Delta*((NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+complex(0,1)*x).^2).^0.5));
    %temp test for density
    %val=exp(NIGdelta*((NIGalpha^2-NIGbeta^2)^0.5-(NIGalpha^2-(NIGbeta+complex(0,1)*x)^2)^0.5)+complex(0,1)*NIGmu*x);
elseif modelType==5 %VG model
    w=log(1-VGtheta*VGnu-0.5*Sigma^2*VGnu)/VGnu;
    val=exp(complex(0,1)*Delta*x*(Rate+w) - (Delta/VGnu) * log(1-complex(0,1)*x*VGnu*VGtheta+0.5*Sigma^2*VGnu*x.^2) );
    %k=0.25;
    %sigma=0.3;
    %theta=0.03;
    %phi=-1/k * log(1+x^2*sigma^2*k/2-complex(0,1)*theta*k*x);
    %val=exp(phi);
elseif modelType==6 %TS model
    w=-gamma(-TSalphap)*TScp*((TSlambdap-1)^TSalphap-TSlambdap^TSalphap) - gamma(-TSalpham)*TScm*( (TSlambdam+1)^TSalpham-TSlambdam^TSalpham);
    val=exp( complex(0,1)*x*Delta*(Rate+w)+Delta*gamma(-TSalphap)*TScp*( (TSlambdap-complex(0,1)*x).^TSalphap-TSlambdap^TSalphap ) + Delta*gamma(-TSalpham)*TScm*(  (TSlambdam+complex(0,1)*x).^TSalpham - TSlambdam^TSalpham) );
elseif modelType==7
    val=exp(-IGSdelta* ( (-2*complex(0,1)*x+IGSgamma^2).^0.5-IGSgamma));
end