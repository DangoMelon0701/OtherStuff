clear,clc

%TODAS LAS UNIDADES DE: DISTANCIA ESTAN EN cm, DENSIDAD EN g/cm3, COEF. ABSORCION EN cm2/g 

co=0; %corrimiento del eje
fprintf('Calculo de la eficiencia de una fuente puntual ubicada a una distancia %f respecto al eje central del detector en funcion de la energia y la distancia fuente detector\n\n',co);

E=0.1:0.01:1.2;   %0.03:0.02:1.5; %energia de la radiaci?n (MeV)
d=0:1:3; %1.4:0.5:10.4  distancia fuente-detector inicial 

%DETECTOR
r=2.54; %radio detector
l=5.08; %largo detector
xv=0.0508; %grosor ventana
denv=2.6984; %densidad ventana (aluminio)
den=3.67; %densidad NaI

ni=length(E);
j=length(d);
e=zeros(j,ni);
n=128; %numero de intervalos 64 bueno  96, 128 muy bueno
%COEFICIENTES ATENUACION
nai=fopen('coenai.txt');
tabnai=textscan(nai, '%f %f %f');
fclose(nai);
alum=fopen('coeal.txt');
tabal=textscan(alum, '%f %f %f');
fclose(alum);
car=fopen('coecar.txt');
tabcar=textscan(car, '%f %f %f');
fclose(car);

eang=1/(2*pi); %factor de angulo solido 1/2pi     exp(-uventana*denv*xv)
inicio=clock;
%INICIO DE CALCULO
            for k=1:j
                d1=d(k);
                for i=1:ni
                    Et=E(i);
                    e12=0;
                    uventana=u(Et,tabal,1); %coef atenuaci?n ventana(Al)
                    udetector=u(Et,tabnai,2); %coef absorci?n NaI
%                     eang=1/(2*pi); %factor de angulo solido 1/2pi    
                    if co>=0 && co<=r
                        hfi=pi/n;
                        fi=(0:n)*hfi;
                        ec=(0:n)*0;
                        for q=1:(n+1)
                              g=(co*cos(fi(q))+sqrt(r^2-co^2*sin(fi(q))^2));
                              a=atan(g/(d1+l)); %primer angulo de integracion
                              b=atan(g/d1); %segundo angulo de integracion
                              e1=0;
                              e2=0;
                              if 0<a
                                h1=a/n;
                                te=(0:n)*h1;
                                x=l./cos(te);
                                %ead=exp(-uventana*denv*xv/cos(te)); % autoabsorcion del detector
                                %f1=(1-exp(-udetector*den.*x)).*sin(te).*ead;
                                f1=(1-exp(-udetector*den.*x)).*sin(te).*exp(-uventana*denv*xv./cos(te));
                                e1=trapecio(f1,h1);
                              end
                              if a<b
                                h2=(b-a)/n;
                                te=a+(0:n)*h2;
                                x=g./sin(te)-d1./cos(te);
                                f2=(1-exp(-udetector*den.*x)).*sin(te).*exp(-uventana*denv*xv./cos(te));
                                e2=trapecio(f2,h2);
                              end
                              ec(q)=e1+e2;
                        end
                        e12=trapecio(ec,hfi);
                    elseif co>r 
                        hfi=asin(r/co)/n;
                        fi=(0:n)*hfi;
                        ec=(0:n)*0;
                        for q=1:(n+1)
                            g2=(co*cos(fi(q))-sqrt(r^2-co^2*sin(fi(q))^2));
                            g=(co*cos(fi(q))+sqrt(r^2-co^2*sin(fi(q))^2));
                            a=atan(g2/d1); %angulo de partida
                            b=atan(g/(d1+l)); %primer angulo de integracion
                            c=atan(g/d1); %segundo angulo de integracion
                            e1=0; %valor de arranque eficiencia de la primera integral
                            e2=0; %valor de arranque eficiencia de la segunda integral
                            if a<b 
                                h1=(b-a)/n;
                                te=a+(0:n)*h1;
                                x=l./cos(te);
                                f1=(1-exp(-udetector*den.*x)).*sin(te).*exp(-uventana*denv*xv./cos(te));
                                e1=trapecio(f1,h1);
                            end
                            if b<c && a<b
                                h2=(c-b)/n;
                                te=b+(0:n)*h2;
                                x=g./sin(te)-d1./cos(te);
                                f2=(1-exp(-udetector*den.*x)).*sin(te).*exp(-uventana*denv*xv./cos(te));
                                e2=trapecio(f2,h2);
                            end
                            if a>b && a<c
                                h2=(c-a)/n;
                                te=a+(0:n)*h2;
                                x=g./sin(te)-d1./cos(te);
                                f2=(1-exp(-udetector*den.*x)).*sin(te).*exp(-uventana*denv*xv./cos(te));
                                e2=trapecio(f2,h2);
                            end
                            ec(q)=e1+e2;
                        end
                        e12=trapecio(ec,hfi);
                    end
                    e(k,i)=e12*eang*100;
                end
            end
 final=clock;
%IMPRESION DE LOS RESULTADOS EN LA VENTANA DE MATLAB Y GRAFICO
 sec=final(6)-inicio(6);
    min=final(5)-inicio(5);
    if sec<0
        sec=final(6)+60-inicio(6);
        min=min-1;
    end

fprintf('Las eficiencias calculadas son:\n');
for k=1:j
    fprintf('Para una distancia Fuente-Detector de %3.1f cm se tiene: \n',d(k));
    for i=1:ni
        fprintf('\tPara  %4.3f MeV: e= %5.3f %%\n',E(i),e(k,i));
    end
end

if min>1
        fprintf('Tiempo de procesado: %d minutos con %3.1f segundos\n',min,sec);
    elseif min==0
        fprintf('Tiempo de procesado: %3.1f segundos\n',sec);
    elseif min==1
        fprintf('Tiempo de procesado: %d minuto con %3.1f segundos\n',min,sec);
    end
for k=1:j
    plot(E,e(k,:))
    hold on;
end
grid on;
hold off;
