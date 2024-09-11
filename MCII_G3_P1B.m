clearvars, close all

ncoef = 20;

% tomando n = 1
k = 2:2:ncoef; % los k impares dan coef = 0
Hpkn = sin((k-1)/2*pi)./(k-1)-sin((k+1)/2*pi)./(k+1); % elementos H'

% se muestran los sucesivos valores de H'_k1
figure
plot(k,Hpkn,'-*')
xlim([0 ncoef])
xlabel('k')
ylabel('<k⁰|Hp|1⁰>')

suma = [];
% se calcula la suma de |H'k1|²/(1-k²) hasta cada valor de k
for i = 1:length(Hpkn)
    suma(i) = sum(Hpkn(1:i).^2./(1-k(1:i).^2));
    disp(['coef: ',num2str(Hpkn(i).^2./(1-k(i).^2)),' - suma: ',num2str(suma(i))])
end

% se muestra la evolución de la suma
figure
plot(k,suma)
xlabel('k')
ylabel('suma')

% tomando n = 2
k = 1:2:ncoef; % los k impares dan coef = 0
Hpkn = sin((k-2)/2*pi)./(k-2)-sin((k+2)/2*pi)./(k+2); % elementos H'

% se muestran los sucesivos valores de H'_k1
figure
plot(k,Hpkn,'-*')
% xlim([0 500])
xlabel('k')
ylabel('<k⁰|Hp|2⁰>')

suma = [];
% se calcula la suma de |H'k2|²/(4-k²) hasta cada valor de k
for i = 1:length(Hpkn)
    suma(i) = sum(Hpkn(1:i).^2./(4-k(1:i).^2));
    disp(['coef: ',num2str(Hpkn(i).^2./(4-k(i).^2)),' - suma: ',num2str(suma(i))])
end

% se muestra la evolución de la suma
figure
plot(k,suma)
xlabel('k')
ylabel('suma')
