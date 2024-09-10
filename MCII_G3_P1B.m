clearvars, close all

k = 2:2:200; % los k impares dan coef = 0
Hpkn = sin((k-1)/2*pi)./(k-1)-sin((k+1)/2*pi)./(k+1); % elementos H'

% se muestran los sucesivos valores de H'_k1
figure
plot(k,Hpkn,'-*')
xlim([0 100])
xlabel('k')
ylabel('<k⁰|Hp|1⁰>')

% se calcula la suma de |H'k1|²/(1-k²) hasta cada valor de k
for i = 1:length(Hpkn)
    suma(i) = sum(Hpkn(1:i).^2./(1-k(1:i).^2));
end

% se muestra la evolución de la suma
figure
plot(k,suma)
xlabel('k')
ylabel('suma')
