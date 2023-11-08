%% Definição do plano de medição 

% Tipos de medições: 1 -  Vk
%                    2 -  Th_Vk
%                    3 -  Pk
%                    4 -  Qk
%                    5 -  Pkm
%                    6 -  Qkm
%                    7 -  Ikm
%                    8 -  Th_Ikm

%% Plano de medição

%---------------- Número de grandezas mensuradas -----------------

N_medidas_Vk = 2;
N_medidas_Thk = 2;
N_medidas_Pk = 0;
N_medidas_Qk = 0;
N_medidas_Pkm = 4;
N_medidas_Qkm = 4;

%---------------------- Tipo de medição -------------------------

tipo = [ones(1,N_medidas_Vk) 2*ones(1,N_medidas_Thk) 5*ones(1,N_medidas_Pkm) 6*ones(1,N_medidas_Qkm)];

%---------------------- Barras de medição ------------------------

fr = [1 2 1 2 1 2 2 2 1 2 2 2];
to = [1 2 1 2 2 1 3 4 2 1 3 4];
%---------------------- Desvio da medição ------------------------

desvio = zeros(1,size(tipo,2));
for k=1:size(tipo,2)
    if tipo(1,k)==1 || tipo(1,k)==2
        desvio(1,k) = 0.004;
    elseif tipo(1,k)==3 || tipo(1,k)==4
        desvio(1,k) = 0.010;
    else 
        desvio(1,k) = 0.008;
    end
end


%---------------------- Plano de medição --------------------------

%                      Tipo          Barra DE     Barra PARA       VALOR            DESVIO
Plano_medicao = [transpose(tipo) transpose(fr) transpose(to) zeros(size(tipo,2),1) transpose(desvio)];

for k=1:size(Plano_medicao,1)
    ruido = (1+0.001*randn);
    if Plano_medicao(k,1)==1
        Plano_medicao(k,4) = V(Plano_medicao(k,2))*ruido;
    elseif Plano_medicao(k,1)==2
        Plano_medicao(k,4) = Th(Plano_medicao(k,2))*ruido;
    elseif Plano_medicao(k,1)==3
        Plano_medicao(k,4) = Pk(Plano_medicao(k,2))*ruido;
    elseif Plano_medicao(k,1)==4
        Plano_medicao(k,4) = Qk(Plano_medicao(k,2))*ruido;
    elseif Plano_medicao(k,1)==5
        for s=1:size(Pkm,1)
            if Pkm(s,1)==Plano_medicao(k,2) && Pkm(s,2)==Plano_medicao(k,3)
                Plano_medicao(k,4) = Pkm(s,3)*ruido;
            end
        end
    elseif Plano_medicao(k,1)==6
        for s=1:size(Qkm,1)
            if Qkm(s,1)==Plano_medicao(k,2) && Qkm(s,2)==Plano_medicao(k,3)
                Plano_medicao(k,4) = Qkm(s,3)*ruido;
            end
        end
    elseif Plano_medicao(k,1)==7
        for s=1:size(Ikm,1)
            if Ikm(s,1)==Plano_medicao(k,2) && Ikm(s,2)==Plano_medicao(k,3)
                Plano_medicao(k,4) = Ikm(s,3)*ruido;
            end
        end        
    elseif Plano_medicao(k,1)==8
        for s=1:size(Del_Ikm,1)
            if Del_Ikm(s,1)==Plano_medicao(k,2) && Del_Ikm(s,2)==Plano_medicao(k,3)
                Plano_medicao(k,4) = Del_Ikm(s,3)*ruido;
            end
        end 
    end
end

% Matriz variância das medições
R = zeros(size(Plano_medicao,1),size(Plano_medicao,1));
for k=1:size(Plano_medicao,1)
    R(k,k) = Plano_medicao(k,5)^2;
end

% Parte para colocar erro grosseiro 

% for k=1:size(Plano_medicao,1)
%    if (Plano_medicao(k,1)==1 && Plano_medicao(k,2)==17)
%        Plano_medicao(k,4) = 0.8391*1.2*Plano_medicao(k,4);
%    elseif (Plano_medicao(k,1)==5 && Plano_medicao(k,2)==17 && Plano_medicao(k,3)==16) 
%        Plano_medicao(k,4) = 0.1322*5.0*Plano_medicao(k,4);
%    elseif (Plano_medicao(k,1)==5 && Plano_medicao(k,2)==17 && Plano_medicao(k,3)==10) 
%        Plano_medicao(k,4) = 0.1704*5.0*Plano_medicao(k,4);
%    elseif (Plano_medicao(k,1)==5 && Plano_medicao(k,2)==12 && Plano_medicao(k,3)==16) 
%        Plano_medicao(k,4) = 1.1171*1.0*Plano_medicao(k,4);
%    elseif (Plano_medicao(k,1)==5 && Plano_medicao(k,2)==12 && Plano_medicao(k,3)==16) 
%        Plano_medicao(k,4) = 1.0585*1.0*Plano_medicao(k,4);     
%    end
% end