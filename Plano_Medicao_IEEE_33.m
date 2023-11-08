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

N_medidas_Vk = 16;
N_medidas_Thk = 16;
N_medidas_Pk = 17;
N_medidas_Qk = 17;
N_medidas_Pkm = 32;
N_medidas_Qkm = 32;

%---------------------- Tipo de medição -------------------------

tipo = [ones(1,N_medidas_Vk) 2*ones(1,N_medidas_Thk) 3*ones(1,N_medidas_Pk) 4*ones(1,N_medidas_Qk) 5*ones(1,N_medidas_Pkm) 6*ones(1,N_medidas_Qkm)];

%---------------------- Barras de medição ------------------------

fr = [1 3 5 7 9 11 13 15 17 19 21 24 26 28 30 32 1 3 5 7 9 11 13 15 17 19 21 24 26 28 30 32 2 4 6 8 10 12 14 16 18 20 22 23 25 27 29 31 33 2 4 6 8 10 12 14 16 18 20 22 23 25 27 29 31 33 1 2 3 3 4 5 6 7 8 9 11 10 13 12 15 14 17 16 2 19 20 21 23 24 26 6 27 28 29 30 31 32 1 2 3 3 4 5 6 7 8 9 11 10 13 12 15 14 17 16 2 19 20 21 23 24 26 6 27 28 29 30 31 32];
to = [1 3 5 7 9 11 13 15 17 19 21 24 26 28 30 32 1 3 5 7 9 11 13 15 17 19 21 24 26 28 30 32 2 4 6 8 10 12 14 16 18 20 22 23 25 27 29 31 33 2 4 6 8 10 12 14 16 18 20 22 23 25 27 29 31 33 2 3 4 23 5 6 7 8 9 10 12 11 14 13 16 15 18 17 19 20 21 22 24 25 27 26 28 29 30 31 32 33 2 3 4 23 5 6 7 8 9 10 12 11 14 13 16 15 18 17 19 20 21 22 24 25 27 26 28 29 30 31 32 33];

%---------------------- Desvio da medição ------------------------

desvio = [0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008,0.008];



%---------------------- Plano de medição --------------------------

%                      Tipo          Barra DE     Barra PARA       VALOR            DESVIO
Plano_medicao = [transpose(tipo) transpose(fr) transpose(to) zeros(size(tipo,2),1) transpose(desvio)];

for k=1:size(Plano_medicao,1)
    ruido = (1+0.00015*randn);
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

for k=1:size(Plano_medicao,1)
    if Plano_medicao(k,1)==1 && Plano_medicao(k,2)==9
        Plano_medicao(k,4) = 1.3*Plano_medicao(k,4);
    end
end