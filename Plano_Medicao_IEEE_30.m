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

% N_medidas_Vk = 21;
% N_medidas_Thk = 21;
% N_medidas_Pk = 0;
% N_medidas_Qk = 0;
% N_medidas_Pkm = 56;
% N_medidas_Qkm = 56;

% N_medidas_Vk = 10;
% N_medidas_Thk = 10;
% N_medidas_Pk = 20;
% N_medidas_Qk = 20;
% N_medidas_Pkm = 25;
% N_medidas_Qkm = 25;

% N_medidas_Vk = 20;
% N_medidas_Thk = 20;
% N_medidas_Pk = 0;
% N_medidas_Qk = 0;
% N_medidas_Pkm = 49;
% N_medidas_Qkm = 49;

%--------- Alocação feita --------
N_medidas_Vk = 13;
N_medidas_Thk = 13;
N_medidas_Pk = 0;
N_medidas_Qk = 0;
N_medidas_Pkm = 35;
N_medidas_Qkm = 35;

%---------------------- Tipo de medição -------------------------

% tipo = [ones(1,N_medidas_Vk) 2*ones(1,N_medidas_Thk) 3*ones(1,N_medidas_Pk) 4*ones(1,N_medidas_Qk) 5*ones(1,N_medidas_Pkm) 6*ones(1,N_medidas_Qkm)];
tipo = [ones(1,N_medidas_Vk) 2*ones(1,N_medidas_Thk) 5*ones(1,N_medidas_Pkm) 6*ones(1,N_medidas_Qkm)];

%---------------------- Barras de medição ------------------------

% fr = [1 5 8 10 11 12 18 23 26 30 1 5 8 10 11 12 18 23 26 30 2 3 4 6 7 9 13 14 15 16 17 19 20 21 22 24 25 27 28 29 2 3 4 6 7 9 13 14 15 16 17 19 20 21 22 24 25 27 28 29 1 1 5 5 8 8 10 10 10 10 10 10 11 12 12 12 12 12 18 18 23 23 26 30 30 1 1 5 5 8 8 10 10 10 10 10 10 11 12 12 12 12 12 18 18 23 23 26 30 30];
% to = [1 5 8 10 11 12 18 23 26 30 1 5 8 10 11 12 18 23 26 30 2 3 4 6 7 9 13 14 15 16 17 19 20 21 22 24 25 27 28 29 2 3 4 6 7 9 13 14 15 16 17 19 20 21 22 24 25 27 28 29 2 3 2 7 6 28 6 9 17 20 21 22 9 4 13 14 15 16 15 19 15 24 25 27 29 2 3 2 7 6 28 6 9 17 20 21 22 9 4 13 14 15 16 15 19 15 24 25 27 29];

% fr = [1 3 5 7 8 9 10 11 12 13 15 17 18 20 22 24 25 27 28 29 30 1 3 5 7 8 9 10 11 12 13 15 17 18 20 22 24 25 27 28 29 30 1 1 3 3 5 5 7 7 8 8 9 9 9 10 10 10 10 10 10 11 12 12 12 12 12 13 15 15 15 15 17 17 18 18 20 20 22 22 22 24 24 24 25 25 25 27 27 27 27 28 28 28 29 29 30 30 1 1 3 3 5 5 7 7 8 8 9 9 9 10 10 10 10 10 10 11 12 12 12 12 12 13 15 15 15 15 17 17 18 18 20 20 22 22 22 24 24 24 25 25 25 27 27 27 27 28 28 28 29 29 30 30];
% to = [1 3 5 7 8 9 10 11 12 13 15 17 18 20 22 24 25 27 28 29 30 1 3 5 7 8 9 10 11 12 13 15 17 18 20 22 24 25 27 28 29 30 2 3 1 4 2 7 5 6 6 28 6 10 11 6 9 17 20 21 22 9 4 13 14 15 16 12 12 14 18 23 10 16 15 19 10 19 10 21 24 22 23 25 24 26 27 25 28 29 30 6 8 27 27 30 27 29 2 3 1 4 2 7 5 6 6 28 6 10 11 6 9 17 20 21 22 9 4 13 14 15 16 12 12 14 18 23 10 16 15 19 10 19 10 21 24 22 23 25 24 26 27 25 28 29 30 6 8 27 27 30 27 29];

% fr = [1 3 5 7 8 9 10 11 12 13 15 17 18 20 21 24 26 28 29 30 1 3 5 7 8 9 10 11 12 13 15 17 18 20 21 24 26 28 29 30 1 1 3 3 5 5 7 7 8 8 9 9 9 10 10 10 10 10 10 11 12 12 12 12 12 13 15 15 15 15 17 17 18 18 20 20 21 21 24 24 24 26 28 28 28 29 29 30 30 1 1 3 3 5 5 7 7 8 8 9 9 9 10 10 10 10 10 10 11 12 12 12 12 12 13 15 15 15 15 17 17 18 18 20 20 21 21 24 24 24 26 28 28 28 29 29 30 30];
% to = [1 3 5 7 8 9 10 11 12 13 15 17 18 20 21 24 26 28 29 30 1 3 5 7 8 9 10 11 12 13 15 17 18 20 21 24 26 28 29 30 2 3 1 4 2 7 5 6 6 28 6 10 11 6 9 17 20 21 22 9 4 13 14 15 16 12 12 14 18 23 10 16 15 19 10 19 10 22 22 23 25 25 6 8 27 27 30 27 29 2 3 1 4 2 7 5 6 6 28 6 10 11 6 9 17 20 21 22 9 4 13 14 15 16 12 12 14 18 23 10 16 15 19 10 19 10 22 22 23 25 25 6 8 27 27 30 27 29];

fr = [1 3 5 6 11 12 17 18 20 22 23 26 27 1 3 5 6 11 12 17 18 20 22 23 26 27 1 1 3 3 5 5 6 6 6 6 6 6 6 11 12 12 12 12 12 17 17 18 18 20 20 22 22 22 23 23 26 27 27 27 27 1 1 3 3 5 5 6 6 6 6 6 6 6 11 12 12 12 12 12 17 17 18 18 20 20 22 22 22 23 23 26 27 27 27 27];
to = [1 3 5 6 11 12 17 18 20 22 23 26 27 1 3 5 6 11 12 17 18 20 22 23 26 27 2 3 1 4 2 7 2 4 7 8 9 10 28 9 4 13 14 15 16 10 16 15 19 10 19 10 21 24 15 24 25 25 28 29 30 2 3 1 4 2 7 2 4 7 8 9 10 28 9 4 13 14 15 16 10 16 15 19 10 19 10 21 24 15 24 25 25 28 29 30];
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
    ruido = (1+0.0005*randn);
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

for m=1:size(Plano_medicao,1)
    
    if Plano_medicao(m,1)==6 && (Plano_medicao(m,2)==12)
        
        Plano_medicao(m,4) = (2 + rand/1000)*Plano_medicao(m,4);
        
    end
    
end
