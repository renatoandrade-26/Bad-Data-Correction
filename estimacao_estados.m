function [V_est,Th_est,Pk_est,Qk_est,Pkm_est,Qkm_est,Ikm_est,Del_Ikm_est,J,H,Ganho,erro] = estimacao_estados(G,B,busdata,linedata,Plano_medicao,R)

% Função que faz o cálculo da estimação de estados

%   Dados de entrada: matriz condutância  (G)
%                     matriz susceptância (B)
%                     dados de barra      (busdata)
%                     dados de linha      (linedata)
%                     plano de medição    (Plano_medicao)


%% ------------------ Variáveis estimadas ----------------------------

nbus = size(busdata,1);

Th_est = zeros(nbus,1); 
V_est = ones(nbus,1);
% inclusive a barra de referência
ref = find(busdata(:,2)==1);
V_est(ref,1) = busdata(ref,3);
Th_est(ref,1) = busdata(ref,4)*pi/180;
% V_est_ref,Th_est_ref = valor definido

%% ----------------- Vetor de erro das variáveis estimadas ----------------

DX_est = ones(2*nbus,1);

%% ------------- Processo iterativo da estimação de estados ---------------


tol = 1e-5;

iteracao = 0;


while max(abs(DX_est))>tol
    
    %-------------- Cálculo da Matriz Jacobiana das medições --------------
    
    H = zeros(size(Plano_medicao,1),2*nbus);
    
    % H = [dMedicao/dTh dMedicao/dV]
    
    for k=1:size(Plano_medicao,1)
        if Plano_medicao(k,1)==1
            fr = Plano_medicao(k,2);
            H(k,nbus+fr) = 1;
        elseif Plano_medicao(k,1)==2
            fr = Plano_medicao(k,2);
            H(k,fr) = 1;   
        elseif Plano_medicao(k,1)==3
            fr = Plano_medicao(k,2);
            for j=1:nbus
                if j==fr
                    for p=1:nbus
                        H(k,j) = H(k,j) + V_est(fr,1)*V_est(p,1)*(-G(fr,p)*sin(Th_est(fr,1) - Th_est(p,1)) + B(fr,p)*cos(Th_est(fr,1) - Th_est(p,1)));
                        H(k,nbus+j) = H(k,nbus+j) + V_est(p,1)*(G(fr,p)*cos(Th_est(fr,1) - Th_est(p,1)) + B(fr,p)*sin(Th_est(fr,1) - Th_est(p,1)));
                    end
                    H(k,j) = H(k,j) - (V_est(fr,1)^2)*B(fr,j);
                    H(k,nbus+j) = H(k,nbus+j) + V_est(fr,1)*G(fr,j);
                else
                    H(k,j) = H(k,j) + V_est(fr,1)*V_est(j,1)*(G(fr,j)*sin(Th_est(fr,1) - Th_est(j,1)) - B(fr,j)*cos(Th_est(fr,1) - Th_est(j,1)));
                    H(k,nbus+j) = H(k,nbus+j) + V_est(fr,1)*(G(fr,j)*cos(Th_est(fr,1) - Th_est(j,1)) + B(fr,j)*sin(Th_est(fr,1) - Th_est(j,1)));
                end
            end
        elseif Plano_medicao(k,1)==4
            fr = Plano_medicao(k,2);
            for j=1:nbus
                if j==fr
                    for p=1:nbus
                        H(k,j) = H(k,j) + V_est(fr,1)*V_est(p,1)*(G(fr,p)*cos(Th_est(fr,1) - Th_est(p,1)) + B(fr,p)*sin(Th_est(fr,1) - Th_est(p,1)));
                        H(k,nbus+j) = H(k,nbus+j) + V_est(p,1)*(G(fr,p)*sin(Th_est(fr,1) - Th_est(p,1)) - B(fr,p)*cos(Th_est(fr,1) - Th_est(p,1)));
                    end
                    H(k,j) = H(k,j) - (V_est(fr,1)^2)*G(fr,j);
                    H(k,nbus+j) = H(k,nbus+j) - V_est(fr,1)*B(fr,j);
                else
                    H(k,j) = H(k,j) + V_est(fr,1)*V_est(j,1)*(-G(fr,j)*cos(Th_est(fr,1) - Th_est(j,1)) - B(fr,j)*sin(Th_est(fr,1) - Th_est(j,1)));
                    H(k,nbus+j) = H(k,nbus+j) + V_est(fr,1)*(G(fr,j)*sin(Th_est(fr,1) - Th_est(j,1)) - B(fr,j)*cos(Th_est(fr,1) - Th_est(j,1)));
                end
            end
        elseif Plano_medicao(k,1)==5
            fr = Plano_medicao(k,2);
            to = Plano_medicao(k,3);
            gkm = 0;
            bkm = 0;
            akm = 0;
            thkm = 0;
            linha_data = 0;
            for s=1:size(linedata,1)
                if or(and(linedata(s,1)==fr,linedata(s,2)==to),and(linedata(s,1)==to,linedata(s,2)==fr))
                    gkm = real(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bkm = imag(1/(linedata(s,3) + 1j*linedata(s,4)));
                    akm = linedata(s,6);
                    thkm = linedata(s,7);
                    linha_data = s;
                end
            end
            if fr==linedata(linha_data,1)
                H(k,fr) = H(k,fr) - akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm) + akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,to) = H(k,to) + akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,fr+nbus) = H(k,fr+nbus) + 2*(akm^2)*gkm*V_est(fr,1) - akm*bkm*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*gkm*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,to+nbus) = H(k,to+nbus) - akm*bkm*V_est(fr,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*gkm*V_est(fr,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm);
            else
                H(k,fr) = H(k,fr) - akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm) + akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,to) = H(k,to) + akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,fr+nbus) = H(k,fr+nbus) + 2*gkm*V_est(fr,1) - akm*bkm*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,to+nbus) = H(k,to+nbus) - akm*bkm*V_est(fr,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(fr,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm);
            end
        elseif Plano_medicao(k,1)==6
            fr = Plano_medicao(k,2);
            to = Plano_medicao(k,3);
            gkm = 0;
            bkm = 0;
            bsh = 0;
            akm = 0;
            thkm = 0;
            linha_data = 0;
            for s=1:size(linedata,1)
                if or(and(linedata(s,1)==fr,linedata(s,2)==to),and(linedata(s,1)==to,linedata(s,2)==fr))
                    gkm = real(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bkm = imag(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bsh = linedata(s,5);
                    akm = linedata(s,6);
                    thkm = linedata(s,7);
                    linha_data = s;
                end
            end
            if fr==linedata(linha_data,1)
                H(k,fr) = H(k,fr) - akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,to) = H(k,to) + akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm) + akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,fr+nbus) = H(k,fr+nbus) + 2*V_est(fr,1)*(-(akm^2)*bkm - bsh - busdata(fr,11)) + akm*bkm*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*gkm*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm);
                H(k,to+nbus) = H(k,to+nbus) + akm*bkm*V_est(fr,1)*cos(Th_est(fr,1) - Th_est(to,1) + thkm) - akm*gkm*V_est(fr,1)*sin(Th_est(fr,1) - Th_est(to,1) + thkm);
            else
                H(k,fr) = H(k,fr) - akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,to) = H(k,to) + akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm) + akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,fr+nbus) = H(k,fr+nbus) + 2*V_est(fr,1)*(-bkm - bsh - busdata(fr,11)) + akm*bkm*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm);
                H(k,to+nbus) = H(k,to+nbus) + akm*bkm*V_est(fr,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(fr,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm);
            end 
        end
    end
    
    % eliminação da 1° coluna referente ao ângulo da barra de referência
    
    H(:,1) = [];
    
    % Cálculo da Matriz Ganho
    
    Ganho = transpose(H)*inv(R)*H;
    
    % 'Jacobiana' da estimação para atualização das variáveis estimadas
    
    Jacobiana_estimacao = inv(Ganho)*transpose(H)*inv(R);
    
    % Cálculo do resíduo entre o erro medido e calculado
    
    erro = zeros(size(Plano_medicao,1),1);
    
    for k=1:size(Plano_medicao,1)
        if Plano_medicao(k,1)==1
            fr = Plano_medicao(k,2);
            erro(k,1) = Plano_medicao(k,4) - V_est(fr,1);
        elseif Plano_medicao(k,1)==2
            fr = Plano_medicao(k,2);
            erro(k,1) = Plano_medicao(k,4) - Th_est(fr,1);
        elseif Plano_medicao(k,1)==3
            fr = Plano_medicao(k,2);
            Pk_est = 0;
            for s=1:nbus
                Pk_est = Pk_est + V_est(fr,1)*V_est(s,1)*(G(fr,s)*cos(Th_est(fr,1) - Th_est(s,1)) + B(fr,s)*sin(Th_est(fr,1) - Th_est(s,1)));
            end
            erro(k,1) = Plano_medicao(k,4) - Pk_est;
        elseif Plano_medicao(k,1)==4
            fr = Plano_medicao(k,2);
            Qk_est = 0;
            for s=1:nbus
                Qk_est = Qk_est + V_est(fr,1)*V_est(s,1)*(G(fr,s)*sin(Th_est(fr,1) - Th_est(s,1)) - B(fr,s)*cos(Th_est(fr,1) - Th_est(s,1)));
            end
            erro(k,1) = Plano_medicao(k,4) - Qk_est;
        elseif Plano_medicao(k,1)==5
            fr = Plano_medicao(k,2);
            to = Plano_medicao(k,3);
            gkm = 0;
            bkm = 0;
            akm = 0;
            thkm = 0;
            linha_data = 0;
            for s=1:size(linedata,1)
                if or(and(linedata(s,1)==fr,linedata(s,2)==to),and(linedata(s,1)==to,linedata(s,2)==fr))
                    gkm = real(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bkm = imag(1/(linedata(s,3) + 1j*linedata(s,4)));
                    akm = linedata(s,6);
                    thkm = linedata(s,7);
                    linha_data = s;
                end
            end
            if fr==linedata(linha_data,1)
                erro(k,1) = Plano_medicao(k,4) - ((akm^2)*gkm*V_est(fr,1)^2 - akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1)+thkm) - akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1)+thkm)); 
            else
                erro(k,1) = Plano_medicao(k,4) - (gkm*V_est(fr,1)^2 - akm*bkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm) - akm*gkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm)); 
            end
        elseif Plano_medicao(k,1)==6
            fr = Plano_medicao(k,2);
            to = Plano_medicao(k,3);
            gkm = 0;
            bkm = 0;
            bsh = 0;
            akm = 0;
            thkm = 0;
            linha_data = 0;
            for s=1:size(linedata,1)
                if or(and(linedata(s,1)==fr,linedata(s,2)==to),and(linedata(s,1)==to,linedata(s,2)==fr))
                    gkm = real(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bkm = imag(1/(linedata(s,3) + 1j*linedata(s,4)));
                    bsh = linedata(s,5);
                    akm = linedata(s,6);
                    thkm = linedata(s,7);
                    linha_data = s;
                end
            end
            if fr==linedata(linha_data,1)
                erro(k,1) = Plano_medicao(k,4) - ((-(akm^2)*bkm - bsh - busdata(fr,11))*V_est(fr,1)^2 - akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1)+thkm) + akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1)+thkm)); 
            else
                erro(k,1) = Plano_medicao(k,4) - ((- bkm - bsh - busdata(fr,11))*V_est(fr,1)^2 - akm*gkm*V_est(fr,1)*V_est(to,1)*sin(Th_est(fr,1) - Th_est(to,1) - thkm) + akm*bkm*V_est(fr,1)*V_est(to,1)*cos(Th_est(fr,1) - Th_est(to,1) - thkm)); 
            end
        end
    end
    
    % Somatório dos erros quadráticos ponderados
    
    J = transpose(erro)*inv(R)*erro;
    
    % Variação das variáveis estimadas
    
    DX_est = Jacobiana_estimacao*erro;
    DX_est = [0;DX_est];
    
    % Atualização das variáveis estimadas;
    
    for k=1:nbus
        V_est(k,1) = V_est(k,1) + DX_est(k+nbus,1);
        Th_est(k,1) = Th_est(k,1) + DX_est(k,1);
    end
    
    iteracao = iteracao + 1;
    
    if iteracao>20
        break;
    end
end


%% Cálculo do Estado Estimado da Rede

% 1) Módulo da Tensão e Ângulo Estimados já calculados

% 2) ------------------- Pk e Qk em cada barra ------------------

Pk_est = zeros(nbus,1);
Qk_est = zeros(nbus,1);

for k=1:nbus
    for m=1:nbus
        Pk_est(k,1) = Pk_est(k,1) + V_est(k,1)*V_est(m,1)*(G(k,m)*cos(Th_est(k,1) - Th_est(m,1)) + B(k,m)*sin(Th_est(k,1) - Th_est(m,1)));
        Qk_est(k,1) = Qk_est(k,1) + V_est(k,1)*V_est(m,1)*(G(k,m)*sin(Th_est(k,1) - Th_est(m,1)) - B(k,m)*cos(Th_est(k,1) - Th_est(m,1)));
    end
end

% 2) ------------------- Pkm e Qkm em cada linha ------------------

Pkm_est = zeros(2*size(linedata,1),3);
Qkm_est = zeros(2*size(linedata,1),3);

for k=1:size(linedata,1)
    from = linedata(k,1);
    to = linedata(k,2);
    gkm = real(1/(linedata(k,3) + 1j*linedata(k,4)));        
    bkm = imag(1/(linedata(k,3) + 1j*linedata(k,4))); 
    % a)Pkm
    Pkm_est(2*k-1,1) = from;
    Pkm_est(2*k-1,2) = to;
    Pkm_est(2*k-1,3) = (linedata(k,6)^2)*gkm*(V_est(from,1)^2) - linedata(k,6)*bkm*V_est(from,1)*V_est(to,1)*sin(Th_est(from,1) - Th_est(to,1) + (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V_est(from,1)*V_est(to,1)*cos(Th_est(from,1) - Th_est(to,1) + (linedata(k,7)*pi/180));
    Pkm_est(2*k,1) = to;
    Pkm_est(2*k,2) = from;
    Pkm_est(2*k,3) = gkm*(V_est(to,1)^2) - linedata(k,6)*bkm*V_est(to,1)*V_est(from,1)*sin(Th_est(to,1) - Th_est(from,1) - (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V_est(to,1)*V_est(from,1)*cos(Th_est(to,1) - Th_est(from,1) - (linedata(k,7)*pi/180));
    % b)Qkm
    Qkm_est(2*k-1,1) = from;
    Qkm_est(2*k-1,2) = to;
    Qkm_est(2*k-1,3) = -(linedata(k,6)^2)*bkm*(V_est(from,1)^2) - (linedata(k,5)+busdata(from,11))*(V_est(from,1)^2) - linedata(k,6)*gkm*V_est(from,1)*V_est(to,1)*sin(Th_est(from,1) - Th_est(to,1) + (linedata(k,7)*pi/180)) + linedata(k,6)*bkm*V_est(from,1)*V_est(to,1)*cos(Th_est(from,1) - Th_est(to,1) + (linedata(k,7)*pi/180));
    Qkm_est(2*k,1) = to;
    Qkm_est(2*k,2) = from;
    Qkm_est(2*k,3) = -bkm*(V_est(to,1)^2) - (linedata(k,5)+busdata(to,11))*(V_est(to,1)^2) + linedata(k,6)*bkm*V_est(from,1)*V_est(to,1)*cos(Th_est(to,1) - Th_est(from,1) - (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V_est(from,1)*V_est(to,1)*sin(Th_est(to,1) - Th_est(from,1) - (linedata(k,7)*pi/180));
end

% 3) ------------------- Ikm e dIkm em cada linha ------------------

Ikm_est = zeros(2*size(linedata,1),3);
Del_Ikm_est = zeros(2*size(linedata,1),3);

for k=1:size(Pkm_est,1)
    from = Pkm_est(k,1);
    to = Pkm_est(k,2);
    Skm_est = Pkm_est(k,3) + 1j*Qkm_est(k,3);
    Ikm_est(k,1) = from;
    Ikm_est(k,2) = to;
    Ikm_est(k,3) = abs(conj((Skm_est)/(V_est(from,1)*exp(1j*Th_est(from,1)))));
    Del_Ikm_est(k,1) = from;
    Del_Ikm_est(k,2) = to;
    Del_Ikm_est(k,3) = angle(conj((Skm_est)/(V_est(from,1)*exp(1j*Th_est(from,1)))));
end
end




