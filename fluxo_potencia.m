%% Resolu��o do Fluxo de Pot�ncia via Newton Raphson

%% Defini��o das vari�veis calculadas

% inicializa��o das vari�veis 

%---------------- Pot�ncia ativa e reativa --------------

Pesp = zeros(nbus,1);
Qesp = zeros(nbus,1);

%---------------- M�dulo e �ngulo da tens�o -------------

V = zeros(nbus,1);
Th = zeros(nbus,1);

%---------------- Inicializa��o -------------------------

for k=1:nbus
    if busdata(k,2)==1
        Pesp(k,1) = 0;
        Qesp(k,1) = 0;
        V(k,1) = busdata(k,3);
        Th(k,1) = busdata(k,4)*pi/180;
    elseif busdata(k,2)==2
        Pesp(k,1) = (busdata(k,5) - busdata(k,7))/Sbase;
        Qesp(k,1) = 0;
        V(k,1) = busdata(k,3);
        Th(k,1) = 0;
    else
        Pesp(k,1) = (busdata(k,5) - busdata(k,7))/Sbase;
        Qesp(k,1) = (busdata(k,6) - busdata(k,8))/Sbase;
        V(k,1) = 1;
        Th(k,1) = 0;
    end
end

%% C�lculo do chute inicial 

Pcalc = zeros(nbus,1);
Qcalc = zeros(nbus,1);

for k=1:nbus
    if busdata(k,2)==1
        Pcalc(k,1)=0;
        Qcalc(k,1)=0;
    elseif busdata(k,2)==2
        Qcalc(k,1)=0;
        for m=1:nbus
            Pcalc(k,1) = Pcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
        end
    else
        for m=1:nbus
            Pcalc(k,1) = Pcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
            Qcalc(k,1) = Qcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*sin(Th(k,1) - Th(m,1)) - B(k,m)*cos(Th(k,1) - Th(m,1)));
        end
    end
end

DY = [Pesp-Pcalc;Qesp-Qcalc];

%% Defini��o da toler�ncia

tol = 1e-5;

%% Processo iterativo - Newton Raphson

while(max(abs(DY)))>tol
    
    %----------------- C�lculo da Matriz Jacobiana ----------------------
    
    % Jacobiana = [H N;M L];
    
    %----------------- C�lculo da Submatriz H - (dP/dTh) ----------------
    
    H = zeros(nbus,nbus);
    
    for k=1:nbus
        for m=1:nbus
            if k==m
                for j=1:nbus
                    H(k,m) = H(k,m) + V(k,1)*V(j,1)*(-G(k,j)*sin(Th(k,1) - Th(j,1)) + B(k,j)*cos(Th(k,1) - Th(j,1)));
                end
                H(k,m) = H(k,m) - (V(k,1)^2)*B(k,m);
            else
                H(k,m) = H(k,m) + V(k,1)*V(m,1)*(G(k,m)*sin(Th(k,1) - Th(m,1)) - B(k,m)*cos(Th(k,1) - Th(m,1)));
            end
        end
    end
    
    %----------------- C�lculo da Submatriz N - (dP/dV) ----------------
    
    N = zeros(nbus,nbus);
    
    for k=1:nbus
        for m=1:nbus
            if k==m
                for j=1:nbus
                    N(k,m) = N(k,m) + V(j,1)*(G(k,j)*cos(Th(k,1) - Th(j,1)) + B(k,j)*sin(Th(k,1) - Th(j,1)));
                end
                N(k,m) = N(k,m) + (V(k,1))*G(k,m);
            else
                N(k,m) = N(k,m) + V(k,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
            end
        end
    end
 
    %----------------- C�lculo da Submatriz M - (dQ/dTh) ----------------
    
    M = zeros(nbus,nbus);
    
    for k=1:nbus
        for m=1:nbus
            if k==m
                for j=1:nbus
                    M(k,m) = M(k,m) + V(k,1)*V(j,1)*(G(k,j)*cos(Th(k,1) - Th(j,1)) + B(k,j)*sin(Th(k,1) - Th(j,1)));
                end
                M(k,m) = M(k,m) - (V(k,1)^2)*G(k,m);
            else
                M(k,m) = M(k,m) + V(k,1)*V(m,1)*(-G(k,m)*cos(Th(k,1) - Th(m,1)) - B(k,m)*sin(Th(k,1) - Th(m,1)));
            end
        end
    end
    
    
    %----------------- C�lculo da Submatriz N - (dP/dV) ----------------
    
    L = zeros(nbus,nbus);
    
    for k=1:nbus
        for m=1:nbus
            if k==m
                for j=1:nbus
                    L(k,m) = L(k,m) + V(j,1)*(G(k,j)*sin(Th(k,1) - Th(j,1)) - B(k,j)*cos(Th(k,1) - Th(j,1)));
                end
                L(k,m) = L(k,m) - (V(k,1))*B(k,m);
            else
                L(k,m) = L(k,m) + V(k,1)*(G(k,m)*sin(Th(k,1) - Th(m,1)) - B(k,m)*cos(Th(k,1) - Th(m,1)));
            end
        end
    end
    
    Jacobiana = [H N;M L];
    
    % acr�scimo do Big Number
    
    for k=1:nbus
        if busdata(k,2)==1
            Jacobiana(k,k) = 1e10;
            Jacobiana(k+nbus,k+nbus) = 1e10;
        elseif busdata(k,2)==2
            Jacobiana(k+nbus,k+nbus) = 1e10;
        end
    end
    
    %-------------------- Resolu��o do sistema linear -------------------
    
    DX = Jacobiana\DY;
    
    %-------------------- Atualiza��o das vari�veis ---------------------
    
    for k=1:nbus
        Th(k,1) = Th(k,1) + DX(k,1);
        V(k,1) = V(k,1) + DX(k+nbus,1);
    end
    
    %------------------- C�lculo das pot�ncia injetadas -----------------
    
    Pcalc = zeros(nbus,1);
    Qcalc = zeros(nbus,1);

    for k=1:nbus
        if busdata(k,2)==1
            Pcalc(k,1)=0;
            Qcalc(k,1)=0;
        elseif busdata(k,2)==2
            Qcalc(k,1)=0;
            for m=1:nbus
                Pcalc(k,1) = Pcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
            end
        else
            for m=1:nbus
                Pcalc(k,1) = Pcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
                Qcalc(k,1) = Qcalc(k,1) + V(k,1)*V(m,1)*(G(k,m)*sin(Th(k,1) - Th(m,1)) - B(k,m)*cos(Th(k,1) - Th(m,1)));
            end
        end
    end

    DY = [Pesp-Pcalc;Qesp-Qcalc];
    
end

%% C�lculo do Estado Operativo da Rede

% 1) M�dulo da Tens�o e �ngulo j� calculados

% 2) ------------------- Pk e Qk em cada barra ------------------

Pk = zeros(nbus,1);
Qk = zeros(nbus,1);

for k=1:nbus
    for m=1:nbus
        Pk(k,1) = Pk(k,1) + V(k,1)*V(m,1)*(G(k,m)*cos(Th(k,1) - Th(m,1)) + B(k,m)*sin(Th(k,1) - Th(m,1)));
        Qk(k,1) = Qk(k,1) + V(k,1)*V(m,1)*(G(k,m)*sin(Th(k,1) - Th(m,1)) - B(k,m)*cos(Th(k,1) - Th(m,1)));
    end
end

% 2) ------------------- Pkm e Qkm em cada linha ------------------

Pkm = zeros(2*size(linedata,1),3);
Qkm = zeros(2*size(linedata,1),3);

for k=1:size(linedata,1)
    from = linedata(k,1);
    to = linedata(k,2);
    gkm = real(1/(linedata(k,3) + 1j*linedata(k,4)));        
    bkm = imag(1/(linedata(k,3) + 1j*linedata(k,4))); 
    % a)Pkm
    Pkm(2*k-1,1) = from;
    Pkm(2*k-1,2) = to;
    Pkm(2*k-1,3) = (linedata(k,6)^2)*gkm*(V(from,1)^2) - linedata(k,6)*bkm*V(from,1)*V(to,1)*sin(Th(from,1) - Th(to,1) + (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V(from,1)*V(to,1)*cos(Th(from,1) - Th(to,1) + (linedata(k,7)*pi/180));
    Pkm(2*k,1) = to;
    Pkm(2*k,2) = from;
    Pkm(2*k,3) = gkm*(V(to,1)^2) - linedata(k,6)*bkm*V(to,1)*V(from,1)*sin(Th(to,1) - Th(from,1) - (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V(to,1)*V(from,1)*cos(Th(to,1) - Th(from,1) - (linedata(k,7)*pi/180));
    % b)Qkm
    Qkm(2*k-1,1) = from;
    Qkm(2*k-1,2) = to;
    Qkm(2*k-1,3) = -(linedata(k,6)^2)*bkm*(V(from,1)^2) - (linedata(k,5)+busdata(from,11))*(V(from,1)^2) - linedata(k,6)*gkm*V(from,1)*V(to,1)*sin(Th(from,1) - Th(to,1) + (linedata(k,7)*pi/180)) + linedata(k,6)*bkm*V(from,1)*V(to,1)*cos(Th(from,1) - Th(to,1) + (linedata(k,7)*pi/180));
    Qkm(2*k,1) = to;
    Qkm(2*k,2) = from;
    Qkm(2*k,3) = -bkm*(V(to,1)^2) - (linedata(k,5)+busdata(to,11))*(V(to,1)^2) + linedata(k,6)*bkm*V(from,1)*V(to,1)*cos(Th(to,1) - Th(from,1) - (linedata(k,7)*pi/180)) - linedata(k,6)*gkm*V(from,1)*V(to,1)*sin(Th(to,1) - Th(from,1) - (linedata(k,7)*pi/180));
end

% 3) ------------------- Ikm e dIkm em cada linha ------------------

Ikm = zeros(2*size(linedata,1),3);
Del_Ikm = zeros(2*size(linedata,1),3);

for k=1:size(Pkm,1)
    from = Pkm(k,1);
    to = Pkm(k,2);
    Skm = Pkm(k,3) + 1j*Qkm(k,3);
    Ikm(k,1) = from;
    Ikm(k,2) = to;
    Ikm(k,3) = abs(conj((Skm)/(V(from,1)*exp(1j*Th(from,1)))));
    Del_Ikm(k,1) = from;
    Del_Ikm(k,2) = to;
    Del_Ikm(k,3) = angle(conj((Skm)/(V(from,1)*exp(1j*Th(from,1)))));
end


