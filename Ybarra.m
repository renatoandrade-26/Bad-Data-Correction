%% Montagem da Matriz Ybarra - Modelo completo - Tap + Defasador + Linha

%% Cálculo

Ybus = zeros(nbus,nbus);
for i=1:nbus
    for j=1:nbus
        if i==j
            for k=1:size(linedata,1)
                if linedata(k,1)==i
                    Ybus(i,j) = Ybus(i,j) + ((linedata(k,6)^2)*(1/((linedata(k,3)) + 1j*linedata(k,4))));
                elseif linedata(k,2)==i
                    Ybus(i,j) = Ybus(i,j) + (1)*(1/(linedata(k,3) + 1j*linedata(k,4)));
                end
            end
        else
            for k=1:size(linedata,1)
                if linedata(k,1)==i && linedata(k,2)==j
                    Ybus(i,j) = Ybus(i,j) - ((linedata(k,6))*(1/(linedata(k,3) + 1j*linedata(k,4))))*exp(-1j*linedata(k,7));
                elseif linedata(k,1)==j && linedata(k,2)==i
                    Ybus(i,j) = Ybus(i,j) - ((linedata(k,6))*(1/(linedata(k,3) + 1j*linedata(k,4))))*exp(1j*linedata(k,7));
                end
            end
        end
    end
end

% Elementos shunts de barra

for k=1:nbus
    Ybus(k,k) = Ybus(k,k) + 1j*busdata(k,11)/100;
end

% Elementos shunts de linha

for k=1:size(linedata,1)
    Ybus(linedata(k,1),linedata(k,1)) = Ybus(linedata(k,1),linedata(k,1)) + 1j*linedata(k,5);
    Ybus(linedata(k,2),linedata(k,2)) = Ybus(linedata(k,2),linedata(k,2)) + 1j*linedata(k,5);
end
%% Matriz Condutância e Matriz Susceptância

G = real(Ybus);
B = imag(Ybus);

%% Matriz Conexão

Conexao = zeros(size(Ybus,1),size(Ybus,1));

for k=1:size(linedata,1)
    de = linedata(k,1);
    para = linedata(k,2);
    Conexao(de,de) = 1;
    Conexao(de,para) = 1;
    Conexao(para,de) = 1;
    Conexao(para,para) = 1;
end
