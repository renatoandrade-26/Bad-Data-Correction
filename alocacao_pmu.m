%% Alocação ótima de PMU baseado na maximização de medidas estimadas e conexões de barras

%% Determinação do n° de conexões das barras

Num_conexao = zeros(size(Conexao,1),2);

for k=1:size(Conexao,1)
    Num_conexao(k,1) = k;
    Num_conexao(k,2) = sum(Conexao(k,:)) - 1;
end

Num_conexao = sortrows(Num_conexao,2,'descend');

%% Enumeração exaustiva

bus_estimadas = [];
bus_pmu = [];

for k=1:size(Num_conexao,1)
    
    %------------------ Possível barra para a PMU ser alocada -------------
    
    possivel_barra = Num_conexao(k,1);
    
    %-------------- Verificação se a barra já pode ser estimada -----------
    
    if possivel_barra~=find(busdata(:,2)==1)
    
        flag = size(find(bus_estimadas==possivel_barra),2);

        if flag==0
            bus_pmu(size(bus_pmu,2)+1) = possivel_barra;
            for j=1:size(linedata,1)
                if linedata(j,1)==possivel_barra 
                    possivel_est = linedata(j,2);
                    if size(find(bus_estimadas==possivel_est),2)==0
                        bus_estimadas(1,size(bus_estimadas,2)+1) = possivel_est;
                    end
                elseif linedata(j,2)==possivel_barra 
                    possivel_est = linedata(j,1);
                     if size(find(bus_estimadas==possivel_est),2)==0
                        bus_estimadas(1,size(bus_estimadas,2)+1) = possivel_est;
                     end
                end
            end
        end  
    end
end

bus_pmu(size(bus_pmu,2)+1) = find(busdata(:,2)==1);