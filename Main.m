%% CORREÇÃO DE ERROS GROSSEIROS NA ESTIMAÇÃO DE ESTADOS

%  Autor: Renato Andrade Mosqueira Furtado

%  Orientadores: Ivo Chaves da Silva Junior e Igor Delgado de Melo

%  Organização: Universidade Federal de Juiz de Fora

close all
clear all
clc
tic
%% Dados de entrada
Dados_Sistema_IEEE_30;

%% Cálculo da Matriz Ybarra
Ybarra;

%% Resolução do Fluxo de Potência via Newton Raphson
fluxo_potencia;

%% Alocação PMU
alocacao_pmu;

%% Plano de Medição 
Plano_Medicao_IEEE_30;

%% Função para Estimação de Estados

[V_est_antes,Th_est_antes,Pk_est_antes,Qk_est_antes,Pkm_est_antes,Qkm_est_antes,Ikm_est_antes,Del_Ikm_est_antes,J_antes,H_antes,Ganho_antes,erro_antes] = estimacao_estados(G,B,busdata,linedata,Plano_medicao,R);

%% Figuras obtidas após a primeira estimação de estados

figure();
b=bar(1:1:nbus,[V V_est_antes],1.05,'grouped');
max_V = max(V);
max_V_est = max(V_est_antes);
set(b(1),'FaceColor','blue');
set(b(2),'FaceColor',[1 0.8 0]);
set(gca,'FontSize',9);
ylim([0.95 max(max_V,max_V_est)+0.01]);
legend('Real','Estimated');
xlabel('Bus');
xticks(1:nbus);
yticks(linspace(0.95,max(max_V,max_V_est)+0.01,15));
xlim([0 nbus+1]);
ylabel('V (pu)')


figure('Renderer','painters','Position',[100 250 450 250])
er=bar(1:nbus, 100*abs((V-V_est_antes)./V),'k');
erro_pre = 100*abs((V-V_est_antes)./V);
set(gca,'FontSize',9);
ylabel('Relative voltage error (%)');
xlim([0 nbus+1]);
ylim([0 max(erro_pre)+ 0.01]);
xlabel('Bus');
xticks(1:nbus);
yticks(linspace(0,max(erro_pre)+ 0.01,15));

%% Detecção e Identificação de erros grosseiros

J = J_antes;
H = H_antes;
Ganho = Ganho_antes;
erro = erro_antes;

K = H*inv(Ganho)*H'*inv(R); % matriz chapéu
Omega = (eye(size(K,1)) - K)*R; % matriz de covariância residual
Diag_Omega = diag(Omega);
for s=1:size(Diag_Omega,1)
    if Diag_Omega(s,1)<=10^(-8)
        Diag_Omega(s,1) = 10^(-8);
    end
end
rn = abs(erro)./sqrt(Diag_Omega); % resíduos normalizados de cada medição
caso = 1;
suspeitas = [];
%----------------- Caso 1 -----------------------

if caso==1
    % Todos os tipos de medições podem ser levados para correção
    suspeitas = find(rn>=(mean(rn) + std(rn)));

%----------------- Caso 2 -----------------------
elseif caso==2
% Somente Vk,Thk,Pk e Qk podem ser levados para correção
contador = 1;
  for k=1:size(rn,1)
      if (rn(k,1)>=mean(rn) + std(rn)) && Plano_medicao(k,1)~=5 && Plano_medicao(k,1)~=6 && Plano_medicao(k,1)~=7 && Plano_medicao(k,1)~=8
        suspeitas(contador,1) = k;
        contador = contador + 1;
      end
  end

%----------------- Caso 3 -----------------------
elseif caso==3
% Todo plano de medição vai receber correção
    suspeitas = transpose([1:size(Plano_medicao,1)]);
end
 
if size(suspeitas,1)~=0
    %% Correção de erros grosseiros via Algoritmo de Otimização Aritmética
    
    N_simulacoes = 1;

    AOA = 1;

    J_inicial = J;

    data_simu = zeros(N_simulacoes,1);

    for simu=1:N_simulacoes
        AOA = 1;
        %Opção a ser seguida

        % Definindo os parâmetros mi,alpha e epsilon
        mi = 0.1; %constante definida no início do processo
        alpha = 11.0; %constante definida no início do processo
        epsilon = 0.000001; %constante para divisão

        % Valores máximos e mínimos de cada variável inteira
        x_min = 0.1*ones(1,size(suspeitas,1));
        x_max = 10.0*ones(1,size(suspeitas,1));

        % Lista gerada randomicamente para as possíveis soluções inteiras (inicial)
        N_solucoes = 5*size(suspeitas,1);
        N_variaveis = size(suspeitas,1);

        x_inicial = zeros(N_solucoes,N_variaveis);

        for i=1:N_solucoes
            for j=1:N_variaveis
                 x_inicial(i,j) = rand*(x_max(j));
            end
        end

        x_sol = x_inicial;  

        % Vetor solução 
        F_sol = zeros(size(x_sol,1),1);
        F_sol_new = zeros(size(x_sol,1),1);

        % Inicialização da função objetivo (iteração 1)
        for a=1:size(x_sol,1)
              Plano_medicao_new = Plano_medicao;
              for s=1:size(suspeitas,1)
                  Plano_medicao_new(suspeitas(s,1),4) = x_sol(1,s)*Plano_medicao_new(suspeitas(s,1),4);
              end
              [V_est,Th_est,Pk_est,Qk_est,Pkm_est,Qkm_est,Ikm_est,Del_Ikm_est,J,H,Ganho,erro] = estimacao_estados(G,B,busdata,linedata,Plano_medicao_new,R);
              F_sol(a,1) = J;
              if a==1
                Best = F_sol(a,1);
                linha = a;
              end
              if F_sol(a,1)<=Best
                Best = F_sol(a,1);
                linha = a;
              end
        end

        x_inicial_best = x_sol(linha,:);

        % Número de iterações
        N_iter = 200; 

        % Contador de iterações 
        iteracao = 0;

        % Vetores para visualização da convergência
        x_label = zeros(1,N_iter+1);
        y_label = zeros(1,N_iter+1);

        % Matriz de atualização
        x_new = zeros(size(x_sol,1),size(x_sol,2));

        % Definição da função objetivo
        %fprintf('Best (0) = %d \n',Best)

        %-------------------------- OTIMIZAÇÃO ARITMÉTICA -------------------------

        while iteracao<=N_iter

          x_label(iteracao+1) = iteracao;
          y_label(iteracao+1) = Best;
          fprintf('Best (%d) = %d \n',iteracao,Best)

          x_best = (x_sol(linha,:));

          AOM = 0.0 + (iteracao*(1.0 - 0.0)/N_iter);

          OMP = 1 - ((iteracao^(1/alpha))/(N_iter^(1/alpha)));

          for a=1:size(x_sol,1)

              for b=1:size(x_sol,2)

                r1 = rand(); %sorteando um valor para cada solução
                r2 = rand(); %sorteando um valor para cada solução
                r3 = rand(); %sorteando um valor para cada solução

                if r1<AOM

                   % ------------------- Busca Local ------------------------------
                  if r3<0.5
                    x_new(a,b) = (((x_best(b) - OMP*((x_max(b) - x_min(b))*mi + x_min(b)))));
                  else
                    x_new(a,b) = (((x_best(b) + OMP*((x_max(b) - x_min(b))*mi + x_min(b)))));
                  end

                else

                  % --------------------- Busca Global ----------------------------
                  if r2<0.5
                    x_new(a,b) = (((x_best(b)/(OMP + epsilon)*((x_max(b) - x_min(b))*mi + x_min(b)))));
                  else
                    x_new(a,b) = (((x_best(b)*(((OMP)*((x_max(b) - x_min(b))*mi + x_min(b)))))));
                  end
                end

                if x_new(a,b)>x_max(b)
                   x_new(a,b) = x_best(b);
                end
                if x_new(a,b)<x_min(b)
                   x_new(a,b) = 1.0;
                end


              end

              Plano_medicao_new = Plano_medicao;
              for s=1:size(suspeitas,1)
                  Plano_medicao_new(suspeitas(s,1),4) = x_new(1,s)*Plano_medicao_new(suspeitas(s,1),4);
              end
              [V_est,Th_est,Pk_est,Qk_est,Pkm_est,Qkm_est,Ikm_est,Del_Ikm_est,J,H,Ganho,erro] = estimacao_estados(G,B,busdata,linedata,Plano_medicao_new,R);
              F_sol_new(a,1) = J;

              if F_sol_new(a,1)<F_sol(a,1)
                  F_sol(a,1) = F_sol_new(a,1);
                  x_sol(a,:) = x_new(a,:);
              end
          end

          for k=1:size(x_sol,1)
              if F_sol(k,1)< Best
                  Best = F_sol(k,1);
                  linha = k;
              end
          end

          iteracao = iteracao + 1;
        end
        data_simu(1,simu) = Best;
        fprintf('Best - Simulação (%d) = %d \n',simu,Best)
    end

    Fatores_otimos = x_best;

    Plano_medicao_corrigido = Plano_medicao;
    for s=1:size(suspeitas,1)
      Plano_medicao_corrigido(suspeitas(s,1),4) = Fatores_otimos(1,s)*Plano_medicao_corrigido(suspeitas(s,1),4);
    end
    [V_est,Th_est,Pk_est,Qk_est,Pkm_est,Qkm_est,Ikm_est,Del_Ikm_est,J,H,Ganho,erro] = estimacao_estados(G,B,busdata,linedata,Plano_medicao_corrigido,R);
    clc;
    toc

    %% FIGURAS obtidas DEPOIS do(s) erro(s) grosseiro(s) tratado (s) usando AOA

    figure();
    c=bar(1:1:nbus,[V V_est],1.05,'grouped');
    max_V = max(V);
    max_V_est = max(V_est);
    set(c(1),'FaceColor','blue');
    set(c(2),'FaceColor',[1 0.8 0]);
    set(gca,'FontSize',9);
    ylim([0.95 max(max_V,max_V_est)+0.01]);
    legend('Real','Estimated');
    xlabel('Barra');
    xticks(1:nbus);
    yticks(linspace(0.95,max(max_V,max_V_est)+0.01,15));
    xlim([0 nbus+1]);
    ylabel('Voltage (pu)')


    figure('Renderer','painters','Position',[100 250 450 250])
    er=bar(1:nbus, 100*abs((V-V_est)./V),'k');
    erro_pos = 100*abs((V-V_est)./V);
    set(gca,'FontSize',9);
    ylabel('Relative voltage error (%)');
    xlim([0 nbus+1]);
    ylim([0 max(erro_pos)+ 0.01]);
    xlabel('Bus');
    xticks(1:nbus);
    yticks(linspace(0,max(erro_pos)+ 0.01,15));

    figure('Renderer','painters','Position',[100 250 450 250])
    er2=bar([1:nbus], [erro_pre,erro_pos],1.0,'grouped');
    set(er2(1),'FaceColor','red');
    set(er2(2),'FaceColor','green');
    set(gca,'FontSize',9)
    ylabel('Relative voltage error (%)')
    xlim([0 nbus+1])
    ylim([0 max(erro_pre)+ 0.1])
    %title('Erro relativo da tensão na estimação de estados')
    legend('With bad data','Bad data corrected')
    xlabel('Bus')
    xticks([1:nbus])
    yticks([linspace(0,max(erro_pre)+ 0.1,15)])


    %% Visualização da convergência
    figure();
    plot(x_label,y_label,'Color','r','linewidth',2.5)
    xlabel('Iteration')
    ylabel('Objective Function Value - J(x)')
    %legend({'Processo de convergência'})
    %title('Algoritmo de Otimização Aritmética (AOA) para correção de erros grosseiros','FontSize',9)
    xlim([1 N_iter])
    ylim([0 max(y_label)+1])
    yticks(linspace(0,max(y_label)+1,15))
    grid on

else
    
    print('Estimação de estados sem dados espúrios detectados foi finalizada com sucesso');
    
end

