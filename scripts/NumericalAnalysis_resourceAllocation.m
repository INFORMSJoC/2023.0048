% more streamlined version. updates done more efficiently.

tic
N = 25;
m = 5;
BT = 45;
% k = 9;
k_max = 10;
iters = 1000;
% eps_max = .5;
% eps_steps = 200;



get_index = @(head,tail) (head-1)*N - (head-1)*head/2 + tail - head;
% eps_vals = linspace(eps_max,0,eps_steps);
n_edges = N*(N-1)/2;
s = zeros(1,n_edges);
t = zeros(1,n_edges);
i1 = 1;
i2 = N-1;

for i=1:(N-1)
    s(i1:i2) = i;
    t(i1:i2) = (i+1):N;
    i1 = i1 + N - i;
    i2 = i1 + N - 2 - i;
end
n = max(t);

results = cell(k_max,iters);
% fprintf('Progress:\n');
% fprintf(['\n' repmat('.',1,iters) '\n\n']);
parfor_progress(iters);
parfor itr = 1:iters
%     fprintf('\b|\n');
    pause(1e-10);
    parfor_progress;
    K = 3*rand(m,1)';
    prev = cell(m,2);
    for i = 1:m
        [prev{i,1}, prev{i,2}] = generateRandomPrevalenceFunctions(N);
    end

    weights = zeros(1,n_edges);

    W = zeros(N,N)-1;

    for i = 1:(N-1)

        p_min = zeros(m,1);
        p_max = zeros(m,1);

        for j = 1:m
            p_min(j) = min(prev{j,1}(i:(i+1)));
            p_max(j) = max(prev{j,1}(i:(i+1)));
        end

        MUs = [p_min, p_max];
        [r,~] = minMaxRegretV2(MUs,BT,K);
        W(i,i+1) = max(r);
    end

    for i = 1:n_edges
        i_start = s(i);
        i_end = t(i);

        mmr = 0;
        for j = i_start:(i_end-1)
            mmr = max(mmr,W(j,j+1));
        end
        weights(i) = mmr;
    end

    W_init = W;
    weights_init = weights;


    % E = [];

    % results = zeros(k_max,eps_steps);


    % algoProgress = {};

    for k = 1:k_max

        W = W_init;
        weights = weights_init;
        %     counter = 1;
        error = inf;
        eps = inf;
        progress = [];

        while error > 0

            [p,d] = BellmanFord_V2_maxSpread(s,t,weights,k,1,n);
            %     algoProgress = [algoProgress;p];

            real_d = 0;

            for i = 1:(length(p)-1)
                %         i
                if W(p(i),p(i+1)) ~= -1 %if already calculated then update width and SKIP
                    real_d = max(real_d,W(p(i),p(i+1)));
                    continue
                end

                p_min = zeros(m,1);
                p_max = zeros(m,1);

                for j = 1:m
                    p_min(j) = min(prev{j,1}(p(i):p(i+1)));
                    p_max(j) = max(prev{j,1}(p(i):p(i+1)));
                end



                MUs = [p_min, p_max ];

                [r,~] = minMaxRegretV2(MUs,BT,K);
                mmr = max(r);


                for j = 1:p(i)
                    for l = p(i+1):N

                        if W(j,l) ~= -1
                            continue
                        end

                        id = get_index(j,l);
                        r_max = max([mmr,weights(id)]);
                        weights(id)=r_max;

                    end
                end
                W(p(i),p(i+1)) = mmr;
                real_d = max(real_d,mmr);

            end

            error = (real_d-d);%/real_d; %CHECK!!!

            if error < eps
                progress = [real_d, d, sum(sum(W>-1))/n_edges; progress];
                eps = error;
            end


            %         if error < eps_vals(counter)
            %
            %             results(k, counter) = sum(sum(W>-1))/n_edges;
            %             counter = counter + 1;
            %         end


            %         E = [E ; real_d, d];


        end
        %     d_end = progress(1,1)
        %     progress(:,1) = (progress(:,1)-d_end) / d_end;
        results{k,itr} = progress;


    end

    % p
    % d
end
parfor_progress(0);
toc
% save("results25.mat","results")



algoProgress = cell(1,10,1000);

Nlist = [25];

for k = 1:10
    for s = 1:1000
        algoProgress{1,k,s} = [(results{k,s}(:,1)-results{k,s}(:,2))./results{k,s}(:,1), results{k,s}(:,3)];


        for N = 1:1
            algoProgress{N,k,s} = [algoProgress{N,k,s}; 1, 2/(Nlist(N)-1)];
            [~,ia] = unique(algoProgress{N,k,s}(:,1));
            algoProgress{N,k,s} = algoProgress{N,k,s}(ia,:);
        end


    end
end

x = linspace(0,1,100);

y = cell(1,10,1000);



for N = 1:1
    for k = 1:10
        for s = 1:1000
            if size(algoProgress{N,k,s},1)<2
                continue
            end
            y{N,k,s} = interp1(algoProgress{N,k,s}(:,1),algoProgress{N,k,s}(:,2),x);

        end
    end
end

close all
Yexport = cell(1,1);

for N = 1:1
    figure;hold on;
    Y = 0;
    Ytmp = zeros(length(x),10);
    for k = 1:10
        for i = 1:1000
            Y = Y + y{N,k,i};
        end
        Y = Y/1000;
        Ytmp(:,k) = Y;
        plot(x,Y)
    end
    Yexport{N} = Ytmp;
    leg = cell(10,1);
    for k = 1:10
        leg{k} = ['k = ', num2str(k)];
    end
    legend(leg)
end