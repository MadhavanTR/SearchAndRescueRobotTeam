function L_mat= findLaplacian(N, followers, leader_agent)
    L = zeros(N, N);
    if(leader_agent==1)
        L(2:N, 2:N) = followers;
        L(2, 2) = L(2, 2) + 1;  % 2nd alone is connected to da leader
        L(2, 1) = -1;
    elseif(leader_agent==N)
        L(1:N-1, 1:N-1) = followers;
        L(N-1, N-1) = L(N-1, N-1) + 1;
        L(N-1, N) = -1;
    else
        L(1:leader_agent-1, 1:leader_agent-1) = followers(1:leader_agent-1, 1:leader_agent-1);
        L(leader_agent, :) = zeros(1, N);
        L(leader_agent+1, leader_agent) = -1;
        L(leader_agent+1:N, leader_agent+1:N) = followers(leader_agent:N-1, leader_agent:N-1);
        L(1:leader_agent-1, leader_agent+1:N) = followers(1:leader_agent-1, leader_agent:N-1);
        L(leader_agent+1:N, 1:leader_agent-1) = followers(leader_agent:N-1, 1:leader_agent-1);
        L(leader_agent+1, leader_agent+1)=L(leader_agent+1, leader_agent+1)+1;
    end
    L_mat=L;
end