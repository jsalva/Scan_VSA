function [output dP_vec] = interpz(dP,error_vec,dproj_vec)
for i = 1:length(dproj_vec)
   proj_vec(i) = sum(dproj_vec(1:i)); 
end

inchworm = 0;
ideal_traj_full_dist = sum(abs(dproj_vec));
output = zeros(ceil(ideal_traj_full_dist / dP),1);
dP_vec = zeros(ceil(ideal_traj_full_dist / dP),1);

E1 = 0;
E2 = 0;
P1 = 0;
P2 = 0;

for i = 1:length(output)
    inchworm = inchworm + dP;
    idx1 = find(proj_vec <= inchworm);
    idx2 = find(proj_vec >= inchworm);
    if isempty(idx1)
        E1 = 0;
        E2 = error_vec(1);
        P1 = 0;
        P2 = proj_vec(1);
        output(i) = E1 + (E2-E1)/(P2-P1)*(inchworm-P1);
        dP_vec(i) = inchworm;
    elseif isempty(idx2)
        output(i) = error_vec(end);
        dP_vec(i) = inchworm;
    else
       E1 = error_vec(idx1(end));
       E2 = error_vec(idx2(1));
       P1 = proj_vec(idx1(end));
       P2 = proj_vec(idx2(1));
       output(i) = E1 + (E2-E1)/(P2-P1)*(inchworm-P1);
       dP_vec(i) = inchworm;

    end
    
    
end
disp('error:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(length(output))
disp('dP:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(length(dP_vec))



end