function [mat,vars,c_idx,rm_idx] = generateInputMat(input_mat,input_vars,input_ind,w,tp)
    mat = input_mat(input_ind==1); mat = [mat{:}];
    vars = input_vars(input_ind==1); vars = [vars{:}];
    c_idx = []; rm_idx = [];
    if w>0
        mat = mat(:,linspace(w,(size(mat,2)+w-tp),(size(mat,2)/tp)));
    elseif w == 0
        z_idx = mat~=0;
        c = sum(z_idx,1);
        rm_idx = find(c<3);
        c_idx = rm_idx;
    end
%     if w > 0

%         z_idx = mat~=0;
%         c = sum(z_idx,1);
%         c_idx = find(c<3);
%         mat(:,c_idx) = [];
%         vars(:,c_idx) = [];
%         rm_idx = [];
%     elseif ~isempty(w)

%     end
end