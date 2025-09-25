for i=1:N
if i == choos
            % 攻击节点: 固定异常值
            for j = 1:numF
                c_p(i,j,t) = -5;
                c_v(i,j,t) = 0;
            end
            p(i,t) = -5;
            vn(i,t) = 0;
else
            % 对每个子集 S 更新
            for j = 1:numF
                S = F{j};
                Nn = current_neighbors{i}; % 使用当前拓扑的邻居
                % 去掉 S 中的节点
                N_j = setdiff(Nn, S);
                % 计算控制输入
                % sum_p_diff = 0;
                % sum_v_diff = 0;
                for k = N_j
                    sum_p_diff(i,j,t-1) = sum_p_diff(i,j,t-1) + (c_p(i,j,t-1) - c_p(k,j,t-1));%
                    % sum_v_diff = sum_v_diff + (c_v(u,j,t) - c_v(v,j,t));
                end
                u_val(i,j,t) = -0.1 * sum_p_diff(i,j,t-1) - 1 * c_v(i,j,t-1);%sum_v_diff
                % 更新速度和位置
                c_v(i,j,t) = c_v(i,j,t-1) + u_val(i,j,t);
                c_p(i,j,t) = c_p(i,j,t-1) + c_v(i,j,t-1);   % 使用 t 时刻的速度
            end

            % Step 3: 选择子集
            chosen = 1; % 默认 ∅
            base_p = c_p(i,1,t);
            for j2 = 2:numF
                if abs(base_p - c_p(i,j2,t)) > epsilon
                    chosen = j2;
                    break;
                end
            end
            p(i,t) = c_p(i,chosen,t);
            vn(i,t) = c_v(i,chosen,t);
end
end