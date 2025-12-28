classdef RunningMeanStd < handle
    % 描述:
    % 使用 Welford's online algorithm 实现移动平均值和方差的计算。
    % https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
    
    properties
        Mean
        Variance
        Count
        Shape
    end
    
    methods
        function this = RunningMeanStd(shape)
            if nargin > 0 && ~isempty(shape)
                this.Shape = shape;
                this.Mean = zeros(shape);
                this.Variance = ones(shape);
                this.Count = 1e-4; % 避免除零
            else
                this.Shape = [];
                this.Mean = 0;
                this.Variance = 1;
                this.Count = 1e-4;
            end
        end
        
        function update(this, x)
            batch_mean = mean(x, ndims(x));
            batch_var = var(x, 0, ndims(x));
            batch_count = size(x, ndims(x));
            
            this.update_from_moments(batch_mean, batch_var, batch_count);
        end
        
        function update_from_moments(this, batch_mean, batch_var, batch_count)
            delta = batch_mean - this.Mean;
            tot_count = this.Count + batch_count;
            
            new_mean = this.Mean + delta .* batch_count ./ tot_count;
            m_a = this.Variance .* this.Count;
            m_b = batch_var .* batch_count;
            M2 = m_a + m_b + (delta.^2) .* this.Count .* batch_count ./ tot_count;
            new_var = M2 ./ tot_count;
            
            this.Mean = new_mean;
            this.Variance = new_var;
            this.Count = tot_count;
        end
    end
end