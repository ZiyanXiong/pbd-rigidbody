classdef (Abstract) ConBase < handle
    %Con Superclass for Constraints
    
    properties
		m       % Number of constraint rows
        C       % Constraint value
        lambda  % Lagrange multiplier
        ground  % If this constraint contains ground
    end

    %%
    methods
        %%
		function this = ConBase(m)
			this.m = m;
			this.C = zeros(m,1);
			this.lambda = zeros(m,1);
            this.ground = false;
		end

		%%
		function clear(this)
			this.C = zeros(this.m,1);
			this.lambda = zeros(this.m,1);
		end

		%%
		function V = computeEnergy(this) %#ok<MANU>
			V = 0;
		end
	end

	%%
	methods (Abstract)
		%%
		init(this);

		%%
		solve(this,h);

		%%
		draw(this);
	end
end