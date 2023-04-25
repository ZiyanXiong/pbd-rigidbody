classdef (Abstract) Constraint < handle
    %CONSTRAINTPOSITION Super class for Constraints
    
    properties
        bodies
        C
        dC
    end

    %%
    methods
        %%
        function this = Constraint(bodies_)
            this.bodies = bodies_;
        end

        %%
        function solvePositions(this)
            nbodies = length(this.bodies);
            if nbodies == 1
                body = this.bodies{1};
		        lambda = -this.C/(this.dC*(this.dC'./body.M));
		        body.q = body.q + lambda*(this.dC'./body.M);
            end
        end
    end
end

