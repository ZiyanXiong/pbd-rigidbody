classdef (Abstract) Constraint < handle
    %CONSTRAINTPOSITION Super class for Constraints
    
    properties
        bodies
    end
    
    methods
        function this = Constraint(bodies_)
            this.bodies = bodies_;
        end
        
        function solvePositions(this)
            this.solvePositions_();
        end
    end
end

