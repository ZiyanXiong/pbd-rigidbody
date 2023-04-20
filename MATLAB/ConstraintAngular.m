classdef ConstraintAngular < Constraint
    %CONSTRAINTPOSITION Positional Constraints
    
    properties
        n1    %In local frame
        n2    %In local frame
        theta
    end
    
    methods
        function this = ConstraintAngular(bodies_)
            this = this@Constraint(bodies_);
        end
        
        function solvePositions_(this)
            nbodies = length(this.bodies);
			if nbodies == 1
                body1 = this.bodies{1};
                w1 = this.n1' * body1.I_inv * this.n1;

                delta_lambda = -this.theta / w1;

                p = delta_lambda * this.n1;
                temp = body1.I_inv * p;
                body1.q = body1.q + 0.5 * quaternion([0 temp']) * body1.q;
            elseif nbodies == 2
                body1 = this.bodies{1};
                body2 = this.bodies{2};
                w1 = this.n1' * body1.I_inv * this.n1;
                w2 = this.n2' * body2.I_inv * this.n2;

                delta_lambda = -this.theta / (w1 + w2);

                p = delta_lambda * this.n1;
                temp = body1.I_inv * p;
                body1.q = body1.q + 0.5 * quaternion([0 temp']) * body1.q;

                p = delta_lambda * this.n2;
                temp = body2.I_inv * p;
                body2.q = body2.q - 0.5 * quaternion([0 temp']) * body2.q;
            else
                fprintf("Invalid body number in ConstraintPosition\n");
            end
        end
    end
end

