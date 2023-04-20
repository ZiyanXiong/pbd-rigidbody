classdef ConstraintPosition < Constraint
    %CONSTRAINTPOSITION Positional Constraints
    
    properties
        n1    %In local frame
        n2    %In local frame
        c
        r1    %In local frame
        r2    %In local frame
    end
    
    methods
        function this = ConstraintPosition(bodies_)
            this = this@Constraint(bodies_);
        end
        
        function solvePositions_(this)
            nbodies = length(this.bodies);
			if nbodies == 1
                body = this.bodies{1};
                temp = cross(this.r1, this.n1);
                w1 = 1 / body.mass + temp' * body.I_inv * temp;
                delta_lambda = -this.c / w1;
                p = delta_lambda * this.n1;
                body.x = body.x + p / body.mass;
                temp = body.I_inv * cross(this.r1, p);
                body.q = body.q + 0.5 * quaternion([0 temp']) * body.q;
            elseif nbodies == 2
                body1 = this.bodies{1};
                body2 = this.bodies{2};
                temp = cross(this.r1, this.n1);
                w1 = 1 / body1.mass + temp' * body1.I_inv * temp;
                temp = cross(this.r2, this.n2);
                w2 = 1 / body2.mass + temp' * body2.I_inv * temp;

                delta_lambda = -this.c / (w1 + w2);

                p = delta_lambda * this.n1;
                body1.x = body1.x + p / body1.mass;
                temp = body1.I_inv * cross(this.r1, p);
                body1.q = body1.q + 0.5 * quaternion([0 temp']) * body1.q;

                p = delta_lambda * this.n2;
                body2.x = body2.x + p / body2.mass;
                temp = body2.I_inv * cross(this.r2, p);
                body2.q = body2.q + 0.5 * quaternion([0 temp']) * body2.q;
            else
                fprintf("Invalid body number in ConstraintPosition\n");
            end
        end
    end
end

