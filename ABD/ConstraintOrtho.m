classdef ConstraintOrtho < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
    end
    
    methods
        %%
        function this = ConstraintOrtho(body_)
            this = this@Constraint(body_);
        end

        %%
        function solvePositions(this)
            body = this.bodies{1};
            A = reshape(body.q(4:12),3,3)';
            C(1,1) = A(1,1)^2 + A(2,1)^2 + A(3,1)^2 - 1;
            C(2,1) = A(1,2)^2 + A(2,2)^2 + A(3,2)^2 - 1;
            C(3,1) = A(1,3)^2 + A(2,3)^2 + A(3,3)^2 - 1;
            C(4,1) = A(1,1)*A(1,2) + A(2,1)*A(2,2) + A(3,1)*A(3,2);
            C(5,1) = A(1,2)*A(1,3) + A(2,2)*A(2,3) + A(3,2)*A(3,3);
            C(6,1) = A(1,1)*A(1,3) + A(2,1)*A(2,3) + A(3,1)*A(3,3); 
            dC = [
                0, 0, 0, 2*A(1,1),        0,        0, 2*A(2,1),        0,        0, 2*A(3,1),        0,        0
                0, 0, 0,        0, 2*A(1,2),        0,        0, 2*A(2,2),        0,        0, 2*A(3,2),        0
                0, 0, 0,        0,        0, 2*A(1,3),        0,        0, 2*A(2,3),        0,        0, 2*A(3,3)
                0, 0, 0,   A(1,2),   A(1,1),        0,   A(2,2),   A(2,1),        0,   A(3,2),   A(3,1),        0
                0, 0, 0,        0,   A(1,3),   A(1,2),        0,   A(2,3),   A(2,2),        0,   A(3,3),   A(3,2)
                0, 0, 0,   A(1,3),        0,   A(1,1),   A(2,3),        0,   A(2,1),   A(3,3),        0,   A(3,1)
            ];
            for i = 1 : length(C)
                this.C = C(i);
                this.dC = dC(i,:);
                solvePositions@Constraint(this);
            end
            this.bodies{1}.updateE();
        end
    end
end

