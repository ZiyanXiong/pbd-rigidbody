classdef ConstraintSphericalJoint < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for spherical joint
    
    properties
        x_joint_w %world position for the joint
        x_joint_l1 %local position for the joint in body1's frame
        x_joint_l2 %local position for the joint in body1's frame

    end
    
    methods
        function this = ConstraintSphericalJoint(bodies_, x_joint_local)
            this = this@Constraint(bodies_);
            nbodies = length(this.bodies);
			if nbodies == 1
                this.x_joint_l1 = x_joint_local;
                this.x_joint_w = bodies_{1}.E_wi * [x_joint_local; 1];
                this.x_joint_w = this.x_joint_w(1:3);
            elseif nbodies == 2
                this.x_joint_l1 = x_joint_local;
                this.x_joint_w = bodies_{1}.E_wi * [x_joint_local; 1];
                this.x_joint_w(4) = 1;
                this.x_joint_l2 = se3.inv(bodies_{2}.E_wi) * this.x_joint_w;
                this.x_joint_l2 = this.x_joint_l2(1:3);
                this.x_joint_w = this.x_joint_w(1:3);
            end
        end
        
        function solvePositions(this)
            nbodies = length(this.bodies);
			if nbodies == 1
                body = this.bodies{1};
                x_joint_l1_w = body.E_wi * [this.x_joint_l1; 1];
                x_joint_l1_w = x_joint_l1_w(1:3);
                delta_x = x_joint_l1_w - this.x_joint_w;
                this.C = norm(delta_x);
                delta_x = delta_x ./ this.C;
                this.dC = delta_x'* se3.getJacobian(this.x_joint_l1);
                if this.C > 1e-12
                    solvePositions@Constraint(this);
                end
            elseif nbodies == 2
                body1 = this.bodies{1};
                x_joint_l1_w = body1.E_wi * [this.x_joint_l1; 1];
                x_joint_l1_w = x_joint_l1_w(1:3);
                body2 = this.bodies{2};
                x_joint_l2_w = body2.E_wi * [this.x_joint_l2; 1];
                x_joint_l2_w = x_joint_l2_w(1:3);
                delta_x = x_joint_l1_w - x_joint_l2_w;
                this.C = norm(delta_x);
                delta_x = delta_x ./ this.C;
                this.dC = delta_x'* se3.getJacobian(this.x_joint_l1);
                this.dC2 = delta_x'* se3.getJacobian(this.x_joint_l2);
                if this.C > 1e-12
                    solvePositions@Constraint(this);
                end
            end
        end
    end
end