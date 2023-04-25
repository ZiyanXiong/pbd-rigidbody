classdef ConstraintSphericalJoint < ConstraintPosition
    %CONSTRAINTGROUNDCONTACT Constraint for spherical joint
    
    properties
        x_joint_w %world position for the joint
        x_joint_l1 %local position for the joint in body1's frame
        x_joint_l2 %local position for the joint in body1's frame

    end
    
    methods
        function this = ConstraintSphericalJoint(bodies_, x_joint_local)
            this = this@ConstraintPosition(bodies_);
            this.x_joint_l1 = x_joint_local;
            this.x_joint_w = bodies_{1}.E_wi * [x_joint_local; 1];
            this.x_joint_w = this.x_joint_w(1:3);
        end
        
        function solvePositions_(this)
            nbodies = length(this.bodies);
			if nbodies == 1
                body = this.bodies{1};
                x_joint_w_current = body.E_wi * [this.x_joint_l1; 1];
                x_joint_w_current = x_joint_w_current(1:3);
                delta_x = x_joint_w_current - this.x_joint_w;
                this.r1 = this.x_joint_l1;
                this.n1 = body.E_iw(1:3,1:3) * delta_x ./ norm(delta_x);
                this.c = norm(delta_x);
                if this.c > 1e-12
                    solvePositions_@ConstraintPosition(this);
                end
            end
        end
    end
end