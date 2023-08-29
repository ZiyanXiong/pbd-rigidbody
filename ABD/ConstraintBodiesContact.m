classdef ConstraintBodiesContact < Constraint
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
        collisions
        planar % 2D planar scene
    end
    
    methods
        function this = ConstraintBodiesContact(bodies_)
            this = this@Constraint(bodies_);
            this.planar = false;
        end
        
        function collide(this)
			% Clear collisions
			this.collisions = {};
			nbodies = length(this.bodies);
			% Bodies collisions (only if there are only 2 bodies)
			if nbodies == 2
				this.collisions = this.bodies{1}.collideBody(this.bodies{2}, this.planar, this.collisions);
            end
        end

        function solvePositions(this)
            this.collide();
            nbodies = length(this.bodies);
            if nbodies == 2
                for i = 1 : length(this.collisions)
				    c = this.collisions{i};
		            %this.C = c.nw'*c.xw;
		            %this.dC = c.nw'* se3.getJacobian(c.xl{1});
                    %this.dC2 = -c.nw'* se3.getJacobian(c.xl{2});
		            this.C = c.d;
		            this.dC = c.nw'* se3.getJacobian(c.xl{1});
                    this.dC2 = c.nw'* se3.getJacobian(c.xl{2});
                    solvePositions@Constraint(this);
                end
            end
        end
    end
end

