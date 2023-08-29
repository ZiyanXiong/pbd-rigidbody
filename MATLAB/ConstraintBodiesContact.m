classdef ConstraintBodiesContact < ConstraintPosition
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
        collisions
        planar % 2D planar scene
    end
    
    methods
        function this = ConstraintBodiesContact(bodies_)
            this = this@ConstraintPosition(bodies_);
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

        function solvePositions_(this)
            this.collide();
            nbodies = length(this.bodies);
            body1 = this.bodies{1};
            body2 = this.bodies{2};
            if nbodies == 2
                for i = 1 : length(this.collisions)
				    collision = this.collisions{i};
                    this.r1 = collision.xl{1};
                    this.n1 = body1.E_iw(1:3,1:3) * (collision.nw);
                    this.c = collision.d;
                    this.r2 = collision.xl{2};
                    this.n2 = body2.E_iw(1:3,1:3) * (-collision.nw);
                    solvePositions_@ConstraintPosition(this);
                end
            end
        end
    end
end

