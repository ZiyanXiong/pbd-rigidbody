classdef ConstraintGroundContact < ConstraintPosition
    %CONSTRAINTGROUNDCONTACT Constraint for ground contact
    
    properties
        collisions
        E_g
        planar % 2D planar scene
    end
    
    methods
        function this = ConstraintGroundContact(body_, E_g_)
            this = this@ConstraintPosition(body_);
            this.E_g = E_g_;
            this.planar = false;
        end
        
        function collide(this)
			% Clear collisions
			this.collisions = {};
			
			% Ground collisions (only if ground.E isn't the zero matrix)
			if this.E_g(4,4) == 1
				this.collisions = this.bodies{1}.collideGround(this.E_g, this.planar, this.collisions);
            end
        end

        function solvePositions_(this)
            this.collide();
            for i = 1 : length(this.collisions)
				collision = this.collisions{i};
                body = collision.bodies{1};
                this.r1 = collision.xl{1};
                this.n1 = body.E_iw(1:3,1:3) * collision.nw;
                this.c = collision.d;
                solvePositions_@ConstraintPosition(this);
            end
        end
    end
end

