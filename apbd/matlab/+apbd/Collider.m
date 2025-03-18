classdef Collider < handle
	%Collider Collision handler
	%   Calls broadphase and narrowphase

	properties
		model
        bpList1
        bpList2
        groundCollNum
        bodyCollNum
		collisions
        activeCollisions
        bodyNum
	end

	methods
		%%
        function this = Collider(model, groundCollisionList, bodyCollisionList)
			this.model = model;
            this.groundCollNum = length(groundCollisionList);
            this.bodyCollNum = size(bodyCollisionList,1);
            this.bodyNum = length(model.bodies);
            this.activeCollisions = [];
			this.collisions = cell(1,this.groundCollNum + this.bodyCollNum);
            groundBody = apbd.BodyRigid(apbd.ShapeCuboid([1 1 0.1]),Inf);
            groundBody.layer = 0;
            groundBody.index = 0;
            for i = 1 : this.groundCollNum
                this.collisions{i} = Collision(model.bodies{groundCollisionList(i)}, groundBody, true);
            end
            for i = 1 : this.bodyCollNum
                this.collisions{i+this.groundCollNum} = Collision(model.bodies{bodyCollisionList(i,1)}, model.bodies{bodyCollisionList(i,2)}, false);
            end
		end

		%%
		function run(this)
            this.activeCollisions = [];
            this.bpList1 = [];
            this.bpList2 = [];
			this.broadphase();
			this.narrowphase();
        end

		%%
		function broadphase(this)
			% Body-ground collisions
			for i = 1 : this.groundCollNum
				body = this.collisions{i}.body1;
				if body.collide
					if body.broadphaseGround(this.model.ground.E)
						this.bpList1(end+1) = i;
					end
				end
			end

			% Body-body collisions
			for i = this.groundCollNum +1 : this.groundCollNum + this.bodyCollNum
                body1 = this.collisions{i}.body1;
                body2 = this.collisions{i}.body2;
				if body1.collide && body2.collide
                    if body1.broadphaseRigid(body2)
	                    this.bpList2(end+1) = i;
                    end
				end
			end
		end

		%%
		function narrowphase(this)

			% Body-ground collisions
			Eg = this.model.ground.E;
            for i = this.bpList1
				body = this.collisions{i}.body1;
			    cdata = body.narrowphaseGround(Eg);
                if(~isempty(cdata))
                    this.collisions{i}.setContacts(cdata);
                    this.collisions{i}.getConstraints();
                    this.activeCollisions(end+1) = i;
                end
            end

			% Body-body collisions
            for i = this.bpList2
                body1 = this.collisions{i}.body1;
                body2 = this.collisions{i}.body2;
			    cdata = body1.narrowphaseRigid(body2);
                if(~isempty(cdata))
                    this.collisions{i}.setContacts(cdata);
                    this.collisions{i}.getConstraints();
                    this.activeCollisions(end+1) = i;
                end
            end
        end

		%%
		function flag = isempty(this)
			flag = isempty(this.collisions);
		end
	end
end
