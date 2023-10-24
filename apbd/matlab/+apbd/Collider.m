classdef Collider < handle
	%Collider Collision handler
	%   Calls broadphase and narrowphase

	properties
		model
		bpList1
		bpList2
		collisions
	end

	methods
		%%
		function this = Collider(model)
			this.model = model;
			this.bpList1 = {};
			this.bpList2 = {};
			this.collisions = {};
		end

		%%
		function run(this)
			this.bpList1 = {};
			this.bpList2 = {};
			this.collisions = {};
			this.broadphase();
			this.narrowphase();
		end

		%%
		function broadphase(this)
			bodies = this.model.bodies;

			% Body-ground collisions
			for i = 1 : length(bodies)
				body = bodies{i};
				if body.collide
					if body.broadphaseGround(this.model.ground.E)
						this.bpList1{end+1} = body;
					end
				end
			end

			% Body-body collisions
			for i = 1 : length(bodies)
				if bodies{i}.collide
					for j = i+1 : length(bodies)
						if bodies{j}.collide
							if bodies{i}.broadphaseRigid(bodies{j})
								this.bpList2{end+1} = {bodies{i},bodies{j}};
							end
						end
					end
				end
			end
		end

		%%
		function narrowphase(this)
			% Body-ground collisions
			Eg = this.model.ground.E;
			for i = 1 : length(this.bpList1)
				body = this.bpList1{i};
				cdata = body.narrowphaseGround(Eg);
				for k = 1 : length(cdata)
					c = cdata(k);
					if isa(body,'apbd.BodyRigid')
						this.collisions{end+1} = apbd.ConCollGroundRigid(body,Eg);
					elseif isa(body,'apbd.BodyRigid2d')
						this.collisions{end+1} = apbd.ConCollGroundRigid2d(body,Eg);
					elseif isa(body,'apbd.BodyAffine')
						this.collisions{end+1} = apbd.ConCollGroundAffineVal(body,Eg);
					end
					this.collisions{end}.d = c.d;
					this.collisions{end}.xl = c.xl;
					this.collisions{end}.xw = c.xw;
					this.collisions{end}.nw = c.nw;
					this.collisions{end}.vw = c.vw;
				end
			end

			% Body-body collisions
			for i = 1 : length(this.bpList2)
				body1 = this.bpList2{i}{1};
				body2 = this.bpList2{i}{2};
				cdata = body1.narrowphaseRigid(body2);
				for k = 1 : length(cdata)
					c = cdata(k);
					if isa(body1,'apbd.BodyRigid') && isa(body2,'apbd.BodyRigid')
							this.collisions{end+1} = apbd.ConCollRigidRigid(body1,body2);
                    end
					if isa(body1,'apbd.BodyRigid2d') && isa(body2,'apbd.BodyRigid2d')
					    this.collisions{end+1} = apbd.ConCollRigidRigid2d(body1,body2);
                    end
					if isa(body1,'apbd.BodyAffine') && isa(body2,'apbd.BodyAffine')
					    this.collisions{end+1} = apbd.ConCollAffineAffineVal(body1,body2);
					end
					this.collisions{end}.setData(c);
				end
			end
		end

		%%
		function flag = isempty(this)
			flag = isempty(this.collisions);
		end
	end
end
