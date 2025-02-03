classdef Collider < handle
	%Collider Collision handler
	%   Calls broadphase and narrowphase

	properties
		model
		bpList1
		bpList2
		collisions
        activeCollisions
        groundBodyIndex
        bodyNum
	end

	methods
		%%
		function this = Collider(model)
			this.model = model;
			this.bpList1 = {};
			this.bpList2 = {};
            this.bodyNum = length(model.bodies);
            this.activeCollisions = [];
            this.groundBodyIndex = [];
			this.collisions = cell(1,(this.bodyNum+1)*this.bodyNum/2);
            groundBody = apbd.BodyRigid(apbd.ShapeCuboid([1 1 0.1]),Inf);
            groundBody.layer = 0;
            groundBody.index = 0;
            for i = 1 : this.bodyNum
                for j = i : this.bodyNum
                    index = (2*this.bodyNum+2-i)*(i-1)/2 + j - i + 1;
                    %disp(index);
                    if i == 1
                        this.collisions{index} = Collision(model.bodies{j}, groundBody, true);
                    else
                        this.collisions{index} = Collision(model.bodies{i-1}, model.bodies{j}, false);
                    end
                end
            end
		end

		%%
		function run(this)
			this.bpList1 = {};
			this.bpList2 = {};
            this.activeCollisions = [];
            this.groundBodyIndex = [];
            for i = 1:this.bodyNum
                if(isinf(this.model.bodies{i}.Mp))
                    this.groundBodyIndex(end+1) = i;
                    this.model.bodies{i}.layer = 1;
                end
            end

			this.broadphase();
			this.narrowphase();
            this.constructCollisionOrder();
		end

        %%
        function constructCollisionOrder(this)
            for i = 1 : length(this.groundBodyIndex)
                groundIndex = this.groundBodyIndex(i);
                nextBodyQueue = groundIndex;
                while ~isempty(nextBodyQueue)
                    currentIndex = nextBodyQueue(1);
                    nextBodyQueue(1) = [];
                    for j = 1 : length(this.model.bodies{currentIndex}.neighbors)
                        neighborIndex = this.model.bodies{currentIndex}.neighbors(j);
                        if(this.model.bodies{neighborIndex}.layer > this.model.bodies{currentIndex}.layer + 1)
                            this.model.bodies{neighborIndex}.layer = this.model.bodies{currentIndex}.layer + 1;
                            nextBodyQueue(end + 1) = neighborIndex;
                        end
                    end
                end
            end

            %conLayers = arrayfun(@(conIndex) this.collisions{conIndex}.body1.layer + this.collisions{conIndex}.body2.layer, this.activeCollisions);
            conLayers = zeros(length(this.activeCollisions),1);
            for i = 1:length(this.activeCollisions)
                conIndex = this.activeCollisions(i);
                this.collisions{conIndex}.layer = this.collisions{conIndex}.body1.layer + this.collisions{conIndex}.body2.layer;
                conLayers(i) = this.collisions{conIndex}.layer;
            end

            for i = 1 : length(this.model.bodies)
                indices = find(ismember(this.activeCollisions, this.model.bodies{i}.collisions));
                [~, order] = sort(conLayers(indices));
                this.model.bodies{i}.collisions = this.model.bodies{i}.collisions(order);
            end

            %[~, idx] = sort(conLayers);
            layers = unique(conLayers);
            sortedCollisions = {};
            for i = 1:length(layers)
                sortedCollisions{end+1} = this.activeCollisions(conLayers==layers(i));
            end
            this.activeCollisions = sortedCollisions;
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
                if(this.collisions{body.index}.broken)
				    cdata = body.narrowphaseGround(Eg);
                    this.collisions{body.index}.setContacts(cdata);
                end
                this.collisions{body.index}.getConstraints();
                if(this.collisions{body.index}.contactNum ~= 0)
                    this.groundBodyIndex(end+1) = body.index;
                    body.layer = 1;
                    this.activeCollisions(end+1) = body.index;
                    if(this.model.useContactCaching)
                        this.collisions{body.index}.broken = false;
                    else
                        this.collisions{body.index}.broken = true;
                    end
                    body.collisions(end+1) = body.index;
                end
            end

			% Body-body collisions
			for i = 1 : length(this.bpList2)
				body1 = this.bpList2{i}{1};
				body2 = this.bpList2{i}{2};
                l = min([body1.index body2.index]) + 1;
                h = max([body1.index body2.index]);
                index = (2*this.bodyNum+2-l)*(l-1)/2 + h - l + 1;
                if(this.collisions{index}.broken)
				    cdata = body1.narrowphaseRigid(body2);
                    this.collisions{index}.setContacts(cdata);
                end
                this.collisions{index}.body1 = body1;
                this.collisions{index}.body2 = body2;
                this.collisions{index}.getConstraints();
                if(this.collisions{index}.contactNum ~= 0)
                    body1.neighbors(end+1) = body2.index;
                    body2.neighbors(end+1) = body1.index;
                    this.activeCollisions(end+1) = index;
                    if(this.model.useContactCaching)
                        this.collisions{index}.broken = false;
                    else
                        this.collisions{index}.broken = true;
                    end
                    body1.collisions(end+1) = index;
                    body2.collisions(end+1) = index;
                end
            end
        end

		%%
		function flag = isempty(this)
			flag = isempty(this.collisions);
		end
	end
end
