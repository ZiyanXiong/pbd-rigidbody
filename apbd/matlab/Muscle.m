classdef Muscle < handle
    properties
        bodies
        viaPointNum
        constraintNum
        constraints
        compliance
        lambdaLen
        ground  %If this is a ground joint
        index   % Global begining index for each collision 
        mIndces % Indices in the matrix
        JIs
        b
        d
    end

    methods
        function this = Muscle(bodies, ground)
            this.constraintNum = 0;
            this.lambdaLen = 0;
            this.constraints = {};
            this.bodies = bodies;
            this.ground = ground;
            this.compliance = 0;
            this.index = 0;
        end

        function draw(this)
            for i = 1 : this.constraintNum
                this.constraints{i}.draw();
            end
        end
    end

	methods (Abstract)

		%% Init Joint
		init(this,h,hs,timestep);

        %% Compute J and b
        computeJ_b(this);

        %% Compute b
        compute_b;

        %% Apply force and torque
        applyForceTorque(this,f,t,timestep);

        %% Apply lambda
        applyLambdas(this,lambdas);

    end
end