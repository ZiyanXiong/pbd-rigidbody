classdef Joint < handle
    properties
        body1   
        body2   % If this is a jont connecting to ground, body2 will be null
        constraintNum
        constraints
        lambdaLen
        ground  %If this is a ground joint
        index   % Global begining index for each collision 
        mIndces % Indices in the matrix
        J1I
        J2I     % If this is ground joint, J2 will be null
        b
        d
    end

    methods
        function this = Joint(body1, body2, ground)
            this.constraintNum = 0;
            this.lambdaLen = 0;
            this.constraints = {};
            this.body1 = body1;
            this.body2 = body2;
            this.ground = ground;
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