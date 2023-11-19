classdef Shape < handle
	%%
	properties
	end

	%%
	methods
		%%
		function this = Shape()
		end
	end
	
	%%
	methods (Abstract)
		%%
		I = computeInertia(this)

		%%
		flag = broadphaseGround(this,E,Eg)

		%%
		cdata = narrowphaseGround(this,E,Eg)

		%%
		flag = broadphaseShape(this,Ethis,that,Ethat)

		%%
		cdata = narrowphaseShape(this,Ethis,that,Ethat)

		%%
		[F,V] = draw(this,E,color,axisSize)
	end
end
