classdef ConColl < apbd.ConBase
	%ConColl Collision constraint

	properties
		nw % collision normal wrt world (3x1)
		d  % penetration depth (negative if collision)
		
		dlambdaNor % Lagrange multiplier update for normal force
        dlambdaSP % Lagrange multiplier saved for shock porpagation
		s % scale for display
	end

	methods
		%%
		function this = ConColl()
			this = this@apbd.ConBase(3);
			this.nw = zeros(3,1);
			this.d = 0;
			this.dlambdaNor = zeros(3,1);
            this.dlambdaSP = zeros(3,1);
			this.s = 1;
		end

		%%
		function solve(this,h) %#ok<INUSD>
			% For this class, call solveNormal() and solveTangent() instead
		end
	end

	%%
	methods (Abstract)

		%%
		solveNorPos(this)

		%%
		solveTanVel(this,k,ks,hs)
	end

	%%
	methods (Static)
		%%
		function [tanx,tany] = generateTangents(nor)
			if abs(nor(3)) < 1e-6
				tmp = [0 0 1]';
			else
				tmp = [1 0 0]';
			end
			tany = se3.cross(nor, tmp);
			tany = tany/norm(tany);
			tanx = se3.cross(tany,nor);
			% Debug to make sure that rotating the tangents won't affect
			% the results.
			% R = se3.aaToMat(nor,pi/4);
			% tanx = R*tanx;
			% tany = R*tany;
        end
	end
end
