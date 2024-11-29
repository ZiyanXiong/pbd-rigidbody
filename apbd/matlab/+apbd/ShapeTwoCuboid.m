classdef ShapeTwoCuboid < apbd.Shape
	%%
	properties
		cuboid1
        cuboid2
        E1
        E2
	end
	
	%%
	methods
		%%
        function this = ShapeTwoCuboid(sides1, sides2, halfDis, halfAngle)
			this = this@apbd.Shape();
			this.cuboid1 = apbd.ShapeCuboid(sides1);
            this.cuboid2 = apbd.ShapeCuboid(sides2);
            E = eye(4);
    		R = se3.aaToMat([0 1 0], halfAngle);
            E(1:3,1:3) = R;
			E(1:3,4) = [-halfDis 0 0]';
            this.E1 = E;

            E(1:3,1:3) = R';
			E(1:3,4) = [halfDis 0 0]';
            this.E2 = E;
		end
		
		%%
		function I = computeInertia(this,density)
            % I is the 6x1 diagonal rigid inertia, assuming that the frame origin is at the
            % center of mass and the axes are oriented along the principal axes. We store
            % the rotation on top of translations, so that I(1:3) is the rotational inertia
            % and I(4:6) is the translational inertia.
            I1 = se3.inertiaCuboid(this.cuboid1.sides,density);
            I2 = se3.inertiaCuboid(this.cuboid2.sides,density);
            A1 = se3.Ad(se3.inv(this.E1));
            A2 = se3.Ad(se3.inv(this.E2));
            I = A1' * diag(I1) * A1 + A2' * diag(I2) * A2;
            I = diag(I);
		end

        %%
        function xlc = toCenterLocal(this, E,xl)
            xlc = E * [xl;1];
            xlc = xlc(1:3);
        end

		%%
		function flag = broadphaseGround(this,E,Eg)
            flag1 = this.cuboid1.broadphaseGround(E*this.E1,Eg);
            flag2 = this.cuboid2.broadphaseGround(E*this.E2,Eg);
            flag = flag1 + flag2;
		end

		%%
		function cdata = narrowphaseGround(this,E,Eg)
            cdata1 = this.cuboid1.narrowphaseGround(E*this.E1,Eg);
            cdata2 = this.cuboid1.narrowphaseGround(E*this.E2,Eg);
            for i = 1: length(cdata1)
                cdata1(i).x1 = this.toCenterLocal(this.E1,cdata1(i).x1);
            end
            for i = 1: length(cdata2)
                cdata2(i).x1 = this.toCenterLocal(this.E2,cdata2(i).x1);
            end
            cdata = [cdata1 cdata2];
		end

		%%
		function flag = broadphaseShape(this,E1,that,E2)
			if isa(that,'apbd.ShapeCuboid')
                flag1 = this.cuboid1.broadphaseShape(E1*this.E1,that,E2);
                flag2 = this.cuboid2.broadphaseShape(E1*this.E2,that,E2);
                flag = flag1 + flag2;
            elseif isa(that, 'apbd.ShapeTwoCuboid')
                flag1 = this.cuboid1.broadphaseShape(E1*this.E1,that.cuboid1,E2*that.E1);
                flag2 = this.cuboid1.broadphaseShape(E1*this.E1,that.cuboid2,E2*that.E2);
                flag3 = this.cuboid2.broadphaseShape(E1*this.E2,that.cuboid1,E2*that.E1);
                flag4 = this.cuboid2.broadphaseShape(E1*this.E2,that.cuboid2,E2*that.E2);
                flag = flag1 + flag2 + flag3 + flag4;
            else
				error('Unsupported shape');
			end
		end

		%%
		function cdata = narrowphaseShape(this,E1,that,E2)
			if isa(that,'apbd.ShapeCuboid')
                cdata1 = this.cuboid1.narrowphaseShape(E1*this.E1,that,E2);
                cdata2 = this.cuboid2.narrowphaseShape(E1*this.E2,that,E2);
                for i = 1: length(cdata1)
                    cdata1(i).x1 = this.toCenterLocal(this.E1,cdata1(i).x1);
                end
                for i = 1: length(cdata2)
                    cdata2(i).x1 = this.toCenterLocal(this.E2,cdata2(i).x1);
                end
                cdata = [cdata1 cdata2];
            elseif isa(that, 'apbd.ShapeTwoCuboid')
                cdata1 = this.cuboid1.narrowphaseShape(E1*this.E1,that.cuboid1,E2*that.E1);
                cdata2 = this.cuboid1.narrowphaseShape(E1*this.E1,that.cuboid2,E2*that.E2);
                cdata3 = this.cuboid2.narrowphaseShape(E1*this.E2,that.cuboid1,E2*that.E1);
                cdata4 = this.cuboid2.narrowphaseShape(E1*this.E2,that.cuboid2,E2*that.E2);
                for i = 1: length(cdata1)
                    cdata1(i).x1 = this.toCenterLocal(this.E1,cdata1(i).x1);
                    cdata1(i).x2 = that.toCenterLocal(that.E1,cdata1(i).x2);
                end
                for i = 1: length(cdata2)
                    cdata2(i).x1 = this.toCenterLocal(this.E1,cdata2(i).x1);
                    cdata2(i).x2 = that.toCenterLocal(that.E2,cdata2(i).x2);
                end
                for i = 1: length(cdata3)
                    cdata3(i).x1 = this.toCenterLocal(this.E2,cdata3(i).x1);
                    cdata3(i).x2 = that.toCenterLocal(that.E1,cdata3(i).x2);
                end
                for i = 1: length(cdata4)
                    cdata4(i).x1 = this.toCenterLocal(this.E2,cdata4(i).x1);
                    cdata4(i).x2 = that.toCenterLocal(that.E2,cdata4(i).x2);
                end
                cdata = [cdata1 cdata2 cdata3 cdata4];
            else
				error('Unsupported shape');
            end
		end

		%%
		function [F,V] = draw(this,E,color,axisSize)
			if nargin < 4
				axisSize = 0;
			end
			if nargin < 3
				color = [0.5 0.5 0.5];
            end
            
			[F,V] = se3.patchCuboid(E*this.E1,this.cuboid1.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',color);
			[F,V] = se3.patchCuboid(E*this.E2,this.cuboid2.sides);
			patch('Faces',F,'Vertices',V,'FaceColor',color);
			if axisSize > 0
				se3.drawAxis(E,axisSize);
			end
		end
	end

end
