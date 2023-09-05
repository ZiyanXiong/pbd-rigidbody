classdef ConAffine < apbd.ConBase
	%ConAffine Constraint for orthogonality

	properties
		affine
	end

	methods
		%%
		function this = ConAffine(affine)
			this = this@apbd.ConBase(6);
			this.affine = affine;
		end

		%%
		function init(this) %#ok<MANU>
			% Do nothing
		end

		%%
		function solve(this,h) %#ok<INUSD>
			af = this.affine;
			Wx = af.Wa(1);
			Wy = af.Wa(2);
			Wz = af.Wa(3);
			A11 = af.x(1);
			A12 = af.x(2);
			A13 = af.x(3);
			A21 = af.x(4);
			A22 = af.x(5);
			A23 = af.x(6);
			A31 = af.x(7);
			A32 = af.x(8);
			A33 = af.x(9);
			% XX: dC = [2*A11, 0, 0, 2*A21, 0, 0, 2*A31, 0, 0];
			i = 1;
			this.C(i) = A11^2 + A21^2 + A31^2 - 1;
			numerator = -this.C(i);
			t1 = Wx*2*A11;
			t2 = Wx*2*A21;
			t3 = Wx*2*A31;
			denominator = 2*(A11*t1 + A21*t2 + A31*t3);
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A11 = A11 + dlambda*t1;
			A21 = A21 + dlambda*t2;
			A31 = A31 + dlambda*t3;
			% YY: dC = [0, 2*A12, 0, 0, 2*A22, 0, 0, 2*A32, 0];
			i = 2;
			this.C(i) = A12^2 + A22^2 + A32^2 - 1;
			numerator = -this.C(i);
			t1 = Wy*2*A12;
			t2 = Wy*2*A22;
			t3 = Wy*2*A32;
			denominator = 2*(A12*t1 + A22*t2 + A32*t3);
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A12 = A12 + dlambda*t1;
			A22 = A22 + dlambda*t2;
			A32 = A32 + dlambda*t3;
			% ZZ: dC = [0, 0, 2*A13, 0, 0, 2*A23, 0, 0, 2*A33];
			i = 3;
			this.C(i) = A13^2 + A23^2 + A33^2 - 1;
			numerator = -this.C(i);
			t1 = Wz*2*A13;
			t2 = Wz*2*A23;
			t3 = Wz*2*A33;
			denominator = 2*(A13*t1 + A23*t2 + A33*t3);
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A13 = A13 + dlambda*t1;
			A23 = A23 + dlambda*t2;
			A33 = A33 + dlambda*t3;
			% XY: dC = [A12, A11, 0, A22, A21, 0, A32, A31, 0];
			i = 4;
			this.C(i) = A11*A12 + A21*A22 + A31*A32;
			numerator = -this.C(i);
			t1 = Wx*A12;
			t2 = Wy*A11;
			t3 = Wx*A22;
			t4 = Wy*A21;
			t5 = Wx*A32;
			t6 = Wy*A31;
			denominator = A12*t1 + A11*t2 + A22*t3 + A21*t4 + A32*t5 + A31*t6;
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A11 = A11 + dlambda*t1;
			A12 = A12 + dlambda*t2;
			A21 = A21 + dlambda*t3;
			A22 = A22 + dlambda*t4;
			A31 = A31 + dlambda*t5;
			A32 = A32 + dlambda*t6;
			% YZ: dC = [0, A13, A12, 0, A23, A22, 0, A33, A32];
			i = 5;
			this.C(i) = A12*A13 + A22*A23 + A32*A33;
			numerator = -this.C(i);
			t1 = Wy*A13;
			t2 = Wz*A12;
			t3 = Wy*A23;
			t4 = Wz*A22;
			t5 = Wy*A33;
			t6 = Wz*A32;
			denominator = A13*t1 + A12*t2 + A23*t3 + A22*t4 + A33*t5 + A32*t6;
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A12 = A12 + dlambda*t1;
			A13 = A13 + dlambda*t2;
			A22 = A22 + dlambda*t3;
			A23 = A23 + dlambda*t4;
			A32 = A32 + dlambda*t5;
			A33 = A33 + dlambda*t6;
			% ZX: dC = [A13, 0, A11, A23, 0, A21, A33, 0, A31];
			i = 6;
			this.C(i) = A11*A13 + A21*A23 + A31*A33;
			numerator = -this.C(i);
			t1 = Wx*A13;
			t2 = Wz*A11;
			t3 = Wx*A23;
			t4 = Wz*A21;
			t5 = Wx*A33;
			t6 = Wz*A31;
			denominator = A13*t1 + A11*t2 + A23*t3 + A21*t4 + A33*t5 + A31*t6;
			dlambda = numerator/denominator;
			this.lambda(i) = this.lambda(i) + dlambda;
			A11 = A11 + dlambda*t1;
			A13 = A13 + dlambda*t2;
			A21 = A21 + dlambda*t3;
			A23 = A23 + dlambda*t4;
			A31 = A31 + dlambda*t5;
			A33 = A33 + dlambda*t6;
			% Apply
			af.x(1) = A11;
			af.x(2) = A12;
			af.x(3) = A13;
			af.x(4) = A21;
			af.x(5) = A22;
			af.x(6) = A23;
			af.x(7) = A31;
			af.x(8) = A32;
			af.x(9) = A33;
		end

		%%
		function draw(this) %#ok<MANU>
			% Do nothing
		end
	end
end
