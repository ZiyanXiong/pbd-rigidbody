classdef ConstraintAffineOrtho < apbd.Constraint
	%ConstraintAffineOrtho Constraint for orthogonality

	properties
		affine
	end

	methods
		%%
		function this = ConstraintAffineOrtho(affine)
			this = this@apbd.Constraint(6,12);
			this.affine = affine;
		end

		%%
		function update(this,i)
			af = this.affine;
			Wx = af.W(1);
			Wy = af.W(2);
			Wz = af.W(3);
			A11 = af.q(1);
			A12 = af.q(2);
			A13 = af.q(3);
			A21 = af.q(4);
			A22 = af.q(5);
			A23 = af.q(6);
			A31 = af.q(7);
			A32 = af.q(8);
			A33 = af.q(9);
			switch i
				case 1
					this.C(1) = A11^2 + A21^2 + A31^2 - 1;
					%this.dC(1,:) = [2*A11, 0, 0, 2*A21, 0, 0, 2*A31, 0, 0];
					t1 = Wx*2*A11;
					t2 = Wx*2*A21;
					t3 = Wx*2*A31;
					this.WdC(1,1) = t1;
					this.WdC(1,4) = t2;
					this.WdC(1,7) = t3;
					this.dCWdC(1) = 2*(A11*t1 + A21*t2 + A31*t3);
				case 2
					this.C(2) = A12^2 + A22^2 + A32^2 - 1;
					%this.dC(2,:) = [0, 2*A12, 0, 0, 2*A22, 0, 0, 2*A32, 0];
					t1 = Wy*2*A12;
					t2 = Wy*2*A22;
					t3 = Wy*2*A32;
					this.WdC(2,2) = t1;
					this.WdC(2,5) = t2;
					this.WdC(2,8) = t3;
					this.dCWdC(2) = 2*(A12*t1 + A22*t2 + A32*t3);
				case 3
					this.C(3) = A13^2 + A23^2 + A33^2 - 1;
					%this.dC(3,:) = [0, 0, 2*A13, 0, 0, 2*A23, 0, 0, 2*A33];
					t1 = Wz*2*A13;
					t2 = Wz*2*A23;
					t3 = Wz*2*A33;
					this.WdC(3,3) = t1;
					this.WdC(3,6) = t2;
					this.WdC(3,9) = t3;
					this.dCWdC(3) = 2*(A13*t1 + A23*t2 + A33*t3);
				case 4
					this.C(4) = A11*A12 + A21*A22 + A31*A32;
					%this.dC(4,:) = [A12, A11, 0, A22, A21, 0, A32, A31, 0];
					t1 = Wx*A12;
					t2 = Wy*A11;
					t3 = Wx*A22;
					t4 = Wy*A21;
					t5 = Wx*A32;
					t6 = Wy*A31;
					this.WdC(4,1) = t1;
					this.WdC(4,2) = t2;
					this.WdC(4,4) = t3;
					this.WdC(4,5) = t4;
					this.WdC(4,7) = t5;
					this.WdC(4,8) = t6;
					this.dCWdC(4) = A12*t1 + A11*t2 + A22*t3 + A21*t4 + A32*t5 + A31*t6;
				case 5
					this.C(5) = A12*A13 + A22*A23 + A32*A33;
					%this.dC(5,:) = [0, A13, A12, 0, A23, A22, 0, A33, A32];
					t1 = Wy*A13;
					t2 = Wz*A12;
					t3 = Wy*A23;
					t4 = Wz*A22;
					t5 = Wy*A33;
					t6 = Wz*A32;
					this.WdC(5,2) = t1;
					this.WdC(5,3) = t2;
					this.WdC(5,5) = t3;
					this.WdC(5,6) = t4;
					this.WdC(5,8) = t5;
					this.WdC(5,9) = t6;
					this.dCWdC(5) = A13*t1 + A12*t2 + A23*t3 + A22*t4 + A33*t5 + A32*t6;
				case 6
					this.C(6) = A11*A13 + A21*A23 + A31*A33;
					%this.dC(6,:) = [A13, 0, A11, A23, 0, A21, A33, 0, A31];
					t1 = Wx*A13;
					t2 = Wz*A11;
					t3 = Wx*A23;
					t4 = Wz*A21;
					t5 = Wx*A33;
					t6 = Wz*A31;
					this.WdC(6,1) = t1;
					this.WdC(6,3) = t2;
					this.WdC(6,4) = t3;
					this.WdC(6,6) = t4;
					this.WdC(6,7) = t5;
					this.WdC(6,9) = t6;
					this.dCWdC(6) = A13*t1 + A11*t2 + A23*t3 + A21*t4 + A33*t5 + A31*t6;
			end
		end

		%%
		function apply(this,i,dlambda)
			switch i
				case 1
					j = [1 4 7];
				case 2
					j = [2 5 8];
				case 3
					j = [3 6 9];
				case 4
					j = [1 2 4 5 7 8];
				case 5
					j = [2 3 5 6 8 9];
				case 6
					j = [1 3 4 6 7 9];
			end
			this.affine.q(j) = this.affine.q(j) + dlambda*this.WdC(i,j)';
		end
	end
end
