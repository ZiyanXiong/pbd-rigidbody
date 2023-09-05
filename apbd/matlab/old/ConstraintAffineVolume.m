classdef ConstraintAffineVolume < apbd.Constraint
    %ConstraintAffineVolume Constraint for volume preservation
    
    properties
		affine
    end
    
    methods
        %%
		function this = ConstraintAffineVolume(affine)
            this = this@apbd.Constraint(1,12);
			this.affine = affine;
        end

        %%
        function update(this,~)
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
            this.C = A11*A22*A33 - A11*A23*A32 - A12*A21*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31 - 1;
			%this.dC = [A22*A33 - A23*A32, A23*A31 - A21*A33, A21*A32 - A22*A31, A13*A32 - A12*A33, A11*A33 - A13*A31, A12*A31 - A11*A32, A12*A23 - A13*A22, A13*A21 - A11*A23, A11*A22 - A12*A21];
			t1 = A22*A33 - A23*A32;
			t2 = A23*A31 - A21*A33;
			t3 = A21*A32 - A22*A31;
			t4 = A13*A32 - A12*A33;
			t5 = A11*A33 - A13*A31;
			t6 = A12*A31 - A11*A32;
			t7 = A12*A23 - A13*A22;
			t8 = A13*A21 - A11*A23;
			t9 = A11*A22 - A12*A21;
			Wt1 = Wx*t1;
			Wt2 = Wy*t2;
			Wt3 = Wz*t3;
			Wt4 = Wx*t4;
			Wt5 = Wy*t5;
			Wt6 = Wz*t6;
			Wt7 = Wx*t7;
			Wt8 = Wy*t8;
			Wt9 = Wz*t9;
			this.WdC(1) = Wt1;
			this.WdC(2) = Wt2;
			this.WdC(3) = Wt3;
			this.WdC(4) = Wt4;
			this.WdC(5) = Wt5;
			this.WdC(6) = Wt6;
			this.WdC(7) = Wt7;
			this.WdC(8) = Wt8;
			this.WdC(9) = Wt9;
			this.dCWdC = t1*Wt1 + t2*Wt2 + t3*Wt3 + t4*Wt4 + t5*Wt5 + t6*Wt6 + t7*Wt7 + t8*Wt8 + t9*Wt9;
		end

		%%
		function apply(this,~,dlambda)
			this.affine.q(1:9) = this.affine.q(1:9) + dlambda*this.WdC(1:9)';
		end
    end
end
