classdef Contact < handle
    properties
        nw
        x1
        x2
    end

    methods
        function this = Contact()
			this.nw = zeros(3,1);
			this.x1 = zeros(3,1);
			this.x2 = zeros(3,1);
        end
        function setData(this, cdata)
			this.nw = cdata.nw;
			this.x1 = cdata.x1;
			this.x2 = cdata.x2;
        end
    end

end