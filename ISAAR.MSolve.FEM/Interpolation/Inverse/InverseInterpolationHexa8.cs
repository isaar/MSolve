using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Interpolation.Inverse
{
	/// <summary>
	/// Inverse mapping of the isoparametric interpolation of a hexahedral finite element with 8 nodes. Since the original
	/// mapping is linear, there are analytic formulas, which are presented in 
	/// "The inverse mapping and distortion measures for 8-node hexahedral isoparametric elements", K. -Y. YuanY. -S. HuangH. -T. YangT. H. H. Pian, 1994
	/// https://link.springer.com/article/10.1007/BF00350284
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class InverseInterpolationHexa8:IInverseInterpolation3D
    {
	    public InverseInterpolationHexa8(IReadOnlyList<Node3D> nodes)
	    {
		    
	    }

		public NaturalPoint3D TransformPointCartesianToNatural(CartesianPoint3D point)
		{
			throw new NotImplementedException();
		}
	}
}
