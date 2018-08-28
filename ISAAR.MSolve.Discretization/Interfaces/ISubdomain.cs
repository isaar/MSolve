using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain
    {
		int ID { get; set; }
		int TotalDOFs { get; set; }
		bool MaterialsModified { get; set; }

		Dictionary<int, IElement> ΙElementsDictionary { get; }
		Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; }
		Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary { get; }
		
		double[] Forces { get; }
		void ResetMaterialsModifiedProperty();

	}
}
