using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain
    {
		int ID { get; set; }
		int TotalDOFs { get; set; } //TODO: rename it to NumFreeDofs or even better expose a DofOrderer 
        bool MaterialsModified { get; set; }

		Dictionary<int, IElement> ElementsDictionary { get; }
        Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; }
        Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary { get; }
		
		double[] Forces { get; }
		void ResetMaterialsModifiedProperty();

	}
}
