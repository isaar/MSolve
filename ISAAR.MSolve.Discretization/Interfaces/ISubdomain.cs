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

		Dictionary<int, IElement> ΙElementsDictionary { get; }
		Dictionary<int, Dictionary<DOFType, double>> Constraints { get; }
        IReadOnlyList<INode> Nodes { get; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        IDofOrdering DofOrdering { get; set; }

        Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; }
        Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary { get; }
		
		double[] Forces { get; }
		void ResetMaterialsModifiedProperty();

	}
}
