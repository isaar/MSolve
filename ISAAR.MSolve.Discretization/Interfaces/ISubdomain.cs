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

		Dictionary<int, IElement> ElementsDictionary { get; }
		Dictionary<int, INode> NodesDictionary { get; }
		Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary { get; }
		Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary { get; }

		IList<INode> Nodes { get; }
		double[] Forces { get; }

		void EnumerateDOFs();
		void AssignGlobalNodalDOFsFromModel();
		void BuildNodesDictionary();
	    double[] GetLocalVectorFromGlobal(IElement element, IVector globalVector);
    }
}
