using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain_v2
    {
        int ID { get; set; }
        bool MaterialsModified { get; set; }

        Dictionary<int, IElement> ΙElementsDictionary { get; }
        Dictionary<int, Dictionary<DOFType, double>> Constraints { get; }
        IReadOnlyList<INode> Nodes { get; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        IDofOrdering DofOrdering { get; set; }
        Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary { get; }

        double[] Forces { get; }
        void ResetMaterialsModifiedProperty();
    }
}
