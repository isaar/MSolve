using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain_v2
    {
        int ID { get; }
        bool MaterialsModified { get; set; }

        Dictionary<int, IElement> ΙElementsDictionary { get; }
        Table<INode, DOFType, double> Constraints { get; }
        IReadOnlyList<INode> Nodes { get; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        IDofOrdering DofOrdering { get; set; }

        double[] Forces { get; } //TODO: this should be a Vector or IVector and stored elsewhere.
        void ResetMaterialsModifiedProperty();
    }
}
