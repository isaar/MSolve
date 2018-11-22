using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain_v2
    {
        Table<INode, DOFType, double> Constraints { get; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        ISubdomainFreeDofOrdering DofOrdering { get; set; }

        IReadOnlyList<IElement> Elements { get; } //TODO: perhaps this should be a set

        Vector Forces { get; } //TODO: this should be a Vector or IVector and stored elsewhere.

        int ID { get; }

        bool MaterialsModified { get; set; }

        IReadOnlyList<INode> Nodes { get; } //TODO: perhaps this should be a set

        void ResetMaterialsModifiedProperty();
    }
}
