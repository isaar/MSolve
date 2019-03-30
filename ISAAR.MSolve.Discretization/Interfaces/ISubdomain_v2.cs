using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: tidy up the methods that concern material state
namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface ISubdomain_v2
    {
        Table<INode, DOFType, double> Constraints { get; }

        ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }

        /// <summary>
        /// This should be set when the analyzer decides. E.g. an XFEM or adaptive FEM analyzer would need to create a new dof 
        /// ordering, whenever the crack propagates or the mesh is refined respectively.
        /// </summary>
        ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

        IReadOnlyList<IElement_v2> Elements { get; } //TODO: perhaps this should be a set

        Vector Forces { get; set; } //TODO: this should be a Vector or IVector and stored elsewhere.

        int ID { get; }

        bool MaterialsModified { get; set; }

        IReadOnlyList<INode> Nodes { get; } //TODO: perhaps this should be a set

        double[] CalculateElementIncrementalConstraintDisplacements(IElement_v2 element, double constraintScalingFactor); //TODO: this should be done by a dedicated class instead of the subdomain

        void ClearMaterialStresses();

        IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution); //TODO: this should be done by a dedicated class instead of the subdomain

        void ResetMaterialsModifiedProperty();

        void ScaleConstraints(double scalingFactor); //TODO: this should be done by a dedicated class instead of the subdomain

        void SaveMaterialState();
    }
}
