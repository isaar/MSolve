using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public enum ElementDimensions
    {
        Unknown = 0,
        OneD = 1,
        TwoD = 2,
        ThreeD = 3
    }

    public interface IElementType_v2
    {
        IElementDofEnumerator_v2 DofEnumerator { get; set; }
        IMatrix StiffnessMatrix(IElement_v2 element);
        IMatrix MassMatrix(IElement_v2 element);
        IMatrix DampingMatrix(IElement_v2 element);
        IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element);
    }
}
