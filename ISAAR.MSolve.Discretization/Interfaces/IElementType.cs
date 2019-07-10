using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh;
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

    public interface IElementType
    {
        CellType CellType { get; }
        IElementDofEnumerator DofEnumerator { get; set; }
        IMatrix StiffnessMatrix(IElement element);
        IMatrix MassMatrix(IElement element);
        IMatrix DampingMatrix(IElement element);
        IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element);
    }
}
