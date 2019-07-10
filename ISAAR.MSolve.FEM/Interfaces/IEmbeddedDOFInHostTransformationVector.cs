using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Embedding;


namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedDOFInHostTransformationVector
    {
        IList<IDofType> GetDependentDOFTypes { get; }
        IReadOnlyList<IReadOnlyList<IDofType>> GetDOFTypesOfHost(EmbeddedNode node);
        double[][] GetTransformationVector(EmbeddedNode node);
    }
}
