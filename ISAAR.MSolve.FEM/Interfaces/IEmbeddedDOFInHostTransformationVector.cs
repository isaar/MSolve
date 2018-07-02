using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IEmbeddedDOFInHostTransformationVector
    {
        IList<DOFType> GetDependentDOFTypes { get; }
        IList<IList<DOFType>> GetDOFTypesOfHost(EmbeddedNode node);
        double[][] GetTransformationVector(EmbeddedNode node);
    }
}
