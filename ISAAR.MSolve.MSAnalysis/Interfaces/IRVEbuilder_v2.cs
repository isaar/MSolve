using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates the nesessary methods that should be implemented by builders of models intended to be used as rves in Multiscale Problems
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IRVEbuilder_v2
    {
        Tuple<Model_v2, Dictionary<int, Node_v2>,double> GetModelAndBoundaryNodes();
        IRVEbuilder_v2 Clone(int a);
    }
}
