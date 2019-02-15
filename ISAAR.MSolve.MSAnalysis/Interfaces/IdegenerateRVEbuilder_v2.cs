using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates additional methods that should be implemented by rveBuilders that will be used in FE2 3D to 2D degenerate analysis
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IdegenerateRVEbuilder_v2 : IRVEbuilder_v2
    {
        Dictionary<Node, IList<DOFType>> GetModelRigidBodyNodeConstraints(Model_v2 model);
    }
}
