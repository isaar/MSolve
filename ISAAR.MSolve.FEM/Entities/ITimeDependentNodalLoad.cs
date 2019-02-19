using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public interface ITimeDependentNodalLoad
    {
        Node_v2 Node { get; set; }
        DOFType DOF { get; set; }

        double GetLoadAmount(int timeStep);
    }
}
