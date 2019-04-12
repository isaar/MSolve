using ISAAR.MSolve.Discretization.FreedomDegrees;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public interface ITimeDependentNodalLoad
    {
        Node Node { get; set; }
        IDofType DOF { get; set; }

        double GetLoadAmount(int timeStep);
    }
}
