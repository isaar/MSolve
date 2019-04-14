using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IMassAccelerationHistoryLoad
    {
        IDofType DOF { get; }
        double this[int currentTimeStep] { get; }
    }
}
