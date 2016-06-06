using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Interfaces
{
    public interface IMassAccelerationHistoryLoad
    {
        DOFType DOF { get; }
        double this[int currentTimeStep] { get; }
    }
}
