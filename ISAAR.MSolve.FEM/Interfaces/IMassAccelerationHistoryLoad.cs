using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IMassAccelerationHistoryLoad
    {
        DOFType DOF { get; }
        double this[int currentTimeStep] { get; }
    }
}
